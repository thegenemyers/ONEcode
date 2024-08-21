/*  File: ONElogan.c
 *  Authors: Rayan Chikhi (rayan.chikhi@pasteur.fr)
 *           Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: ONEcode file format conversion for Logan contigs 
 *              (https://github.com/IndexThePlanet/Logan)
 *              Record sequences, abundances, graph links
 * Created: August 15 2024
 *-------------------------------------------------------------------
 */

#include "ONElib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <omp.h>

#define MAX_LINE_LENGTH 1000000
#define MAX_HEADER_LENGTH 1000
#define CHUNK_SIZE 1000  // Number of sequences to process in each chunk

void strip_newline(char *s) {
    char *p = strchr(s, '\n');
    if (p) *p = '\0';
}

typedef struct {
    char *header;
    char *sequence;
    size_t seq_length;
    float coverage;
    char **links;
    int link_count;
} FastaEntry;

void process_chunk(OneFile *vf, FastaEntry *entries, int count) {
    for (int i = 0; i < count; i++) {
        FastaEntry *entry = &entries[i];

        // Write sequence
        oneInt(vf, 0) = strlen(entry->header);
        oneWriteLine(vf, 'S', entry->seq_length, entry->sequence);

        // Write coverage
        oneReal(vf, 0) = entry->coverage;
        oneWriteLine(vf, 'K', 0, NULL);

        // Write links
        for (int j = 0; j < entry->link_count; j++) {
            char *link = entry->links[j];
            char orientation = link[2];
            int64_t id = atoll(link + 4);
            char end = link[strlen(link) - 1];
            oneChar(vf, 0) = orientation;
            oneInt(vf, 1) = id;
            oneChar(vf, 2) = end;
            oneWriteLine(vf, 'L', 0, NULL);
        }

        // Free memory
        free(entry->header);
        free(entry->sequence);
        for (int j = 0; j < entry->link_count; j++) {
            free(entry->links[j]);
        }
        free(entry->links);
    }
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_fasta> <output_one> <num_threads>\n", argv[0]);
        exit(1);
    }

    int num_threads = atoi(argv[3]);
    if (num_threads < 1) {
        fprintf(stderr, "Error: Number of threads must be at least 1\n");
        exit(1);
    }

    FILE *input = fopen(argv[1], "r");
    if (!input) {
        fprintf(stderr, "Error: Cannot open input file %s\n", argv[1]);
        exit(1);
    }

    OneSchema *schema = oneSchemaCreateFromText(
        "1 3 def 2 1\n"
        "P 3 seq\n"
        "S 5 logan\n"
        "O S 1 3 DNA\n"
        "D K 1 4 REAL\n"
        "D L 3 4 CHAR 3 INT 4 CHAR\n"
    );

    OneFile *vf = oneFileOpenWriteNew(argv[2], schema, "logan", true, num_threads);
    if (!vf) {
        fprintf(stderr, "Error: Cannot open output file %s\n", argv[2]);
        fclose(input);
        exit(1);
    }

    oneAddProvenance(vf, "ONElogan", "1.2", "ONElogan input_fasta output_one num_threads");

    char line[MAX_LINE_LENGTH];
    FastaEntry *entries = malloc(CHUNK_SIZE * sizeof(FastaEntry));
    int entry_count = 0;
    int seq_count = 0;
    size_t max_seq_length = 0;
    size_t total_seq_length = 0;

    while (fgets(line, sizeof(line), input)) {
        strip_newline(line);
        if (line[0] == '>' || feof(input)) {
            if (entry_count > 0 && (line[0] == '>' || feof(input))) {
                FastaEntry *entry = &entries[entry_count - 1];
                seq_count++;
                total_seq_length += entry->seq_length;
                if (entry->seq_length > max_seq_length) max_seq_length = entry->seq_length;
            }

            if (entry_count == CHUNK_SIZE || (feof(input) && entry_count > 0)) {
                #pragma omp parallel num_threads(num_threads)
                {
                    int thread_id = omp_get_thread_num();
                    int chunk_size = (entry_count + num_threads - 1) / num_threads;
                    int start = thread_id * chunk_size;
                    int end = (thread_id == num_threads - 1) ? entry_count : (start + chunk_size);
                    process_chunk(&vf[thread_id], &entries[start], end - start);
                }
                entry_count = 0;
            }

            if (!feof(input)) {
                FastaEntry *entry = &entries[entry_count++];
                entry->header = strdup(line + 1);
                entry->sequence = NULL;
                entry->seq_length = 0;
                entry->coverage = 0;
                entry->links = NULL;
                entry->link_count = 0;

                char *token = strtok(entry->header, " \t");
                while (token != NULL) {
                    if (strncmp(token, "ka:f:", 5) == 0) {
                        entry->coverage = atof(token + 5);
                    } else if (strncmp(token, "L:", 2) == 0) {
                        entry->links = realloc(entry->links, (entry->link_count + 1) * sizeof(char*));
                        entry->links[entry->link_count++] = strdup(token);
                    }
                    token = strtok(NULL, " \t");
                }
            }
        } else {
            FastaEntry *entry = &entries[entry_count - 1];
            size_t line_len = strlen(line);
            entry->sequence = realloc(entry->sequence, entry->seq_length + line_len + 1);
            for (size_t i = 0; i < line_len; i++) {
                if (!isspace(line[i])) {
                    entry->sequence[entry->seq_length++] = tolower(line[i]);
                }
            }
            entry->sequence[entry->seq_length] = '\0';
        }
    }

    // Process any remaining entries
    if (entry_count > 0) {
        #pragma omp parallel num_threads(num_threads)
        {
            int thread_id = omp_get_thread_num();
            int chunk_size = (entry_count + num_threads - 1) / num_threads;
            int start = thread_id * chunk_size;
            int end = (thread_id == num_threads - 1) ? entry_count : (start + chunk_size);
            process_chunk(&vf[thread_id], &entries[start], end - start);
        }
    }

    // Update header information
    vf->info['S']->given.count = seq_count;
    vf->info['K']->given.count = seq_count;
    vf->info['L']->given.count = seq_count;
    // Note: We don't know the exact count of 'L' lines, so we don't update it here

    oneFileClose(vf);
    oneSchemaDestroy(schema);
    fclose(input);
    free(entries);

    printf("Conversion complete. Wrote %d sequences.\n", seq_count);
    return 0;
}
