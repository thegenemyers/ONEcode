# Python interface for ONEcode

## Installation and compilation

The ONEcode Python bindings provide a complete Python interface to the ONEcode library using pybind11.

```bash
# First compile the ONElib library
make

# Then compile the Python binding
make python
```

This creates a `.so` file.

## Usage

Run `python` in the working directory:
```Python
import ONEcode

# Read schema from a file
schema = ONEcode.ONEschema(open("...").read())

# Construct a ONEfile object
onefile = ONEcode.ONEfile("...", "r", schema, "...", 1)

# Various methods
onefile.readLine()
```

### Test:
```Python
import ONEcode

schema = ONEcode.ONEschema(
  "P 3 seq                 SEQUENCE\n"
  "S 6 segseq              segment sequences - objects are 1:1 with those in seg file\n"
  "S 7 readseq             read sequences\n"
  "O S 1 3 DNA             sequence: the DNA string\n"
  "D I 1 6 STRING          id - sequence identifier; unnecessary for segments\n"
)

onefile = ONEcode.ONEfile("./TEST/small.seq", "r", schema, "", 1)

print("Opened 1seq one with", onefile.givenCount('S'), "sequences")

while onefile.readLine():
    if onefile.lineType() == 'S':
        print("Sequence length", onefile.length())
```

## API Reference

### ONEschema Class

- `ONEschema(schema_text: str)`: Create a schema from schema text

### ONEfile Class

#### Constructors
- `ONEfile(path: str)`: Open existing file for reading (1 thread)
- `ONEfile(path: str, nthreads: int)`: Open existing file with multiple threads
- `ONEfile(path: str, mode: str, schema: ONEschema, type: str, nthreads: int)`: Full constructor
- `ONEfile(path: str, mode: str, source: ONEfile, nthreads: int)`: Create from another file

#### Reading Methods
- `readLine() -> char`: Read next line, returns line type
- `lineType() -> char`: Get current line type
- `lineNumber() -> int`: Get current line number
- `length() -> int`: Get length of current list
- `getInt(field: int) -> int`: Get integer field value
- `getReal(field: int) -> float`: Get real field value
- `getChar(field: int) -> char`: Get character field value
- `getIntList() -> list[int]`: Get integer list
- `getRealList() -> list[float]`: Get real list
- `getString() -> str`: Get string from list
- `getStringList() -> list[str]`: Get string list
- `getDNAchar() -> str`: Get DNA as character string
- `getDNA2bit() -> bytes`: Get DNA in 2-bit compressed format
- `getComment() -> str`: Read comment line

#### Writing Methods
- `setInt(field: int, value: int)`: Set integer field
- `setReal(field: int, value: float)`: Set real field
- `setChar(field: int, value: char)`: Set character field
- `writeLine(lineType: str)`: Write line with no list
- `writeLine(lineType: str, s: str)`: Write line with string
- `writeLine(lineType: str, slist: list[str])`: Write line with string list
- `writeLineIntList(lineType: str, data: list[int])`: Write line with integer list
- `writeLineRealList(lineType: str, data: list[float])`: Write line with real list
- `writeLineDNA2bit(lineType: str, dnaLen: int, dnaBuf: bytes)`: Write DNA in 2-bit format
- `writeComment(comment: str)`: Write comment line

#### Metadata Methods
- `givenCount(lineType: str) -> int`: Get given count for line type
- `givenMax(lineType: str) -> int`: Get given max for line type
- `givenTotal(lineType: str) -> int`: Get given total for line type
- `currentCount(lineType: str) -> int`: Get current count for line type
- `currentMax(lineType: str) -> int`: Get current max for line type
- `currentTotal(lineType: str) -> int`: Get current total for line type
- `count(lineType: str) -> int`: Get count (write mode uses accum, read mode uses given)
- `max(lineType: str) -> int`: Get max value
- `total(lineType: str) -> int`: Get total value
- `fileName() -> str`: Get file name

#### Provenance Methods
- `addProvenance(prog: str, version: str, commandLine: str) -> bool`: Add provenance
- `inheritProvenance(source: ONEfile) -> bool`: Inherit provenance from another file
- `addReference(filename: str, count: int) -> bool`: Add reference to another file
- `inheritReference(source: ONEfile) -> bool`: Inherit reference from another file

#### Navigation Methods
- `gotoObject(lineType: str, index: int) -> bool`: Jump to specific object

#### Schema Methods
- `checkSchemaText(text: str) -> bool`: Check if file matches schema text
- `checkSchema(schema: ONEschema, isRequired: bool) -> bool`: Check against schema object

## Notes

- All methods that take a `lineType` parameter expect a single character string (e.g., `'S'`, `'I'`)
- The module handles conversion between Python types and C++ types automatically
- Files are automatically closed when the ONEfile object is destroyed
- For writing files, mode `"w"` creates ASCII format, `"wb"` creates binary format
- Thread count parameter allows parallel decompression when reading binary files