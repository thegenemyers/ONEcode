# Python interface for ONEcode

## Installation and compilation

The Python bindings for the C++ interface use [pybind11](https://github.com/pybind/pybind11), which is included as a git submodule in `extern/pybind11`.  Fetch it with
```
git submodule update --init extern/pybind11
```
(or install pybind11 some other way, e.g. with [conda](https://pybind11.readthedocs.io/en/latest/installing.html#include-with-conda-forge), and use `$(python3 -m pybind11 --includes)` in place of `-Iextern/pybind11/include` below).

After installing, first compile the library:
```
make
```
Then compile the Python binding with:
```
g++ -O3 -Wall -shared -std=c++11 -fPIC -Iextern/pybind11/include $(python3-config --includes) pyONElib.cpp -o ONEcode$(python3-config --extension-suffix) ONElib.o
```
On macOS add `-undefined dynamic_lookup`, so that the Python symbols are resolved when the module is loaded.  This should result in a `.so` file, which Python finds if it is in the working directory or on `PYTHONPATH`.

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

`readLine()` returns the type of the line it read as a one-character string, or `""` at the end of the file.  `count(t)`, `max(t)` and `total(t)` give the number of lines of type `t`, their maximum list length and their total list length: as given in the header of a file being read, and as accumulated so far for a file being written.

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

print("Opened 1seq one with", onefile.count('S'), "sequences")

while onefile.readLine():
    if onefile.lineType() == 'S':
        print("Sequence length", onefile.length())
```

## Further work
* Use `setuptools` to compile and install, so the library will be available system-wide and publishable to `pip`.