# Python interface for ONEcode

## Installation and compilation

To set up the Python bindings for the C++ interface, first install [pybind11](https://github.com/pybind/pybind11). The following was done using the [conda installion](https://pybind11.readthedocs.io/en/latest/installing.html#include-with-conda-forge).

After installing, first compile the library:
```
make
```
Then compile the Python binding with:
```
g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) pyONElib.cpp -o ONEcode$(python3-config --extension-suffix) ONElib.o
```
This should result in a `.so` file.

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

## Further work
* Use `setuptools` to compile and install, so the library will be available system-wide and publishable to `pip`.