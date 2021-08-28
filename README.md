# Statislib

A Fortran library for statistical distributions. The library includes various distributions' random number generators, probability density function or probability mass function, and cumulative distribution function. Detail usages for each distribution are shown in corresponding documents.

## Getting started
### Requirements
- A Fortran 2003 or better compliant compiler
- CMake version 3.14 or newer
- The [fypp preprocessor](https://github.com/aradi/fypp.git)
- The [Fortran standard library](https://github.com/fortran-lang/stdlib.git)
### Get code

``` sh
git clone https://github.com/Jim-215-Fisher/statislib.git
cd statislib
```
### Build with CMake

``` sh
cmake -B build
cmake --build build
```

Your compiled statislib will be in `statislib/build/src` as `libstatislib.a`, and you can copy the library to the selected place. To test, run each test file under `tests/` directory.

### Usage
In your program, access the library through `use` statement. You can use local name as short hand for library procedures or functions as follows:

``` fortran
program test
    use statislib, only : uni => uniform_distribution_rvs
    ...
```

## Documentation

Detail usages are listed in files for each distribution in `doc/`

## TODO

More distributions and functions will be implementated gradually.
