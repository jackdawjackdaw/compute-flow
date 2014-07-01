Flow Compute
===========

ccs, July-1st-2014
cec24@phy.duke.edu

Compiled routines to compute ptBinned v2 and v3.


# Requires

* CMAKE
* Check (for tests)
* A C compiler
* libGSL


## Building & Installing

This project uses CMAKE to generate Makefiles, it is canonical to do out of place builds using cmake. An "out of place" builds puts all the temporary files and compiler junk into a directory that is outside the source tree.

From the project root do:

    mkdir ./build
    cd ./build
    cmake ..
    make && make install

Cmake defaults to installing things in /usr/local, if you don't want that you should set invoke cmake as

    cmake -DCMAKE_INSTALL_PREFIX:PATH=/your/install/path ..


