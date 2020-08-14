#!/bin/bash

#> monolis
cd submodule/monolis
make clean
git checkout .
git checkout master
git pull
./install_lib.sh
make FLAGS=MPI,METIS

