#!/bin/bash

#> monolis
cd submodule/monolis
make clean
git checkout .
git checkout master
git pull
./install_lib.sh METIS
make FLAGS=MPI,METIS

