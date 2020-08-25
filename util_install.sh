#!/bin/bash

mkdir bin

cd util/meshgen
make clean
make
cp meshgen_tet ../../bin
cp meshgen_hex ../../bin
cd ../../

cd util/converters
make clean
make
cp vtk2dat_tet ../../bin
cp vtk2dat_hex ../../bin
cd ../../

cd util/surface
make clean
make
cp surf_dbc_all ../../bin
cp surf_nbc ../../bin
cp surf_dbc ../../bin
cp surf_conn ../../bin
cp surf_bc_merge ../../bin
cd ../../

cd util/cmd2cond
make clean
make
cp cmd2cond ../../bin
cd ../../

cd util/mesh
make clean
make
cp mesh_remove ../../bin
cp mesh_extract ../../bin
cp mesh_surf_remove ../../bin
cp mesh_surf_extract ../../bin
cd ../../
