#!/bin/bash

mkdir lib
mkdir include
mkdir include/BB
mkdir include/BBFE

cd libBB
make clean
make
cp *.h ../include/BB
cd ..

cd FE_std
make clean
make
cp *.h ../include/BBFE
cd ..

cd FE_sys
make clean
make
cp *.h ../include/BBFE
cd ..

cd FE_elemmat
make clean
make
cp *.h ../include/BBFE
cd ..

cd FE_manusol
make clean
make
cp *.h ../include/BBFE
cd ..

make clean
make
