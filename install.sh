#!/bin/bash

mkdir lib
mkdir include
mkdir include/BB
mkdir include/BBFE
mkdir include/BBFE/std
mkdir include/BBFE/sys
mkdir include/BBFE/elemmat
mkdir include/BBFE/manusol

cd libBB
make clean
make
cp *.h ../include/BB
cd ..

cd FE_std
make clean
make
cp *.h ../include/BBFE/std
cd ..

cd FE_sys
make clean
make
cp *.h ../include/BBFE/sys
cd ..

cd FE_elemmat
make clean
make
cp *.h ../include/BBFE/elemmat
cd ..

cd FE_manusol
make clean
make
cp *.h ../include/BBFE/manusol
cd ..

make clean
make
