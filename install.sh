#!/bin/bash

mkdir lib
mkdir include
mkdir include/BB
mkdir include/BBFE
mkdir include/BBFE/std
mkdir include/BBFE/sys
mkdir include/BBFE/elemmat
mkdir include/BBFE/manusol

# install Bebops libraries

cd libBB
make clean
make
cp *.h ../include/BB
cp *.a ../lib/
cd ..

# install FE libraries

cd FE_std
make clean
make
cp *.h ../include/BBFE/std
cp *.a ../lib/
cd ..

cd FE_sys
make clean
make
cp *.h ../include/BBFE/sys
cp *.a ../lib/
cd ..

cd FE_elemmat
make clean
make
cp *.h ../include/BBFE/elemmat
cp *.a ../lib/
cd ..

cd FE_manusol
make clean
make
cp *.h ../include/BBFE/manusol
cp *.a ../lib/
cd ..
