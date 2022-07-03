mkdir build
rm -rf build/*
cd build
cmake ..
make -j$(nproc)
./teste michigan 2000 10 10 0.8 0 3 0 0
cd ..
date
