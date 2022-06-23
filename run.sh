mkdir build
rm -rf build/*
cd build
cmake ..
make -j$(nproc)
# dataset > gen number > experiments number > k value, cross value  
./teste michigan 2000 5 15 0.8
cd ..