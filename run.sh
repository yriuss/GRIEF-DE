mkdir build
rm -rf build/*
cd build
cmake ..
make -j$(nproc)
# dataset > gen number > experiments number > k value, cross value  
./teste michigan 10 5 10 0.8
cd ..
date
