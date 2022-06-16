# git pull
rm -rf build/*
cd build
cmake ..
make -j$(nproc)
# dataset > gen number > individuals number > experiments number > k value, cross value  
./teste michigan 3 4 2 10 0.8
cd ..