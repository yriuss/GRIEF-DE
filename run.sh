# git pull
rm -rf build/*
cd build
cmake ..
make 
# dataset > gen number > individuals number > experiments number > k value, cross value  
./teste michigan 5 5 5 10 0.8
cd ..