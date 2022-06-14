# git pull
rm -rf build/*
cd build
cmake ..
make 
# dataset > gen number > individuals number > experiments number > k value, cross value  
./teste michigan 200 10 5 10 0.8
cd ..