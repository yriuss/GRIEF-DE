git pull
rm -rf build/*
cd build
cmake .. >> log
make >> log
./teste michigan 200 8 3 >> log
cd ..