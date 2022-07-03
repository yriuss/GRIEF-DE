mkdir build
rm -rf build/*
cd build
cmake ..
make -j$(nproc)
# Teste 4
./teste michigan 2000 10 10 0.8 1 4 0 0
cd ..
date
