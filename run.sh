mkdir build
rm -rf build/*
cd build
cmake ..
make -j$(nproc)
# Teste 3
./teste michigan 2000 10 10 0.8 0 4 0 0
cd ..
date
