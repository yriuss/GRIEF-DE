mkdir build
rm -rf build/*
cd build
cmake ..
make -j$(nproc)

# Teste 7
./teste michigan 2000 5 10 0.8 0 4 1 0 0

cd ..
date
