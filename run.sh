mkdir build
rm -rf build/*
cd build
cmake ..
make -j$(nproc)

# Teste 10
./teste michigan 2000 5 10 0.8 1 3 1 0 0

cd ..
date
