CC=g++
NVCC=nvcc
SOURCES=$(TARGET).cpp 

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=$(TARGET)
#CFLAGS := -Wall -g -O3 -Wfatal-errors -mssse3 -msse2 -march=native
CFLAGS := -ggdb 
#LDFLAGS :=  grief/grief.cpp -I   -Llib -Wl,rpath='$$ORIGIN/lib' -g `pkg-config opencv --libs`  É O CERTO
LDFLAGS :=  DE/DE.cc GRIEF/grief.cc -I/usr/include/python3.8 -I/home/adriel/.local/lib/python3.8/site-packages/numpy/core/include -Llib -Wl,-g `pkg-config opencv --libs` -lpython3.8
NVCC_RESULT := $(shell which nvcc 2> NULL)
NVCC_TEST := $(notdir $(NVCC_RESULT))


# OpenCV
CFLAGS += `pkg-config opencv --cflags`




	

all: main.cc GRIEF/grief.cc 
ifeq ($(NVCC_TEST), nvcc) 
	$(NVCC) main.cu $(LDFLAGS) -o bin/$@
else
	$(CC) main.cc $(CFLAGS) $(LDFLAGS) -o bin/$@
endif
	


grief: teste_GRIEF/grief.cc GRIEF/grief.cc 
	$(CC) teste_GRIEF/grief.cc $(CFLAGS) $(LDFLAGS) -o $@

grief_de: teste_GRIEF/grief.cc GRIEF/grief.cc 
	$(CC) teste_GRIEF/grief.cc DE.cc $(CFLAGS) $(LDFLAGS) -o $@

run:
	./bin/all GRIEF-datasets/planetarium

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) 
	rm -f match_all evolve_grief annotate 

cleanall:
	rm -f *.o 


 