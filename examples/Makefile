CC=g++

SOURCES=$(TARGET).cpp 

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=$(TARGET)

#CFLAGS := -Wall -g -O3 -Wfatal-errors -mssse3 -msse2 -march=native
CFLAGS := -ggdb 
#LDFLAGS :=  grief/grief.cpp -I   -Llib -Wl,rpath='$$ORIGIN/lib' -g `pkg-config opencv --libs`  É O CERTO
LDFLAGS :=  ../DE/DE.cc ../grief.cc -I/usr/include/python3.8 -I~/.local/lib/python3.8/site-packages/numpy/core/include -Llib -Wl,-g `pkg-config opencv --libs` -lpython3.8



# OpenCV
CFLAGS += `pkg-config opencv --cflags`

all: convergence_test.cc ../grief.cc 
	$(CC) convergence_test.cc $(CFLAGS) $(LDFLAGS) -o $@


grief: teste_grief.cc grief.cc 
	$(CC) teste_grief.cc $(CFLAGS) $(LDFLAGS) -o $@

grief_de: teste_grief.cc grief.cc 
	$(CC) teste_grief.cc DE.cc $(CFLAGS) $(LDFLAGS) -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) 
	rm -f match_all evolve_grief annotate 

cleanall:
	rm -f *.o 


 