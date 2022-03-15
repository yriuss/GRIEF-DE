CC=g++

SOURCES=$(TARGET).cpp 

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=$(TARGET)

#CFLAGS := -Wall -g -O3 -Wfatal-errors -mssse3 -msse2 -march=native
CFLAGS := -ggdb 
#LDFLAGS :=  grief/grief.cpp -I   -Llib -Wl,rpath='$$ORIGIN/lib' -g `pkg-config opencv --libs`  É O CERTO
LDFLAGS :=  grief.cc -I   -Llib -Wl,-g `pkg-config opencv --libs`



# OpenCV
CFLAGS += `pkg-config opencv --cflags`

all: teste_grief.cpp grief.cpp 
	$(CC) teste_grief.cpp $(CFLAGS) $(LDFLAGS) -o $@


evolve_grief: teste_grief.cc grief.cc 
	$(CC) teste_grief.cc DE.cc $(CFLAGS) $(LDFLAGS) -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) 
	rm -f match_all evolve_grief annotate 

cleanall:
	rm -f *.o 


 