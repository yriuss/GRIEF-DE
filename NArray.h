#ifndef NARRAY_H_INCLUDED
#define NARRAY_INCLUDED

#include <iostream>
#include <vector>

class Shape{
public:
    Shape(){
        //allocate 2 elements
        ReAlloc(2);
    }
    
    void push_back(const int value){
        
        if(size > capacity)
            ReAlloc(capacity + capacity / 2);

        data[size++] = value;
    }

    const int operator[](size_t index) const{
        if(index > size){
            std::cout << "Invalid index!" << std::endl;
            exit(-1);
        }
        return data[index];
    }

    //int operator[](size_t index){
    //    std::cout << "Not implemented yet!" << std::endl;
    //    exit(-1);
    //}

    void operator=(std::vector<int> N){
        if(once == 0){
            for(int i = 0; i < N.size(); i++){
                this->push_back(N[i]);
            }
            once = 1;
        }else{
            std::cout << "Cannot reshape!" << std::endl;
            exit(-1);
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const Shape& instance){
        os << "shape: [";
        for(int i = 0; i < instance.size; i++){
            if(i < instance.size - 1)
                os << instance.data[i] << ",";
            else
                os << instance.data[i];
        }
        os << "]";
        return os;
    }

private:
    void ReAlloc(size_t newCapacity){
        //allocate a new block of memory
        int* newBlock = new int[newCapacity];

        if(newCapacity < size)
            size = newCapacity;

        //copy/move old elements into new block
        for(size_t i = 0; i < size; i++)
            newBlock[i] = data[i];
        
        //delete
        delete[] data;
        data = newBlock;
        capacity = newCapacity;
    }

    int* data = nullptr;
    size_t size = 0;
    size_t capacity = 0;
    bool once = 0;
};

class NArray{
public:
    NArray(std::vector<int> N);
    Shape shape;
private:
    
};

#endif