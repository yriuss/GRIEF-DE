#ifndef NARRAY_H_INCLUDED
#define NARRAY_INCLUDED

#include <iostream>
#include <vector>

class Array{
public:
    Array(){
        //allocate 2 elements
        ReAlloc(2);
    }
    
    void push_back(const int value){
        
        if(size > capacity)
            ReAlloc(capacity + capacity / 2);

        data[size++] = value;
    }

    const Array& operator[](Array a) const{
        return *this;
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

    friend std::ostream& operator<<(std::ostream& os, const Array& instance){
        os << "[";
        for(int i = 0; i < instance.size; i++){
            if(i < instance.size - 1)
                os << instance.data[i] << ",";
            else
                os << instance.data[i];
        }
        os << "]";
        return os;
    }

    size_t length(){
        return size;
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
    Array shape;
    //const int operator[](size_t index) const{
    //    if(index > size){
    //        std::cout << "Invalid index!" << std::endl;
    //        exit(-1);
    //    }
    //    return data[index];
    //}

    NArray& operator[](size_t index){
        if(!wait_operation){
            times_accessed++;

            access += times_accessed*(index);
            if(times_accessed == shape.length()){
                times_accessed = 0;
                wait_operation = true;
            }
        }
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, NArray& instance){
        if(instance.times_accessed == 0){
            os << instance.data[instance.access];
        }
        instance.access = 0;
        instance.wait_operation = false;
        return os;
    }

    void operator=(int element){
        if(times_accessed != 0){
            std::cout << "Invalid Operation" << std::endl;
            exit(-1);
        }
        
        data[access] = element;
        access = 0;
        wait_operation = false;
        
    }

    void print(){
        std::cout << data[0];
    }

private:
    int* data = nullptr;
    size_t size = 1;
    size_t access = 0;
    size_t times_accessed = 0;
    bool wait_operation = false;
};

#endif