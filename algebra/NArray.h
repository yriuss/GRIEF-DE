#ifndef NARRAY_H_INCLUDED
#define NARRAY_INCLUDED

#include <iostream>
#include <vector>
#include<cmath>

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

    NArray& operator[](size_t index){
        if(index < shape[times_accessed]){
            if(!wait_operation){
                times_accessed++;
                //std::cout << access << std::endl;
                if(shape.length() > 1){
                    access += index * pow(shape[shape.length() - 1],(shape.length() - times_accessed));
                    //std::cout << access << std::endl;
                }else
                    access += index;
                //std::cout << shape[shape.length() - 1] << " " << index << " " << access << " " << (int)(index * pow(shape[shape.length() - 1],(times_accessed - 1))) << std::endl;
                if(times_accessed == shape.length()){
                    times_accessed = 0;
                    wait_operation = true;
                }
            }else{
                n_access_waiting_mode++;
            }
            return *this;
        }else{
            std::cout << "Invalid index!" << std::endl;
            exit(-1);
        }
    }

    operator double() {
        double a = data[access];
        operation_done();
        return a;
    }

    operator float() {
        float a = data[access];
        operation_done();
        return a;
    }

    operator int() {
        int a = data[access];
        operation_done();
        return a;
    }

    NArray operator*(float number){
        std::vector<int> N;
        for(int i = 0; i < shape.length(); i++)
            N.push_back(shape[i]);
        NArray aux(N);
        aux = *this;
        for(int i = 0; i < size; i++)
            aux.data[i] *= number;
        return aux;
    }

    NArray operator+(const NArray& other){
        std::vector<int> N;
        for(int i = 0; i < shape.length(); i++)
            N.push_back(shape[i]);
        NArray aux(N);
        aux = *this;
        for(int i = 0; i < size; i++)
            aux.data[i] += other.data[i];
        return aux;
    }

    NArray operator-(const NArray& other){
        std::vector<int> N;
        for(int i = 0; i < shape.length(); i++)
            N.push_back(shape[i]);
        NArray aux(N);
        aux = *this;
        for(int i = 0; i < size; i++)
            aux.data[i] -= other.data[i];
        return aux;
    }

    friend std::ostream& operator<<(std::ostream& os, NArray& instance){
        if(instance.n_access_waiting_mode == instance.shape.length() || instance.n_access_waiting_mode == 0){
            if(!instance.wait_operation){
                if(instance.shape.length() <= 4){//needs to be 2d at maximum
                    os<<"[";
                    if(instance.shape.length() > 1){
                        for(int i=0; i < instance.size; i++){
                            os << instance.data[i];
                            if(i < instance.size - 1){
                                if(!((i + 1)%instance.shape[1]))
                                    os << ";";
                                else
                                    os << ",";
                            }
                        }
                        
                    }else{
                        for(int i=0; i < instance.shape[0] - 1; i++){
                            os << instance.data[i] << ",";
                        }
                        os << instance.data[instance.shape[0] - 1];
                    }
                    os<<"]";
                }else{
                    std::cout << "Not implemented yet!" << std::endl;
                }
            }else
                os << instance.data[instance.access];
        }else{
            std::cout << "Invalid format! " << instance.n_access_waiting_mode << std::endl;

            exit(-1);
        }
        instance.operation_done();
        return os;
    }

    void operator=(NArray array){
        
        for(int i = 0; i < array.size; i++){
            //std::cout << data[i + access] << std::endl;
            data[i + access] = array.data[i];
            
        }
        operation_done();
    }

    void operator=(int element){
        //std::cout << element << "for access" << access << std::endl;
        if(n_access_waiting_mode == shape.length() || n_access_waiting_mode == 0){
            //std::cout << access << std::endl;
            data[access] = element;
        }else{
            std::cout << "Invalid Operation " << n_access_waiting_mode << std::endl;
            exit(-1);
        }
        operation_done();
    }

    void print(){
        std::cout << data[0];
    }

    void operation_done(){
        access = 0;
        wait_operation = false;
        n_access_waiting_mode = 0;
    }
private:
    float* data = nullptr;
    int prev_index = -1;
    size_t size = 1;
    size_t access = 0;
    size_t times_accessed = 0;
    size_t n_access_waiting_mode = 0;
    bool wait_operation = false;
};





#endif