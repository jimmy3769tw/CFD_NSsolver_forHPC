#pragma once
#include "../import/stl.hpp"
#include <string>
#include <sstream>

class IOassistant{

    int reserveSize_;
    int prevoiusNextStep_;
    int initStep_;
    std::string fileName_;

    std::vector<std::string> string_;

    public:
    IOassistant(std::string fileName, int resizeSize, int initStep)
    :reserveSize_(resizeSize), initStep_(initStep), fileName_(fileName), prevoiusNextStep_(initStep)
    {
        string_.resize(resizeSize);
    }

    void outPutFile(int current){
        std::ofstream file;
        file.open (fileName_, std::ios::out|ios::app);
        for(int i = prevoiusNextStep_; i <= current ; i++){
            file << string_[i-initStep_];
        }
        file.close();
        prevoiusNextStep_ = current+1;
    }

    std::string& operator[](int current){
        return string_[current-initStep_];
    }

    // initStep_ = 1
    // prevoiusNextStep_ = 1
    // 1 0 
    // 2 1 
    // 3 2
    // 4 3
    // write (4) 1 2 3 4 -> 0 1 2 3
    // prevoiusNextStep_ = 5
    // 5 4
    // 6 5
    // 7 6
    // write (7) 5 6 7 -> 4 5 6
};
