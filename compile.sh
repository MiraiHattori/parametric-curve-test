#!/usr/bin/env bash

#CXX='g++-7'
#CC='gcc-7'
#
#CMAKE_ROOT='$HOME/.future731/share/cmake-3.12'
#CMAKE_BIN='$HOME/.future731/bin/cmake'

if [ $# -ne 1 ]; then
    echo "usage: ./compile.sh [filename]"
    exit -1
fi

# g++-7 -std=c++14 -I$HOME/.future731/include/ -L${HOME}/.future731/lib -o main $1 -I/usr/include/eigen3
g++ -std=c++11 -o main $1 -isystem /usr/include/eigen3
