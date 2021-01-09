#!/usr/bin/env bash

g++ -std=c++14 "$(pkg-config --cflags eigen3)" -o main bspline_plot.cpp

./main > log
for i in `seq 0 7`
do
    cat log | grep "^$i " > $i.log
done

gnuplot -e "
    plot 0;
    replot \"0.log\" using 2:3;
    replot \"1.log\" using 2:3;
    replot \"2.log\" using 2:3;
    replot \"3.log\" using 2:3;
    replot \"4.log\" using 2:3;
    replot \"5.log\" using 2:3;
    replot \"6.log\" using 2:3;
    replot \"7.log\" using 2:3;
    pause mouse, key;
"

for i in `seq 0 7`
do
    rm $i.log
done
