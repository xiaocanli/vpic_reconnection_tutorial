#!/bin/sh

for i in {0..255}
do
    # echo $i
    for j in {0..31}
    do
        rank=$((i * 32 + j))
        echo $rank
        mv restore1/tracer1.$rank restore1/$i
    done
    # mv restore1/$i/* restore1
    # ls $i/tracer1.*
    # mv $i/tracer1* tracer_restart
done
