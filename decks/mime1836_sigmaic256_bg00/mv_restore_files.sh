#!/bin/sh

for i in {0..63}
do
    echo $i
    mv restore1/$i/restore.0* restore1
    # ls $i/tracer1.*
    # mv $i/tracer1* tracer_restart
done
