#!/bin/bash

for i in $(ls); do

for j in $(ls $i); do

for k in $(ls $i/$j); do

cat $i/$j/$k/*/bilayer_position.dat | grep -v \#  | grep -v @ > ~/$i.$j.$k.posi
cat $i/$j/$k/*/bun_tilt.xvg | grep -v \#  | grep -v @ > ~/$i.$j.$k.tilt
cat $i/$j/$k/*/bun_kink.xvg | grep -v \#  | grep -v @ > ~/$i.$j.$k.kink
cat $i/$j/$k/*/rot.dat | grep -v \#  | grep -v @ > ~/$i.$j.$k.rota
cat $i/$j/$k/*/tilt.dat | grep -v \#  | grep -v @ > ~/$i.$j.$k.internal.tilt

done

done

done
