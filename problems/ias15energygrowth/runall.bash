#!/bin/bash

rm energy_*.txt

make ias15
./nbody --epsilon=1e0 --e=0.1
echo " ">> energy_ias15.txt
./nbody --epsilon=1e-4 --e=0.1
mv energy_ias15.txt energy_ias15_variable.txt
./nbody --timesteps=100000 --e=0.1

make wh
./nbody --timesteps=100000 --e=0.1
