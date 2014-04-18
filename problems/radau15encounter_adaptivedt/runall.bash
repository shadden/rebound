#!/bin/bash
rm *.txt
make radau15
./nbody --epsilon=1e-2
mv dt.txt dt_1e-2.txt
./nbody --epsilon=1e-4
mv dt.txt dt_1e-4.txt
./nbody --epsilon=1e-6
mv dt.txt dt_1e-6.txt
./nbody --epsilon=1e-8
mv dt.txt dt_1e-8.txt
