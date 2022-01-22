#!/bin/bash

make

./Tracer $1 $2 "../data/small.csv"

python3 imagize.py
