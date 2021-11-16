#!/bin/bash

make

./Tracer $1 $2 $3

python3 imagize.py
