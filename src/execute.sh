#!/bin/bash

make

./Tracer $1 $2 "../data/P.ppm"

open ../data/P.ppm
