#!/bin/bash

cd ~/Documents/MPhysProject
echo "plotting!"

gfortranDev ChargeDensitySampler.f90
./a.out
python plot.py

echo "done!"
