#!/usr/bin/env python3


import sys
import numpy as np
from PIL import Image
#import numdifftools as nd
import matplotlib.pyplot as plt
import pandas as pd

print("Rendering image...")
raw = pd.read_csv("../data/small.csv", delimiter=";", header = None)

data =[]

print(len(raw))

step = (int)(len(raw)/100)
perc = 0
print(step)
for x in range(0, len(raw)):
    thisrow = []
    if (step > 1 and x%step == 0): 
        print(str(perc) + '%')
        perc += 1
    for y in range(len(raw.iloc[x]) - 1):
        stri = raw.iloc[x].iloc[y]
        pixel = []
        for y in stri.split(','):
            if y =='nan':
                y = 0
            pixel.append(int(float(y)))
        thisrow.append(pixel)
    data.append(np.array(thisrow))

data = np.array(data, dtype=np.uint8)
img = Image.fromarray(data, 'RGB')
img.save('../data/LastRender.png')
print("Voila")
img.show()
