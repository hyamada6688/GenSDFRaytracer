#!/usr/bin/env python3

import numpy as np
from PIL import Image
#import numdifftools as nd
import matplotlib.pyplot as plt
import pandas as pd

raw = pd.read_csv("../data/small.csv", delimiter=";", header = None)

data =[]
for x in range(0, len(raw)):
    thisrow = []
    for y in range(len(raw.iloc[x]) - 1):
        #print(raw.iloc[x].iloc[y])
        stri = raw.iloc[x].iloc[y]
        #print(stri)
        #print(stri.split(', '))
        pixel = []
        for y in stri.split(','):
            if y =='nan':
                y = 0
            pixel.append(int(float(y)))
        thisrow.append(pixel)
   # print(thisrow)
    data.append(np.array(thisrow))

data = np.array(data, dtype=np.uint8)
img = Image.fromarray(data, 'RGB')
img.save('../data/LastRender.png')
img.show()
