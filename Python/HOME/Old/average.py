#!/usr/bin/python

import numpy as np

# os.chdir("/home/eixeres/Downloads/Berendsen")
input = np.loadtxt('Contact_Angles.txt')
tinput = np.transpose(input)
y2 = tinput[2]
x2 = tinput[3]
i2 = len(y2)
x = np.zeros(i2)
y = np.zeros(i2)



for i in range(0,i2):
    x[i] = (x2[i])**(-1)
    y[i] = np.cos(y2[i])
    