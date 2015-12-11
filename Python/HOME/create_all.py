#!/usr/bin/env python

import numpy as np
import os

StartingPercentage = 0
EndingPercentage = 30

for CurPercentage in np.arange(StartingPercentage,EndingPercentage+1):
  print "now doing " + str(CurPercentage) + "%"
  os.system('python laila.py ' + str(CurPercentage))

