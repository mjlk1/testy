#!/usr/bin/env python3

import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec

if len(sys.argv)<=1:
	print("usage: bulletplot.py <file> ...")
	exit(1)

fig = plt.figure()
gs = gridspec.GridSpec(2,2)
axXY = fig.add_subplot(gs[:,0])
axXY.set_title("y(x)")
axTV = fig.add_subplot(gs[:,1])
axTV.set_title("v(t)")

for i in range(1,len(sys.argv)):
	data = np.loadtxt(sys.argv[i])
	t = data[:,0]
	x = data[:,1]
	y = data[:,2]
	vx = data[:,3]
	vy = data[:,4]
	v = np.sqrt(vx**2+vy**2)
	axXY.plot(x,y,label=sys.argv[i],marker='x')
	axTV.plot(t,v,label=sys.argv[i],marker='x')

axXY.legend()
axTV.legend()
plt.show()

