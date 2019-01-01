#!/usr/bin/env python3

import sys
import numpy as np

if len(sys.argv)<=2:
	print("usage: cmp.py <reference> <test0> ...")
	exit(1)

refxy = None

testxy = []

for i in range(1,len(sys.argv)):
	data = np.loadtxt(sys.argv[i])
	t = data[:,0]
	x = data[:,1]
	y = data[:,2]
	vx = data[:,3]
	vy = data[:,4]
	v = np.sqrt(vx**2+vy**2)
	if i==1:
		refxy = np.array([x[-1],y[-1]])
	else:
		testxy.append([x[-1],y[-1]])

testxy = np.array(testxy)

err = testxy-refxy
err = np.sum(err**2,axis=1)

for x in err:
	print(x)
