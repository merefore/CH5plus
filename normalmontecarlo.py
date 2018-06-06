import numpy as np
import timeit
from numpy.random import uniform as unirand

START = timeit.default_timer()

p=2

if p==1: ## Cs MIN
	w = [3132.850210,3001.025576,3224.281611,2707.756165,2417.807304,839.630358,1586.862768,1477.703750,1296.722907,1303.486805,1499.919440,199.670987,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000]

if p==2: ## C2V
	w = [3239.048896,3121.942709,2847.899544,2685.012264,2639.841589,497.745170,0.000000,1332.553861,1261.151760,1393.160933,1476.316601,1452.062281,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000]

if p==3: ## Cs SP
	w = [3236.059015,3037.762896,2728.879739,3090.466356,2400.988472,1603.353909,995.837684,1340.351157,1507.094234,1151.168976,1482.097119,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000]

lim = []
for i in xrange(0,len(w)):
	if w[i] == 0:
		twosig = 0.02
	else:
		twosig = 2/(np.sqrt(2*w[i]))
	lim.append(twosig)

def P(qs):
	product = 1
	for i in xrange(0,len(qs)):
		phi = np.exp(-w[i]*(qs[i]**2))
		product = product * phi
	return product

totaltime = 2000000
keepers = []

for t in xrange(0,totaltime):
	displace = np.zeros(18)
	for i in xrange(0,len(displace)):
		displace[i] = unirand(-lim[i],lim[i])	
	prob = P(displace)
	key = unirand()
	if prob>key:
		keepers.append(displace)

#output = open("MCkeepersC2V.txt", "a+")

#for i in xrange(0,len(keepers)):
#	for j in xrange(0,len(keepers[0])):
#		if j < (len(keepers[0])-1):
#			output.write("%f       " % keepers[i][j])
#		else: 
#			output.write("%f\n" % keepers[i][j])

#output.close()

print len(keepers)

STOP = timeit.default_timer()
print "Time to calculate: ", STOP - START, " seconds"
