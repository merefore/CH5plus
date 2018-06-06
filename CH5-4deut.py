import timeit
import math
import itertools
import random
import numpy as np
# import matplotlib.pyplot as mpl
import CH5pot

## - IS THIS FOR MOX??? No: 1, Yes: 2
mox = 2

START = timeit.default_timer()

n0 = 20000 #number of initial walkers
dtau = 10 #time step
m1 = 21874.658 #C
m2 = 1837.152 #H
m3 = 3671.482 #D
sigma1 = np.sqrt(dtau / m1)
sigma2 = np.sqrt(dtau / m2)
sigma3 = np.sqrt(dtau / m3)
cycles = 10000
alpha = 0.5 / dtau

def V(whereUat): #-----------our potential function
	v = np.zeros(len(whereUat))
	emin = -40.65276470207075 
	v = CH5pot.mycalcpot(whereUat,len(whereUat))
	v = np.array(v)
	return v # must return a 1D numpy array
	
def avg(potentials, init): #-------------our vref function
	potTotal = potentials.sum()
	vref = potTotal / float(len(potentials)) - alpha*(len(potentials)-init)/float(init)
	return vref
	
def pcalc(potentials, vref, dtau): #------------to calculate pvalues
	a = -(potentials - vref) * dtau
	pvalues = np.exp(a)
	return pvalues
	
def Fate(whereUat,whereUfrom,potentials,pvalues): #--------will return our new collection of walkers 
	positions2 = []
	potentials2 = []
	elders2 = []

	whole = pvalues.astype(int)
	frac = pvalues - whole
	Nr = np.random.random(len(whereUat))
	Nr2 = np.random.random(len(whereUat))

	for i in range(0,len(whole)):
		if potentials[i] > 0:
			if whole[i] != 0:
				for h in range(0,whole[i]):
					positions2.append(whereUat[i])
					potentials2.append(potentials[i])
					elders2.append(whereUfrom[i])
			if Nr[i] < frac[i]:
				positions2.append(whereUat[i])
				potentials2.append(potentials[i])
				elders2.append(whereUfrom[i])
	
	positions2 = np.array(positions2)
	potentials2 = np.array(potentials2)
	elders2 = np.array(elders2)
	return positions2,potentials2,elders2
	
#Initialize Walkers 


# - MIN ENERGY CONFIG
beginning = np.array([[0.000000000000000, 0.000000000000000, 0.000000000000000],
[0.1318851447521099, 2.088940054609643, 0.000000000000000],
[1.786540362044548, -1.386051328559878,   0.000000000000000],
[2.233806981137821, 0.3567096955165336, 0.000000000000000],
[-0.8247121421923925, -0.6295306113384560, -1.775332267901544],
[-0.8247121421923925, -0.6295306113384560, 1.775332267901544]])

# - C2V CONFIG
#beginning = np.array([[0.000000000000000, 0.000000000000000, 0.3869923621587414],
#[0.000000000000000, 0.000000000000000, -1.810066283748844],
#[1.797239666982623, 0.000000000000000,   1.381637275550612],
#[-1.797239666982623, 0.000000000000000, 1.381637275550612],
#[0.000000000000000, -1.895858229423645, -0.6415748897955779],
#[0.000000000000000, 1.895858229423645, -0.6415748897955779]])

whereUat = np.zeros((n0,beginning.shape[0],beginning.shape[1]))

for i in range(0,n0):
        whereUat[i] = beginning * 1.1
		
# potentials = V(whereUat)	

#------    FILENAMES   ---------
#######                                          #########
#######                                          #########
#######                                          #########
if mox == 1:
	fileEnergy = open("5deut_energies.txt", "a+")
	fileXYZ = open("5deut_ancestors.xyz", "a+")
timestep = 0
energies = []
population = []
ancestorPile = []

#BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM

for c in range(0,cycles):    

	if timestep%500 == 0:
		ancestors = whereUat
		whereUfrom = np.arange(0,len(whereUat)) #initialize the "index array", which will reference indexes in "ancestors"
		
	# Create Random Displacements # (Remember the mass!)
	displace = np.random.normal(0,sigma3,whereUat.shape)
	for i in range(0,len(displace)): 
		displace[i,0,:] = np.random.normal(0,sigma1,(1,3))
		displace[i,1,:] = np.random.normal(0,sigma2,(1,3))
		
	
	# Displace the Walkers #
	whereUat = whereUat + displace
	
	# Calculate the Potential #
	potentials = V(whereUat) #using my function V 
	
	# Take Vref #
	vref = avg(potentials, n0) #umf avg
	energies.append(vref)
	if mox == 1:
		fileEnergy.write("%f\n" % vref)
	population.append(len(whereUat))
	
	# Calculate PValues #
	pvalues = pcalc(potentials, vref, dtau) #umf pcalc
	whereUat,potentials,whereUfrom = Fate(whereUat,whereUfrom,potentials,pvalues) #umf Fate
	
	if timestep%500 == 0:
		print vref, ", ", timestep, ", ", len(whereUat)
	
	# Descendant Weighting #	
	if timestep%500 == 50 and timestep > 4000:
		for i in whereUfrom:
			ancestorPile.append(ancestors[i])

	timestep += 1
ancestorPile = np.array(ancestorPile)	

b=0
convert = 0.52917725
v = np.zeros(len(ancestorPile))
v = CH5pot.mycalcpot(ancestorPile,len(ancestorPile))

for item in ancestors: #-----------for jmol xyz file
	if mox == 1:	
		fileXYZ.write("6\n")
		fileXYZ.write("CH5+ coordinates         %f" % v[b])
		fileXYZ.write("%i\n" % b)
		fileXYZ.write("C       %f       %f       %f\n" % (item[0,0]*convert, item[0,1]*convert, item[0,2]*convert))
		fileXYZ.write("H       %f       %f       %f\n" % (item[1,0]*convert, item[1,1]*convert, item[1,2]*convert))
		fileXYZ.write("H       %f       %f       %f\n" % (item[2,0]*convert, item[2,1]*convert, item[2,2]*convert))
		fileXYZ.write("H       %f       %f       %f\n" % (item[3,0]*convert, item[3,1]*convert, item[3,2]*convert))
		fileXYZ.write("H       %f       %f       %f\n" % (item[4,0]*convert, item[4,1]*convert, item[4,2]*convert))
		fileXYZ.write("H       %f       %f       %f\n\n" % (item[5,0]*convert, item[5,1]*convert, item[5,2]*convert))
		b+=1
	else:
		print "6\n"
		print "CH5+ coordinates         %f" % v[b]
		print "%i\n" % b
		print "C       %f       %f       %f\n" % (item[0,0]*convert, item[0,1]*convert, item[0,2]*convert)
		print "H       %f       %f       %f\n" % (item[1,0]*convert, item[1,1]*convert, item[1,2]*convert)
		print "H       %f       %f       %f\n" % (item[2,0]*convert, item[2,1]*convert, item[2,2]*convert)
		print "H       %f       %f       %f\n" % (item[3,0]*convert, item[3,1]*convert, item[3,2]*convert)
		print "H       %f       %f       %f\n" % (item[4,0]*convert, item[4,1]*convert, item[4,2]*convert)
		print "H       %f       %f       %f\n\n" % (item[5,0]*convert, item[5,1]*convert, item[5,2]*convert)
		b+=1

if mox == 1:
	fileEnergy.close()
	fileXYZ.close()

STOP = timeit.default_timer()
print "Time to calculate: ", STOP - START, " seconds"
avgEn =np.average(energies)
print "Energy: ", avgEn

# --------Energy and Population Tracking
# t = range(0,cycles)
# mpl.subplot(211)	
# mpl.plot(t, energies)
# mpl.subplot(212)
# mpl.plot(t, population, 'r')
# mpl.show()
