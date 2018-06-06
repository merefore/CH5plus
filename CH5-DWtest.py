import timeit
import math
import itertools
import random
import numpy as np
# import matplotlib.pyplot as mpl
import CH5pot

START = timeit.default_timer()

n0 = 20000 #number of initial walkers
dtau = 10 #time step
m1 = 21894.687
m2 = 1837.515
sigma1 = np.sqrt(dtau / m1)
sigma2 = np.sqrt(dtau / m2)
cycles = 10000
alpha = 0.5 / dtau

deltaT = 1000

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
beginning = np.array([[0.000000000000000, 0.000000000000000, 0.3869923621587414],
[0.000000000000000, 0.000000000000000, -1.810066283748844],
[1.797239666982623, 0.000000000000000,   1.381637275550612],
[-1.797239666982623, 0.000000000000000, 1.381637275550612],
[0.000000000000000, -1.895858229423645, -0.6415748897955779],
[0.000000000000000, 1.895858229423645, -0.6415748897955779]])

whereUat = np.zeros((n0,beginning.shape[0],beginning.shape[1]))

for i in range(0,n0):
        whereUat[i] = beginning * 1.1	

#------    FILENAMES   ---------
#######                                          #########
#######                                          #########
#######                                          #########
#fileEnergy = open("energies.txt", "a+")
fileXYZ = open("DWtest1000.xyz", "a+")
#fileXYZ1 = open("dT1000-10ancestorsDC.xyz", "a+")
#fileXYZ2 = open("dT1000-25ancestorsDC.xyz", "a+")
#fileXYZ3 = open("dT1000-50ancestorsDC.xyz", "a+")
#fileXYZ4 = open("dT1000-100ancestorsDC.xyz", "a+")
timestep = 0
energies = []
population = []
ancestorPile = []
descendantPile1 = []
descendantPile2 = []
descendantPile3 = []
descendantPile4 = []


#BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM BEGIN PROGRAM

for c in range(0,cycles):    

	if timestep%deltaT == 0:
		ancestors = whereUat*1.0 #pointer in the right place? deep copy? 
		whereUfrom = np.arange(0,len(ancestors)) #initialize the "index array", which will reference indexes in "ancestors"
		if timestep >= .4*cycles:
			descendants = np.zeros(len(ancestors))				
			for i in range(0,len(ancestors)):
				ancestorPile.append(ancestors[i])
		
	# Create Random Displacements # (Remember the mass!)
	displace = np.random.normal(0,sigma2,whereUat.shape)
	for i in range(0,len(displace)): 
		displace[i,0,:] = np.random.normal(0,sigma1,(1,3))
		
	# Displace the Walkers #
	whereUat = whereUat + displace
	
	# Calculate the Potential #
	potentials = V(whereUat) #using my function V 
	
	# Take Vref #
	vref = avg(potentials, n0) #umf avg
	energies.append(vref)
	#fileEnergy.write("%f\n" % vref)
	population.append(len(whereUat))
	
	# Calculate PValues #
	pvalues = pcalc(potentials, vref, dtau) #umf pcalc
	whereUat,potentials,whereUfrom = Fate(whereUat,whereUfrom,potentials,pvalues) #umf Fate
	
	if timestep%(cycles/20) == 0:
		print vref, ", ", timestep, ", ", len(whereUat)

	if timestep%deltaT == 10 and timestep > .4*cycles:
		for index in whereUfrom:
			descendants[index] = descendants[index] + 1
		descendantPile1 = np.concatenate((descendantPile1, descendants))
		descendants = np.zeros(len(ancestors))
	if timestep%deltaT == 25 and timestep > .4*cycles:
		for index in whereUfrom:
			descendants[index] = descendants[index] + 1
		descendantPile2 = np.concatenate((descendantPile2, descendants))
		descendants = np.zeros(len(ancestors))
	if timestep%deltaT == 50 and timestep > .4*cycles:
		for index in whereUfrom:
			descendants[index] = descendants[index] + 1
		descendantPile3 = np.concatenate((descendantPile3, descendants))
		descendants = np.zeros(len(ancestors))
	if timestep%deltaT == 99 and timestep > .4*cycles:
		for index in whereUfrom:
			descendants[index] = descendants[index] + 1
		descendantPile4 = np.concatenate((descendantPile4, descendants))

	timestep += 1

ancestorPile = np.array(ancestorPile)

b=0
convert = 0.52917725
#fileXYZ = fileXYZ1
v = np.zeros(len(ancestorPile))
v = CH5pot.mycalcpot(ancestorPile,len(ancestorPile))
for i in range(0,len(ancestorPile)): #-----------for jmol xyz file
	fileXYZ.write("6\n")
	fileXYZ.write("Walker#, potential, descendants [10,25,50,100]         %i          %f          %i          %i          %i          %i\n" % (b, v[i], descendantPile1[i], descendantPile2[i], descendantPile3[i], descendantPile4[i]))
	fileXYZ.write("C       %f       %f       %f\n" % (ancestorPile[i,0,0]*convert, ancestorPile[i,0,1]*convert, ancestorPile[i,0,2]*convert))
	fileXYZ.write("H       %f       %f       %f\n" % (ancestorPile[i,1,0]*convert, ancestorPile[i,1,1]*convert, ancestorPile[i,1,2]*convert))
	fileXYZ.write("H       %f       %f       %f\n" % (ancestorPile[i,2,0]*convert, ancestorPile[i,2,1]*convert, ancestorPile[i,2,2]*convert))
	fileXYZ.write("H       %f       %f       %f\n" % (ancestorPile[i,3,0]*convert, ancestorPile[i,3,1]*convert, ancestorPile[i,3,2]*convert))
	fileXYZ.write("H       %f       %f       %f\n" % (ancestorPile[i,4,0]*convert, ancestorPile[i,4,1]*convert, ancestorPile[i,4,2]*convert))
	fileXYZ.write("H       %f       %f       %f\n\n" % (ancestorPile[i,5,0]*convert, ancestorPile[i,5,1]*convert, ancestorPile[i,5,2]*convert))
	b+=1


#fileEnergy.close()
#fileXYZ1.close()
#fileXYZ2.close()
#fileXYZ3.close()
#fileXYZ4.close()
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
