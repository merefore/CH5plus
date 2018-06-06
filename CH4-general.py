import timeit
import math
import itertools
import random
import numpy as np
# import matplotlib.pyplot as mpl
import CH4pot

START = timeit.default_timer()

n0 = 20000 #number of initial walkers
dtau = 1 #time step
m1 = 21874.658
m2 = 1837.152
sigma1 = np.sqrt(dtau / m1)
sigma2 = np.sqrt(dtau / m2)
cycles = 10000
alpha = 0.5 / dtau

convert = 0.52917725

def V(whereUat): #-----------my potential function, must return a 1D numpy array
	newcoord = np.zeros((len(whereUat),10))
	for i in xrange(0,len(whereUat)):
		
		xh = np.zeros(4)
		yh = np.zeros(4)
		zh = np.zeros(4)
		xc = whereUat[i,0,0]
		yc = whereUat[i,0,1]
		zc = whereUat[i,0,2]
		
		for j in xrange(0,4):
			xh[j] = whereUat[i,j+1,0]
			yh[j] = whereUat[i,j+1,1]
			zh[j] = whereUat[i,j+1,2]
			newcoord[i,j] = np.sqrt((xh[j]-xc)**2+(yh[j]-yc)**2+(zh[j]-zc)**2)*convert 
		
		count = 4
		
		for k in range(0,3):
			u = [xh[k], yh[k], zh[k]]
			umag = np.linalg.norm(u)
			for m in range(k+1,4):
				v = [xh[m], yh[m], zh[m]]
				vmag = np.linalg.norm(v)
				dotproduct = np.dot(u,v)
				newcoord[i,count] = (180.00/np.pi)*np.arccos(dotproduct/(umag*vmag))
				count = count + 1
	
	potentials = np.zeros(len(whereUat)) 
	force = np.zeros(289)
	force[0] = 5
	potentials[0] = 5
	potentials,force = CH4pot.walkercalcpot(newcoord,len(newcoord))
	potentials = np.array(potentials)
	potentials = potentials*0.000004556335
	return potentials 
	
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
		if 1==1: #potentials[i] > 0:
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
beginning = np.array([[0.000000000000000, 0.000000000000000, 0.000000000000000],
[-1.6756644725, 0.000000000000000, -1.1848737116],
[1.6756644725, 0.000000000000000,  -1.1848737116],
[0.000000000000000, -1.6756644725, 1.1848737116],
[0.000000000000000, 1.6756644725, 1.1848737116]])

#beginning = np.array([[0.000000000000000, 0.000000000000000, 0.000000000000000],
#[-0.8, 0.000000000000000, -1.5848737116],
#[0.8, 0.000000000000000,  -1.5848737116],
#[0.000000000000000, -0.8, 1.548737116],
#[0.000000000000000, 0.8, 1.548737116]])

whereUat = np.zeros((n0,beginning.shape[0],beginning.shape[1]))

for i in range(0,n0):
        whereUat[i] = beginning * 1.15
		
# potentials = V(whereUat)	

#------    FILENAMES   ---------
#######                                          #########
#######                                          #########
#######                                          #########
fileEnergy = open("nokillCH4energies.txt", "a+")
fileXYZ = open("nokillCH4ancestors.xyz", "a+")
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
	fileEnergy.write("%f\n" % vref)
	population.append(len(whereUat))
	
	# Calculate PValues #
	pvalues = pcalc(potentials, vref, dtau) #umf pcalc
	whereUat,potentials,whereUfrom = Fate(whereUat,whereUfrom,potentials,pvalues) #umf Fate
	
	if timestep%500 == 0:
		print vref, ", ", timestep, ", ", len(whereUat)

	if timestep%500 == 50 and timestep > 4000:
		for i in whereUfrom:
			ancestorPile.append(ancestors[i])

	timestep += 1
ancestorPile = np.array(ancestorPile)

b=0
v = np.zeros(len(ancestorPile))
#v = CH5pot.mycalcpot(ancestorPile,len(ancestorPile))
#for i in range(0,len(ancestorPile)): #-----------for jmol xyz file
#	fileXYZ.write("6\n")
#	fileXYZ.write("CH5+ coordinates         %i          %f\n" % (b, v[i]))
#	fileXYZ.write("C       %f       %f       %f\n" % (ancestorPile[i,0,0]*convert, ancestorPile[i,0,1]*convert, ancestorPile[i,0,2]*convert))
#	fileXYZ.write("H       %f       %f       %f\n" % (ancestorPile[i,1,0]*convert, ancestorPile[i,1,1]*convert, ancestorPile[i,1,2]*convert))
#	fileXYZ.write("H       %f       %f       %f\n" % (ancestorPile[i,2,0]*convert, ancestorPile[i,2,1]*convert, ancestorPile[i,2,2]*convert))
#	fileXYZ.write("H       %f       %f       %f\n" % (ancestorPile[i,3,0]*convert, ancestorPile[i,3,1]*convert, ancestorPile[i,3,2]*convert))
#	fileXYZ.write("H       %f       %f       %f\n" % (ancestorPile[i,4,0]*convert, ancestorPile[i,4,1]*convert, ancestorPile[i,4,2]*convert))
#	fileXYZ.write("H       %f       %f       %f\n\n" % (ancestorPile[i,5,0]*convert, ancestorPile[i,5,1]*convert, ancestorPile[i,5,2]*convert))
#	b+=1

for item in ancestors: #-----------for jmol xyz file
	fileXYZ.write("5\n")
	fileXYZ.write("CH4 coordinates         ")
	fileXYZ.write("%i\n" % b)
	fileXYZ.write("C       %f       %f       %f\n" % (item[0,0]*convert, item[0,1]*convert, item[0,2]*convert))
	fileXYZ.write("H       %f       %f       %f\n" % (item[1,0]*convert, item[1,1]*convert, item[1,2]*convert))
	fileXYZ.write("H       %f       %f       %f\n" % (item[2,0]*convert, item[2,1]*convert, item[2,2]*convert))
	fileXYZ.write("H       %f       %f       %f\n" % (item[3,0]*convert, item[3,1]*convert, item[3,2]*convert))
	fileXYZ.write("H       %f       %f       %f\n" % (item[4,0]*convert, item[4,1]*convert, item[4,2]*convert))
	b+=1

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
