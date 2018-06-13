import timeit
import math
import itertools
import random
import numpy as np
import scipy
# import matplotlib.pyplot as mpl
import CH5pot
import numpy.linalg
from numpy.linalg import eig, norm 
from scipy.linalg import sqrtm

## - IS THIS FOR MOX??? No: 1, Yes: 2
mox = 1

START = timeit.default_timer()

n0 = 20000 #number of initial walkers
dtau = 10 #time step
m1 = 21874.658
m2 = 1837.152
massarray = [m1,m2,m2,m2,m2,m2]
sigma1 = np.sqrt(dtau / m1)
sigma2 = np.sqrt(dtau / m2)
alpha = 0.5 / dtau
cycles = 10000
timestep = 0


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

def conavg(potentials,trueweights,init):
	potTotal = potentials.sum()
	vref = potTotal / float(len(potentials)) - alpha*(trueweights.sum()-init)/float(init)
	return vref
	
def pcalc(potentials, vref, dtau): #------------to calculate pvalues
	a = -(potentials - vref) * dtau
	pvalues = np.exp(a.astype(float))
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

def conFate(walkers,trueweights):
	minind = np.argmin(trueweights)
	maxind = np.argmax(trueweights)
	if trueweights[minind] <= 1.0/len(walkers):
		walkers = np.concatenate((walkers,np.array([walkers[maxind]])),axis=0)
		walkers = np.delete(walkers,minind,axis=0)
		trueweights = np.append(trueweights,trueweights[maxind]/2)
		trueweights = np.delete(trueweights,minind)
		trueweights[maxind] = trueweights[maxind]/2
	return walkers,trueweights

def COM(walkers): # puts walkers in a COM-shifted [[principal axis]] system
	COMwalkers = []
	
	M = sum(massarray)

	for l in xrange(len(walkers)):
		walker = np.array(walkers[l])
		COMorigin = np.zeros(3)
		#COM shift:
		for i in xrange(0,3):
			for j in xrange(0,6):
				COMorigin[i] += massarray[j]*walker[j][i]/M
		COMshift = [COMorigin,COMorigin,COMorigin,COMorigin,COMorigin,COMorigin]
		walker = walker - COMshift

		# #I-matrix calculation:
		# moi, Bvalue = MomInertia([walker],massarray)

		# w,v = eig(moi[0])
		# princ = [v[0]/norm(v[0]),v[1]/norm(v[1]),v[2]/norm(v[2])]
		# princT = np.transpose(princ)
		# COMwalker = []
		# for i in xrange(len(walker)):
		# 	COMwalker.append(np.dot(walker[i],princT))
		COMwalker = walker
		COMwalkers.append(COMwalker)
	COMwalkers = np.array(COMwalkers)
	return COMwalkers

def Eck(walkers,ref): # must input walkers already in COM frame
	Fvec = np.zeros((3,3))
	Fdot = np.zeros((3,3))
	newwalkers = np.zeros((len(walkers),6,3))
	for l in xrange(len(walkers)):
		x=0
		walker = np.array(walkers[l])
		for i in xrange(0,3):
			for j in xrange(0,6):
				Fvec[i] += massarray[j]*ref[j][i]*walker[j]
		for k in xrange(0,3):
			for m in xrange(0,3):
				Fdot[k,m] = np.dot(Fvec[k],Fvec[m])
		if np.linalg.det(Fdot)==0:
			print "singular matrix!"
			x=1
			Fdotinv = np.linalg.inv(Fdot+perturb)
		else:
			Fdotinv = np.linalg.inv(Fdot)
		Fminushalf = scipy.linalg.sqrtm(Fdotinv)
		if x==1:
			print "perturb success!"
		Eckaxes = np.dot(Fvec,Fminushalf)
		for j in xrange(0,6):
			newwalkers[l][j] = np.dot(walker[j],np.transpose(Eckaxes))
	return newwalkers

def MomInertia(walkers,massarray):
	Bvalues = []
	mois = []
	conv = 1 #5.84946*(10**-34)
	for walker in walkers:
		walker = np.array(walker)
		moi = np.zeros((3,3))
		for j in xrange(0,3):
			for k in xrange(0,3):
				for i in xrange(len(walker)):
					if j==0 and k==0:
						moi[0][0] += massarray[i]*(walker[i][1]**2+walker[i][2]**2)*conv
					if j==1 and k==1:
						moi[1][1] += massarray[i]*(walker[i][0]**2+walker[i][2]**2)*conv
					if j==2 and k==2:
						moi[2][2] += massarray[i]*(walker[i][0]**2+walker[i][1]**2)*conv
					else:
						moi[j][k] += -1*massarray[i]*walker[i][j]*walker[i][k]*conv
		Bvalue = 3/(moi[1][1]+moi[2][2]+moi[0][0]) #(1/moi[1][1]+1/moi[2][2]+1/moi[0][0])/3
		# if timestep>90:
		# 	print Bvalue*220000
		Bvalues.append(Bvalue)

		mois.append(moi)

	return mois, Bvalues

perturb = np.zeros((3,3))
perturb[1][2] = 1
perturb[0][2] = 1
perturb[2][2] = 1


#Initialize Walkers
beginning = np.array([[0.000000000000000, 0.000000000000000, 0.3869923621587414],
[0.000000000000000, 0.000000000000000, -1.810066283748844],
[1.797239666982623, 0.000000000000000,   1.381637275550612],
[-1.797239666982623, 0.000000000000000, 1.381637275550612],
[0.000000000000000, -1.895858229423645, -0.6415748897955779],
[0.000000000000000, 1.895858229423645, -0.6415748897955779]])

c2v = np.array([[0.000000000000000, 0.000000000000000, 0.3869923621587414],
[0.000000000000000, 0.000000000000000, -1.810066283748844],
[1.797239666982623, 0.000000000000000,   1.381637275550612],
[-1.797239666982623, 0.000000000000000, 1.381637275550612],
[0.000000000000000, -1.895858229423645, -0.6415748897955779],
[0.000000000000000, 1.895858229423645, -0.6415748897955779]])


#I-matrix calculation:
moi, Bvalue = MomInertia([c2v],massarray)

w,v = eig(moi[0])
princ = [v[0]/norm(v[0]),v[1]/norm(v[1]),v[2]/norm(v[2])]
princT = np.transpose(princ)
COMwalker = []
for i in xrange(len(c2v)):
	COMwalker.append(np.dot(c2v[i],princT))
c2vREF = COMwalker#COM([c2v])[0]

whereUat = np.zeros((n0,beginning.shape[0],beginning.shape[1]))

for i in range(0,n0):
    whereUat[i] = beginning * 1.1
		
# potentials = V(whereUat)	

#------    FILENAMES   ---------
#######                                          #########
#######                                          #########
#######                                          #########
if mox == 1:
	#fileEnergy = open("5deut_energies.txt", "a+")
	fileXYZ = open("CH5firstrot3.xyz", "a+")
energies = []
population = []
ancestorPile = []


equiltime = 4000
for c in range(0,cycles):    

	if timestep%500 == 0:
		ancestors = whereUat
		whereUfrom = np.arange(0,len(whereUat)) #initialize the "index array", which will reference indexes in "ancestors"
		
	if timestep < equiltime:
		displace = np.random.normal(0,sigma2,whereUat.shape)
		for i in range(0,len(displace)): 
			displace[i,0,:] = np.random.normal(0,sigma1,(1,3))
		whereUat = whereUat + displace 
		potentials = V(whereUat)  #using my function V
		vref = avg(potentials, n0) #umf avg
		pvalues = pcalc(potentials, vref, dtau) #umf pcalc
		whereUat,potentials,whereUfrom = Fate(whereUat,whereUfrom,potentials,pvalues)

	if timestep == equiltime:
		conwalkers = whereUat
		trueweights = np.ones(len(whereUat))

	if timestep >= equiltime:
		displace = np.random.normal(0,sigma2,conwalkers.shape)
		for i in range(0,len(displace)): 
			displace[i,0,:] = np.random.normal(0,sigma1,(1,3))
		conwalkers = conwalkers + displace
		conwalkers = COM(conwalkers)
		conwalkers = Eck(conwalkers,c2vREF)
		mois, Bvalues = np.array(MomInertia(conwalkers,massarray))
		conpotentials = np.array(V(conwalkers)) #+ Bvalues)
		convref = conavg(conpotentials,trueweights,n0)
		energies.append(convref)
		conpvalues = pcalc(conpotentials, convref, dtau)
		trueweights = trueweights * conpvalues
		conwalkers,trueweights = conFate(conwalkers,trueweights)
	
	# Descendant Weighting #	
	# if timestep%500 == 50 and timestep > 4000:
	# 	for i in whereUfrom:
	# 		ancestorPile.append(ancestors[i])

	timestep += 1

convert = 0.52917725
v = np.zeros(len(conwalkers))
v = CH5pot.mycalcpot(conwalkers,len(conwalkers))
b=0
for item in conwalkers: # for jmol xyz file
	if mox == 1:	
		fileXYZ.write("6\n")
		fileXYZ.write("weight/potential/Bvalue 		%f         %f         %f" % (trueweights[b],v[b],Bvalues[b]))
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
		print "CH5+ coordinates         %f         %f" % (v[b],Bvalues[b])
		print "%i\n" % b
		print "C       %f       %f       %f\n" % (item[0,0]*convert, item[0,1]*convert, item[0,2]*convert)
		print "H       %f       %f       %f\n" % (item[1,0]*convert, item[1,1]*convert, item[1,2]*convert)
		print "H       %f       %f       %f\n" % (item[2,0]*convert, item[2,1]*convert, item[2,2]*convert)
		print "H       %f       %f       %f\n" % (item[3,0]*convert, item[3,1]*convert, item[3,2]*convert)
		print "H       %f       %f       %f\n" % (item[4,0]*convert, item[4,1]*convert, item[4,2]*convert)
		print "H       %f       %f       %f\n\n" % (item[5,0]*convert, item[5,1]*convert, item[5,2]*convert)
		b+=1

if mox == 1:
	#fileEnergy.close()
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
