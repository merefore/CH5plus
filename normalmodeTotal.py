import numpy as np
import CH5pot
import sys
import timeit
from numpy.lib.scimath import sqrt as csqrt
from numpy.linalg import eig
from numpy.linalg import inv
from numpy.random import uniform as unirand

START = timeit.default_timer()

#p = int(sys.argv[1])
p = 2
N = 6
delta = 0.001
mc = 21874.65821
mh = 1837.15234
md = 3671.48233
mass = [mc, mh, mh, mh, mh, mh]
#mass = [mc, md, md, md, md, md]
for i in xrange(1,len(sys.argv)):
	deutindex = int(sys.argv[i])
	#mass[deutindex] = md
	mass[deutindex] = mh

longmass = [mass[0], mass[0], mass[0], mass[1], mass[1], mass[1], mass[2], mass[2], mass[2], mass[3], mass[3], mass[3], mass[4], mass[4], mass[4], mass[5], mass[5], mass[5]]

filename = "CH5_C2V_imagreal.xyz"
#filename = "CD1H4_CsSP_%i.xyz" % int(sys.argv[1])
#filename = "CD3H2_CsMIN_%i%i%i.xyz" % (int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))

def shapen(matrix):
	walker = np.zeros((6,3))
	rowcount = 0
	for i in xrange(0,18):
		ind = i%3
		walker[rowcount][ind]=matrix[i]
		if ind == 2:
			rowcount+=1
	return walker


#######      Cs MIN structure
if p==1:
	beg = np.array([[0.000000000000000, 0.000000000000000, 0.000000000000000],
[0.1318851447521099, 2.088940054609643, 0.000000000000000],
[1.786540362044548, -1.386051328559878, 0.000000000000000],
[2.233806981137821, 0.3567096955165336, 0.000000000000000],
[-0.8247121421923925, -0.6295306113384560, -1.775332267901544],
[-0.8247121421923925, -0.6295306113384560, 1.775332267901544]])


#######      C2V structure
if p==2:
	beg = np.array([[0.000000000000000, 0.000000000000000, 0.3869923621587414],
[0.000000000000000, 0.000000000000000, -1.810066283748844],
[1.797239666982623, 0.000000000000000,   1.381637275550612],
[-1.797239666982623, 0.000000000000000, 1.381637275550612],
[0.000000000000000, -1.895858229423645, -0.6415748897955779],
[0.000000000000000, 1.895858229423645, -0.6415748897955779]])


#######      Cs SP structure
if p==3:
	beg = np.array([[0.000000000000000, 0.000000000000000, 0.000000000000000],
[1.931652478009080, -4.5126502395556294E-008,  -0.6830921182334913],
[5.4640011799588715E-017, 0.8923685824271653, 2.083855680290835],
[-5.4640011799588715E-017, -0.8923685824271653, 2.083855680290835],
[-1.145620108130841, -1.659539840225091, -0.4971351597887673],
[-1.145620108130841, 1.659539840225091, -0.4971351597887673]])

####### First calculation of F matrix
F = np.zeros((18,18))
Ffirst = np.zeros(18)
row=0
col=0
track = 0

for z in xrange(0,N):
	for y in xrange(0,3):
		shift = np.zeros((6,3))
		shift[z][y] = delta
		shiftarr = [beg+shift, beg-shift]
		Vfirst = CH5pot.mycalcpot(shiftarr,len(shiftarr))
		Ffirst[track] = (1/(2*delta))*(Vfirst[0]-Vfirst[1])
		track += 1
		

for i in xrange(0,N):
	for j in xrange(0,3):

		shift1 = np.zeros((6,3))
		shift1[i][j] = delta

		for k in xrange(0,N):
			for m in xrange(0,3):

				shift2 = np.zeros((6,3))
				shift2[k][m] = delta
				shiftarray = [beg+shift1+shift2, beg+shift1-shift2, beg-shift1+shift2, beg-shift1-shift2]
				V = CH5pot.mycalcpot(shiftarray,len(shiftarray)) 
				#print V
				F[row][col]=(1/(4*delta*delta))*(V[0]-V[1]-V[2]+V[3])/(np.sqrt(mass[i])*np.sqrt(mass[k]))
				col+=1

		col=0
		row+=1


w,v = eig(F)
w = np.array(w)
wsq = csqrt(w)
wcm = wsq*219474.6567

Ybeg = np.zeros((6,3))
for i in xrange(0,N):
	Ybeg[i][0] = beg[i][0]*np.sqrt(mass[i])
	Ybeg[i][1] = beg[i][1]*np.sqrt(mass[i])
	Ybeg[i][2] = beg[i][2]*np.sqrt(mass[i])

Ys = Ybeg.flatten()
Qbeg = np.dot(inv(v),Ys)

qone = np.zeros(18)

for n in xrange(0,18):
	add = 0
	for i in xrange(0,18):
		add+=v[i][n]*Ffirst[i]/(np.sqrt(longmass[i]))
	qone[n] = add

diffq = -qone/w

diffy = np.dot(v,diffq)

diffx = diffy/np.sqrt(longmass)

diffx = shapen(diffx)

newbeg = beg + diffx

####### Now use newbeg (shifted to minimum) to find better F matrix

track = 0
newfirst = np.zeros(18)
for z in xrange(0,N):
	for y in xrange(0,3):
		shift = np.zeros((6,3))
		shift[z][y] = delta
		shiftarr = [newbeg+shift, newbeg-shift]
		Vfirst = CH5pot.mycalcpot(shiftarr,len(shiftarr))
		newfirst[track] = (1/(2*delta))*(Vfirst[0]-Vfirst[1])
		track += 1

newF = np.zeros((18,18))
row=0
col=0
for i in xrange(0,N):
	for j in xrange(0,3):
	
		shift1 = np.zeros((6,3))
		shift1[i][j] = delta

		for k in xrange(0,N):
			for m in xrange(0,3):

				shift2 = np.zeros((6,3))
				shift2[k][m] = delta
				shiftarray = [newbeg+shift1+shift2, newbeg+shift1-shift2, newbeg-shift1+shift2, newbeg-shift1-shift2]
				V = CH5pot.mycalcpot(shiftarray,len(shiftarray)) 
				newF[row][col]=(1/(4*delta*delta))*(V[0]-V[1]-V[2]+V[3])/(np.sqrt(mass[i])*np.sqrt(mass[k]))
				col+=1

		col=0
		row+=1


nw,nv = eig(newF)
nw = np.array(nw)
nwsq = csqrt(nw)
#################   **IF YOU WANT TO USE THE MAGNITUDE OF THE IMAGINARY FREQUENCIES**    ######################
for h in xrange(len(nwsq)):
	if np.real(nwsq[h])==0 and np.imag(nwsq[h])!=0:
		nwsq[h] = np.imag(nwsq[h])
################                           ####
nwsq = nwsq.real
nwcm = nwsq*219474.6567


####### Now we do MC sampling!

lim = []
for i in xrange(0,len(nwsq)):
	if nwsq[i] < 0.001:
		nwsq[i] = 0
		twosig = 1
	else:
		twosig = 3/(np.sqrt(2*nwsq[i]))
	lim.append(twosig)

def P(qs):
	product = 1
	for i in xrange(0,len(qs)):
		phi = np.exp(-nwsq[i]*(qs[i]**2))
		product = product * phi
	return product

totaltime = 1000000000 
keepers = []

for t in xrange(0,totaltime):
	displace = np.zeros(18)
	for i in xrange(0,len(displace)):
		displace[i] = unirand(-lim[i],lim[i])	
	prob = P(displace)
	key = unirand()
	if prob>key:
		keepers.append(displace)

output = open(filename, "a+")
b=0
convert = 0.52917725
for i in xrange(0,len(keepers)):
	keepy = np.dot(v,keepers[i])
	keepx = keepy/np.sqrt(longmass)
	keepx = shapen(keepx)
	item = newbeg + keepx
	output.write("6\n")
	output.write("CH5+ coordinates         ")
	output.write("%i\n" % b)
	output.write("C       %f       %f       %f\n" % (item[0,0]*convert, item[0,1]*convert, item[0,2]*convert))
	output.write("D       %f       %f       %f\n" % (item[1,0]*convert, item[1,1]*convert, item[1,2]*convert))
	output.write("D       %f       %f       %f\n" % (item[2,0]*convert, item[2,1]*convert, item[2,2]*convert))
	output.write("D       %f       %f       %f\n" % (item[3,0]*convert, item[3,1]*convert, item[3,2]*convert))
	output.write("H       %f       %f       %f\n" % (item[4,0]*convert, item[4,1]*convert, item[4,2]*convert))
	output.write("H       %f       %f       %f\n\n" % (item[5,0]*convert, item[5,1]*convert, item[5,2]*convert))
	b+=1

output.close()

print len(keepers)

STOP = timeit.default_timer()
print "Time to calculate: ", STOP - START, " seconds"


