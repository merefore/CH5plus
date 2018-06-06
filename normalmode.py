import numpy as np
import CH5pot
import sys
from numpy.lib.scimath import sqrt as csqrt
from numpy.linalg import eig
from numpy.linalg import inv

p = 1
N = 6
delta = 0.001
mc = 21874.65821
mh = 1837.15234
mass = [mc, mh, mh, mh, mh, mh]

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

print "Eigenfrequencies (cm-1)"
for i in xrange(0,len(wcm)):
	print wcm[i]

print "First derivatives"
for i in xrange(0,len(Ffirst)):
	print '%f' % Ffirst[i]

#print format(w, 'f')
#print np.dot(inv(v),np.dot(F,v))


