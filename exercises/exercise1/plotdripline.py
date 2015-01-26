import matplotlib
from numpy import loadtxt
from pylab import *

z,a=loadtxt('dripline.dat',unpack=True,usecols=[0,1])
ez,ea=loadtxt('bedata.dat',unpack=True,usecols=[0,1])

plot(ea-ez,ez,marker='.',linestyle='none',color='r',label='experimental data')
plot(a-z,z,marker='o',label='neutron drip line')
xlabel('Neutron #')
ylabel('Proton #')
title('Neutron drip line')
savefig('neutrondripline.png')
show()
