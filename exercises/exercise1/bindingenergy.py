import matplotlib
from numpy import loadtxt
from pylab import *


z,a,be=loadtxt('bedata.dat',unpack=True,usecols=[0,1,2])

#plot Binding Energies 
plot(a,be/a,linestyle='none',marker='o')

xlabel('Nucleons (A)')
ylabel('B.E. (MeV/A)')
title('Binding Energy of Nuclei')
savefig("expBE.png")
show()

#calculate separation energy given a nuceli
def separationenergy(neutron=[],nuceli=[],protonnum=0):
    prevbe=0
    for x,y,w in zip(z,a,be):
        if x==protonnum:
            nuceli.append(w-prevbe)#1n separation energy  
            neutron.append(y) #neutron number 
            prevbe=w #set previous binding energy
 
#remove first enerty which does not represent a separation energy
    neutron.pop(0)
    nuceli.pop(0)
    return


#Create plots for Oxygen, Calcium, Nickel, Tin, and lead
oxygen=[]
oxygenneutron=[]
calcium=[]
calciumneutron=[]
nickel=[]
nickelneutron=[]
tin=[]
tinneutron=[]
lead=[]
leadneutron=[]

separationenergy(oxygenneutron,oxygen,8)
plot(oxygenneutron,oxygen)
xlabel("Neutrons")
ylabel('S_n (MeV)')
title('Separation Energy of Oxygen')
savefig("separationenergy_oxygen.png")
show()
    
    
separationenergy(calciumneutron,calcium,20)
plot(calciumneutron,calcium)
xlabel("Neutrons")
ylabel('S_n (MeV)')
title('Separation Energy of Calcium')
savefig("separationenergy_calcium.png")
show()
    
separationenergy(nickelneutron,nickel,28)
plot(nickelneutron,nickel)
xlabel("Neutrons")
ylabel('S_n (MeV)')
title('Separation Energy of Nickel')
savefig("separationenergy_nickel.png")
show()

separationenergy(tinneutron,tin,50)
plot(tinneutron,tin)
xlabel("Neutrons")
ylabel('S_n (MeV)')
title('Separation Energy of Tin')
savefig("separationenergy_tin.png")
show()

separationenergy(leadneutron,lead,82)
plot(leadneutron,lead)
xlabel("Neutrons")
ylabel('S_n (MeV)')
title('Separation Energy of Lead')
savefig("separationenergy_lead.png")
show()
