import matplotlib
from numpy import loadtxt
from pylab import *


z,a,be=loadtxt('bedata.dat',unpack=True,usecols=[0,1,2])

#plot Binding Energies 
plot(a,be/a,linestyle='none',marker='o',label='Exp. Data')

xlabel('Nucleons (A)')
ylabel('B.E. (MeV/A)')
title('Binding Energy of Nuclei')
savefig("expBE.png")
show()


def brownliquiddrop(a,z,bindingenergy=[], a1=15.49,a2=17.23, a3=0.697, a4=22.6):

    for i in range(0,len(a)):
        bindingenergy.append(a1*a[i] - a2*(a[i]**(2./3.)) - a3*(z[i]**2)/(a[i]**(1./3.))-a4*((a[i]-2.*z[i])**2)/a[i])
        #if a[i]==124:
        #   print a[i],z[i],bindingenergy[i],a1*a[i],a2*(a[i]**(2./3.)),a3*(z[i]**2)/(a[i]**(1./3.)),a4*((a[i]-2.*z[i])**2)/a[i]
    #print bindingenergy
    return

brownbe=[]
brownliquiddrop(a,z,brownbe,15.49,0.,0.,0.)
plot(a,brownbe/a,linestyle='none',marker='x',alpha=.6, color='g',label='Volume ON')
brownbe=[]
brownliquiddrop(a,z,brownbe,15.49,17.23,0.,0.)
plot(a,brownbe/a,linestyle='none',marker='^',alpha=.5, color='y',label='Surfuce ON')
brownbe=[]
brownliquiddrop(a,z,brownbe,15.49,17.23,0.697,0.)
plot(a,brownbe/a,linestyle='none',marker='s',alpha=.5, color='c',label='Coulomb ON')
#plot Binding Energies 
plot(a,be/a,linestyle='none',marker='o',label='Exp. Data')
brownbe=[]
brownliquiddrop(a,z,brownbe)
plot(a,brownbe/a,linestyle='none',marker='.',alpha=.5, color='r',label='A.Brown Liquid Drop')
legend(loc='upper right')
savefig("Liquiddrop.png")
show()


fluorinediff=[]
fluorinea=[]

for d,g,e,t in zip(z,a,be,brownbe):
    if d==9.: 
        fluorinediff.append(e-t)
        fluorinea.append(g)

xlabel('Nucleons (A)')
ylabel('BE(exp)-BE(liquid drop) MeV')
title('Fluorine(Z=9) Theoretical residuals in Binding Energy')
plot(fluorinea,fluorinediff,marker='o', color='b')
savefig("fluorineresiduls.png")
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

