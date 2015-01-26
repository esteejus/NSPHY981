def brownliquiddrop_cal(a,z, a1=15.49,a2=17.23, a3=0.697, a4=22.6):
    be = a1*a - a2*(a**(2./3.)) - a3*(z**2)/(a**(1./3.))-a4*((a-2.*z)**2)/a
    return be
#find drip line i.e. where b.e.<0
for i in range(1,121):
# might as well start off at n=z and go to some large neutron number like 2*z
    for j in range(1,2*i):
        a_drip=[]
        z_drip=[]
        i_prev=0
        j_prev=0
      #  print i,j
        b_energy=brownliquiddrop_cal(j,i)
        print b_energy
        if b_energy>0.:
            i_prev=i
            j_prev=j
        else:
            a_drip.append(j_prev)
            z_drip.append(i_prev)
            continue 
            

print a_drip
