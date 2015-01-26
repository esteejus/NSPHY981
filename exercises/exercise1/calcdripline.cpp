#include <iostream>
#include <cmath>
using namespace std;

double calcdripline(int,int);

int main(){
  bool flag=false;
  double deltae=0;
  int i_prev=0,j_prev=0;
  for(int i=1;i<=120;i++){
    flag=false;
    for(int j=1;j<600;j++){
      // cout<<i<<" "<<j<<endl;
      //cout<<be<<endl;
      deltae=calcdripline(i,j+1)-calcdripline(i,j);
      if(deltae<0.){
	cout<<deltae<<endl;
	  cout<<i<<" "<<j-1<<endl;
	  break;
	}
      
    }
  }
  


  
  
  return 0;
}

double calcdripline(int z, int a){
  double be=0;//binding energy
  double a1=15.49,a2=17.23, a3=0.697, a4=22.6;
  be=a1*a - a2*pow(a,(2./3.))-a3*(pow(z,2))/(pow(a,(1./3.)))-a4*pow((a-2.*z),2.)/a;

    return be;
}
