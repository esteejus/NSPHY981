#include <iostream>
#include <armadillo>


using namespace std;
using namespace arma;

void mscheme(int,int,int);
unsigned nChoosek(unsigned,unsigned);
bool equal(rowvec,rowvec);
void onebody(rowvec,mat,int,mat&,int&);
void twobody(rowvec,mat,int,int,mat&,int&);

int main(){
  
  /*
  mat a (3,4);
  a(0,0)=1;
  a(0,1)=2;
  a(0,2)=3;
  a(0,3)=4;

  a(1,0)=1;
  a(1,1)=2;
  a(1,2)=5;
  a(1,3)=6;

  a(2,0)=3;
  a(2,1)=4;
  a(2,2)=5;
  a(2,3)=6;
  */
  mat a(2,2);
  a(0,0)=1;
  a(0,1)=2;
  a(1,0)=3;
  a(1,1)=4;
  
  //a.print();
   
 int npart=4;
 int splevels=6;
 int  nsd=2;   
 mat hamilt(2,2);
  hamilt.zeros();
mscheme(npart,0,splevels);
  for(int i=0;i<nsd;i++){
    //    rowvec wf = a.row(i);//wave function(SD) for bra
    //    onebody(wf,a,nsd,hamilt,i);
    //    twobody(wf,a,nsd,splevels,hamilt,i);

  }  
  
  return 0 ;
  
}


void mscheme(int npart,int mvalue,int splevels){ 
  int nsd=nChoosek(splevels,npart);
   mat tempsd(nsd,npart);//all slater determinents 
   int pos = npart-1;//position for switching particle positions
   tempsd.zeros(); 
   
   for(int i=0;i<nsd;i++){
     for(int j=0;j<npart;j++){      
       

       if(i==0)tempsd(i,j)=j+1;//initial filling
       
       else{
       if(j!=pos && pos!=0)tempsd(i,j)=tempsd(i-1,j);
       else{
	 if(tempsd(i-1,j)<splevels)tempsd(i,j)=tempsd(i-1,j)+1;
	 else{
	   if(tempsd(i-1,j-1)+1!=tempsd(i-1,j)){
	     tempsd(i,j-1)=tempsd(i-1,j-1)+1;
	     tempsd(i,j)=tempsd(i,j-1)+1;
	   }
	   else{
	 
	        pos=pos-1;
	     	     j=j-1;
	   }
	 }

       }//end of first else	 
        if(pos==0 && j==0 && tempsd(i-1,j)<splevels){
	 tempsd(i,j)=tempsd(i-1,j)+1;
       }
       }//end of nsd else
       

   
     }//end of j
     tempsd.print();
   }//end of nsd




   return;
}

void onebody(rowvec wf,mat sd, int nsd, mat &hamilt, int &i){
  int matxel=0; 
 
  for(int j=0;j<nsd;j++){
    rowvec wf2 = sd.row(j);
    if(equal(wf,wf2)==true)hamilt(i,j)=as_scalar(sum(wf));     
      }
  //  hamilt.print();
  return ;
}


void twobody(rowvec wf,mat sd, int nsd,int splevels, mat &hamilt, int &i){ 
  rowvec temp(4);//temp wf for calculating matrix elements
    
    for(int j=0;j<nsd;j++){
      rowvec wf2=sd.row(j);
      
      //a+(p+)a+(p-)a(q+)a(q-)
      for(int r=0;r<wf.n_elem;r++){ //loop over where index is   
	for(int p=1;p<=splevels;p++){//loop over pair creation
	  for(int q=1;q<=splevels;q++){//loop over pair anilhiation
	    if(r%2==0 && !(p%2==0) && !(q%2==0)){//need to look at pairs 
	      //	      cout<<"Before"<<" "<<p<<p+1<<q<<q+1<<endl;
	      temp=wf;
	      //temp.print();
	      //checking to see if pair state exists to be anihilated 
	      //i.e. if temp=|1234>
	      //|1234>-|1200>=|0034>
	      if(temp(r)-q==0 && temp(r+1)-(q+1)==0){
		//next i create a state with p,p+1
		//i.e. |p(p+1)34>
		//if p and p+1 is 34 resulting in a state |3434>=0
		//since the program will not be able to find a SD with |3434> by construction, then the matrix element will be zero as given element wise by the equal function
		temp(r)=p;
		temp(r+1)=p+1;
		if(equal(temp,wf2)==true)hamilt(i,j)+=2;

		//		cout<<"After"<<endl;
		//temp.print();
		hamilt.print();
		cout<<endl<<endl;
	      }
	    }
	  }	  
	}
      }    
    }
  return;
}

  bool equal(rowvec wf, rowvec wf2){
    if(wf.n_elem!=wf2.n_elem)cout<<"Num elements not equal"<<endl; 
    for(int i=0;i<wf.n_elem;i++)
      if(wf(i)==wf2(i))return true;
      else return false;

  }

//barrowed from online
unsigned nChoosek( unsigned n, unsigned k )
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  int result = n;
  for( int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}
