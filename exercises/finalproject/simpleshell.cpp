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
  
  /*
  mat a(2,2);
  a(0,0)=1;
  a(0,1)=2;
  a(1,0)=3;
  a(1,1)=4;
  */

  //a.print();
   
 int npart=4;
 int splevels=6;
 int  nsd=3;   
 mat hamilt(3,3);
  hamilt.zeros();
  mscheme(npart,0,splevels);
  for(int i=0;i<nsd;i++){
       rowvec wf = a.row(i);//wave function(SD) for bra
              onebody(wf,a,nsd,hamilt,i);
       twobody(wf,a,nsd,splevels,hamilt,i);

  }  
    hamilt.print();

  mat v = linspace<mat>(1,4,4);
  /*
  for(int k=0;k<nsd;k++){
  for(int i=0;i<splevels;i++){
    for(int j=0;j<4;j++){
     
      if(v(4-j)<splevels) v(4-j)+=1;

    }
  }
  }
  v.print();
  */ 
  return 0 ;
  
}


void mscheme(int npart,int mvalue,int splevels){ 
  int nsd=nChoosek(splevels,npart);
   mat tempsd(nsd,npart);//all slater determinents 
   int pos = npart-1;//position for switching particle positions
   tempsd.zeros(); 
   bool resum = false;
   for(int i=0;i<nsd;i++){
     for(int j=0;j<npart;j++){      
       

       if(i==0)tempsd(i,j)=j+1;//initial filling
       
       else{
       if(j!=pos)tempsd(i,j)=tempsd(i-1,j);
       else if(j>pos && resum==true){
	 tempsd(i,j)=tempsd(i,j-1)+1;
	 if(j==npart){
	   resum=false;
	   pos=npart-1;
	 }
       }
	 else{
	 if(tempsd(i-1,j)<splevels)tempsd(i,j)=tempsd(i-1,j)+1;
	 else{
	   if(tempsd(i-1,j-1)+1!=tempsd(i-1,j)){
	     tempsd(i,j-1)=tempsd(i-1,j-1)+1;
	     tempsd(i,j)=tempsd(i,j-1)+1;
	   }
	   else{
	     resum =true;
	        pos=pos-2;
	     	     j=j-2;
	   }
	 }

       }//end of first else	 
        if(pos==0 && j==0 && tempsd(i-1,j)<splevels){
	 tempsd(i,j)=tempsd(i-1,j)+1;
       }
       }//end of nsd else
       

   
     }//end of j
     //     tempsd.print();
   }//end of nsd




   return;
}

void onebody(rowvec wf,mat sd, int nsd, mat &hamilt, int &i){
  //int matxel=0; 
 
  for(int j=0;j<nsd;j++){
    rowvec wf2 = sd.row(j);//bra vector
    if(equal(wf,wf2)==true){//bra and ket test
      for(int k=0;k<wf.n_elem-1;k++){//sum over sp states elem of vector
	if((int)wf(k+1)%2==0){
	  //here i count 2,4,6,8,... even levels
	  //from them i evauluate p; p=level/2; as above 1,2,3,....
	  //take 2*(p-1) where as_scalar(wf(k+1)) is p
	  //the 2* is for both particles since 
	  //i am only counting the even labels
	  //	  hamilt(i,j)=hamilt(i,j)+2*(as_scalar(wf(k+1)/2-1));
 
	}
      }
    }
  }
 
  return ;
  }


void twobody(rowvec wf,mat sd, int nsd,int splevels, mat &hamilt, int &i){ 
  int me=-2;//matrix element
  
rowvec temp(wf.n_elem);//temp wf for calculating matrix elements
    
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
		if(equal(temp,wf2)==true)hamilt(i,j)+=me;

		//		cout<<"After"<<endl;
		//temp.print();
		//		hamilt.print();
		//cout<<endl<<endl;
	      }
	    }
	  }	  
	}
      }    
    }
  return;
}

  bool equal(rowvec wf, rowvec wf2){
    bool flag=true;
    if(wf.n_elem!=wf2.n_elem)cout<<"Num elements not equal"<<endl; 
    for(int i=0;i<wf.n_elem;i++){
      if(wf(i)!=wf2(i)) flag=false;
    }
    return flag;
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
