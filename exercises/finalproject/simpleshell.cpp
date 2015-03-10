#include <iostream>
#include <armadillo>
#include <stdio.h>

using namespace std;
using namespace arma;

void mscheme(int,int &,int, mat &);
unsigned nChoosek(unsigned,unsigned);
bool equal(rowvec,rowvec);
void onebody(rowvec,mat,int,mat&,int&);
void twobody(rowvec,mat,int,int,mat&,int&);
void printCombination(int*, int, int, int);
void combinationUtil(int*, int*, int, int,int, int, int, rowvec &sd,int &counter,mat &temp);

int main(){
     
 int npart=4;
 int splevels=6;
 int nsd=nChoosek(splevels,npart);
 mat sd(nsd,npart);//all possible slater det.
 mscheme(npart,nsd,splevels,sd);//sorts only pairs
 mat hamilt(nsd,nsd);
 hamilt.zeros();
 for(int i=0;i<nsd;i++){
   rowvec wf = sd.row(i);//wave function(SD) for bra
   onebody(wf,sd,nsd,hamilt,i);
   twobody(wf,sd,nsd,splevels,hamilt,i);
   
 }  
  hamilt.print();
  

  return 0 ;
  
}


void mscheme(int npart,int &nsd,int splevels, mat &tempsd){
  int arr[splevels];  
  for(int k=1;k<=splevels;k++)arr[k-1]=k;  
  int z=0;
  int data[splevels];
  nsd=nChoosek(splevels,npart);
  mat tempsd1(nsd,npart);//all possible sd
  rowvec sd(npart);
  //run once to get sd number , z here.
  combinationUtil(arr,data,0,splevels-1,0,npart,nsd,sd,z,tempsd1); 
  nsd=z;//set number of sd to total counted states availible with pairing
   z=0;//set counter back to 0 
 tempsd.zeros(); 
 combinationUtil(arr,data,0,splevels-1,0,npart,nsd,sd,z,tempsd); 
 //sd.print(); 
 // tempsd.print();
 

} 


   /* arr[]  ---> Input Array
   data[] ---> Temporary array to store current combination
   start & end ---> Staring and Ending indexes in arr[]
   index  ---> Current index in data[]
   r ---> Size of a combination to be printed */
void combinationUtil(int arr[], int data[], int start, int end, int index, int r, int nsd, rowvec &sd,int &counter,mat &tempwf)
   {
     // Current combination is ready to be printed, print it
     if(index==r){
	 if (index == r)
	   {
	     
	     
	     for (int j=0; j<r; j++){
	       
	       //if((j+1)%2==0 && data[j]==data[j-1]+1 && data[j]%2==0){
		 sd(j)=data[j];
		 //sd(j-1)=data[j-1];
		 //}	 // cout<<"counter"<<counter<<endl;
		 // }
	       //else if((j+1)%2 && (data[j]!=data[j-1]+1)){
		  //if(data[j]==2 && data[j-1]==1)cout<<j<<endl;
		  //sd.zeros();
		 }
	       int product=1;
	       for(int l=0;l<r;l++){
		 product=sd(l)*product;
	       }
	       //	 if(product!=0){
	       //  counter++;
	       //	   tempwf.insert_rows(counter-1,sd);
		   //	   sd.print();
		   //}
	       
	   
	     
	     //printf("\n");
	     sd.print();
	     return;
	   }
   }
     
     // replace index with all possible elements. The condition
     // "end-i+1 >= r-index" makes sure that including one element
     // at index will make a combination with remaining elements
     // at remaining positions
     for (int i=start; i<=end && end-i+1 >= r-index; i++)
       {
	 data[index] = arr[i];
	 combinationUtil(arr, data, i+1, end, index+1, r,nsd,sd,counter,tempwf);
       }
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
