#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>

using namespace std;

void newmatrixelements();

int main(){
  int stat=-999;//status of the script
  stat=system("./test.sh");//execute script 
  if(stat==0){
    //script finished successfully 
    
    
  }
  else{
    cout<<"Script failed"<<endl;
      cout<<"System status "<<stat<<endl;
    
  }
  //newmatrixelements();
  









  return 0;

}


void newmatrixelements(){
ifstream hamilton;
int nlines=7;
string temp;

 int count=0;

 hamilton.open("/home/justin/aaa/nushellx/sps/usdb.int");
if(!hamilton.is_open())cout<<".int file not found"<<endl;
for(int i=0;i<nlines;i++)getline(hamilton,temp);

while(!hamilton.eof()){
int a=-1,b=0,c=0,d=0,j=0,t=0;//abcd are the states j total J and t isospin
 double me=0.0;//matrix element
   hamilton>>a>>b>>c>>d>>j>>t>>me;
   if(a!=-1){//incase there is a blank line at eof prevents double read in
     //   cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<j<<" "<<t<<" "<<me<<endl;
     count++;}
 }
 
// cout<<count<<"total me"<<endl;
 return;
}
