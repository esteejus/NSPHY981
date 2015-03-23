
// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "sdmodel0.h"

#include <BAT/BCMath.h>
#include <cmath>
#include<BAT/BCDataSet.h>
#include<BAT/BCDataPoint.h>
#include<BAT/BCParameter.h>

// ---------------------------------------------------------
sdmodel0::sdmodel0() : BCModel()
{
  // default constructor
   DefineParameters();
};

// ---------------------------------------------------------
sdmodel0::sdmodel0(const char * name) : BCModel(name)
{
   // constructor
   DefineParameters();
};

// ---------------------------------------------------------
sdmodel0::~sdmodel0()
   // default destructor
{
};

// ---------------------------------------------------------
void sdmodel0::DefineParameters()
{
   // Add parameters to your model here.
   // You can then use them in the methods below by calling the
   // parameters.at(i) or parameters[i], where i is the index
   // of the parameter. The indices increase from 0 according to the
   // order of adding the parameters.

//   AddParameter("x", -10.0, 10.0); // index 0
//   AddParameter("y",  -5.0,  5.0); // index 1
  //1= s1/2  3= d3/2 5= d5/2
  //since i have a closed core 16o I'm only considering valence neutrons
  //thus the pauli principle slelects these "good" J
  //if i opened the core i could look at good combinations of T isospin
  //i then the triangle rule would apply and i would have to combine good T
  AddParameter("a",-100,100);
  AddParameter("b",-100,100);
  AddParameter("c",-100,100);
  //SetPriorGauss(0,2,3);
  //  SetPriorGauss(1,3,5);
  //SetPriorGauss(3,5,7);
  
}

// ---------------------------------------------------------
double sdmodel0::LogLikelihood(const std::vector<double> &parameters)
{
   // This methods returns the logarithm of the conditional probability
   // p(data|parameters). This is where you have to define your model.

   double logprob = 0.;

//   double x = parameters.at(0);
//   double y = parameters.at(1);
//   double eps = 0.5;

   // Breit-Wigner distribution of x with nuisance parameter y
//   logprob += BCMath::LogBreitWignerNonRel(x + eps*y, 0.0, 1.0);

   for (int i = 0; i < GetNDataPoints(); ++i)
     {
       // get data point
       // double j    = GetDataPoint(i)->GetValue(0);//total J
       double speExp   = GetDataPoint(i)->GetValue(0);//single particle(sp) Experimental(Exp)energy
       double errspeExp = GetDataPoint(i)->GetValue(1);//sp exp error
       double speTh=0;
       // calculate expected single part energy
       for(int i=0;i<3;i++){
	 speTh+= parameters.at(i)*i;
       }

       // calculate likelihood
       logprob += BCMath::LogGaus(speExp, speTh, errspeExp);
     }
   return logprob;
}


// ---------------------------------------------------------
double sdmodel0::LogAPrioriProbability(const std::vector<double> &parameters)
{
   // This method returns the logarithm of the prior probability for the
   // parameters p(parameters).

   double logprob = 0.;

   // double x = parameters.at(0);
   double y = parameters.at(1);
   double dz = parameters.at(0);
   double dx = parameters.at(2);
   
   logprob += BCMath::LogGaus(dx,0,5);                    // flat prior for x
   logprob += BCMath::LogGaus(dz,3,5);
   logprob += BCMath::LogGaus(y, 1.,5. );   // Gaussian prior for y

   return logprob;
}
// ---------------------------------------------------------


