// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __BAT__CCXSEC__H
#define __BAT__CCXSEC__H

#include <BAT/BCModel.h>

// This is a sdmodel0 header file.
// Model source code is located in file sdmodel0Analysis/sdmodel0.cxx

// ---------------------------------------------------------
class sdmodel0 : public BCModel
{
   public:

      // Constructors and destructor
      sdmodel0();
      sdmodel0(const char * name);
      ~sdmodel0();

      // Methods to overload, see file sdmodel0.cxx
      void DefineParameters();
      double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector<double> &parameters);
      // void MCMCIterationInterface();
};
// ---------------------------------------------------------

#endif

