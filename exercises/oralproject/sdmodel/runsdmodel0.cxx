// ***************************************************************
// This file was created using the CreateProject.sh script
// for project sdmodel0Analysis.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCDataSet.h>
#include<BAT/BCEngineMCMC.h>

#include "sdmodel0.h"

int main()
{

  // BCEngineMCMC::MCMCSetFlagFillHistograms(false);

  // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   
   // read data set from file
   BCDataSet * mydataset = new BCDataSet();
   mydataset->ReadDataFromFile("data.dat",3);
  
   
   
   // create new sdmodel0 object
   sdmodel0 *m = new sdmodel0();
   m->SetDataSet(mydataset);
   

   //      BCEngineMCMC::MCMCSetFlagFillHistograms(false);
   
   BCLog::OutSummary("Test model created");
   
   // create a new summary tool object
   BCSummaryTool * summary = new BCSummaryTool(m);
   
   
   // perform your analysis here

   // normalize the posterior, i.e. integrate posterior
   // over the full parameter space
     m->Normalize();

   // run MCMC and marginalize posterior wrt. all parameters
   // and all combinations of two parameters
     m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
     m->MCMCSetPrecision(BCEngineMCMC::kMedium);     
     m->MarginalizeAll();

   // run mode finding; by default using Minuit
  //  m->FindMode();

   // if MCMC was run before (MarginalizeAll()) it is
   // possible to use the mode found by MCMC as
   // starting point of Minuit minimization
   //m->FindMode( m->GetBestFitParameters() );

   // draw all marginalized distributions into a PostScript file
    m->PrintAllMarginalized("/home/justin/Documents/sdmodel0_plots.ps");

   // print all summary plots
//   summary->PrintParameterPlot("sdmodel0_parameters.eps");
//   summary->PrintCorrelationPlot("sdmodel0_correlation.eps");
//   summary->PrintKnowledgeUpdatePlots("sdmodel0_update.ps");

   // calculate p-value
//   m->CalculatePValue( m->GetBestFitParameters() );

   // print results of the analysis into a text file
   //   m->PrintResults("sdmodel0_results.txt");

   delete m;
   delete summary;

   BCLog::OutSummary("Test program ran successfully");
   BCLog::OutSummary("Exiting");

   // close log file
   BCLog::CloseLog();

   return 0;

}

