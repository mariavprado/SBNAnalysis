#include <stdio.h>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMinuit.h"
#include <errno.h>
#include <limits.h>
#include "TRandom3.h"
#include "TMath.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TTree.h"
#include <math.h>
#include "TLegend.h"
#include "TH3.h"
#include "TGraph2D.h"
#include <fstream>
#include "TF1.h"
#include <stdlib.h>

using namespace std;

//The specific arrangement of this function goes with the TF1 class. *x and *parameters are 1D arrays.
double FunctionOutL( double *x, double *parameters ){

    double L=x[0];
    
    double prob = 0.1*pow(sin(parameters[0]*L/(3*pow(10,8)*4*parameters[1]*6.582119514*pow(10,-16))),2)*100;
    
    return prob;


}

double FunctionOutE( double *x, double *parameters ){
    
    double E=x[0];
    
    double prob = 0.1*pow(sin(parameters[0]*parameters[1]/(3*pow(10,8)*4*E*pow(10,9)*6.582119514*pow(10,-16))),2)*100;
    
    return prob;
    
    
}

void GraphGenerator(double DeltaM41, std::string name,std::string name2,std::string name3,std::string name4){

    TCanvas* canvas1 = new TCanvas(name.c_str(),"3+1 Neutrino Model",200,10,600,900);
    
    TPad* pad1 = new TPad("pad1","Ev = 700 MeV",0.05,0.50,0.95,0.95);
    TPad* pad2 = new TPad("pad2","L = 600 m",0.05,0.05,0.95,0.50);
    pad1->Draw();
    pad2->Draw();
    
    pad1->cd();
    
    //2 is the number of parameters.
    TF1* Prob = new TF1(name2.c_str(),FunctionOutL,0,900,2);
    
    Prob->SetParameters(DeltaM41,700*pow(10,6));
    Prob->SetParNames("DeltaM41","Energy");
    Prob->GetXaxis()->SetTitle("Length of Neutrino Flight (m)");
    Prob->GetYaxis()->SetTitle("Nu mu to Nu e Oscillation Probability (%)");
    Prob->SetTitle(name3.c_str());
    
    Prob->Draw();
    
    //second graph.
    pad2->cd();
    
    TF1* Prob2 = new TF1(name4.c_str(),FunctionOutE,0,3,2);
    
    Prob2->SetParameters(DeltaM41,600);
    Prob2->SetParNames("DeltaM41","Distance");
    Prob2->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    Prob2->GetYaxis()->SetTitle("Nu mu to Nu e Oscillation Probability (%)");
    Prob2->SetLineColor(kBlue);
    
    
    Prob2->Draw();
    

}


//main function
void SBNProb(){
    
 
    GraphGenerator(0.4,"canvas1","graph1A","DeltaM41 = 0.4 eV^2","graph1B");

    GraphGenerator(1.1,"canvas2","graph2A","DeltaM41 = 1.1 eV^2","graph2B");

    GraphGenerator(6.0,"canvas3","graph3A","DeltaM41 = 6.0 eV^2","graph3B");

    GraphGenerator(20.0,"canvas4","graph4A","DeltaM41 = 20.0 eV^2","graph4B");
    
    
    
    

}
