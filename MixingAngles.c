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
#include "TF2.h"
#include <stdlib.h>

using namespace std;

//Function without any approximations
double FunctionOne(double *x, double *par){

    double thetamumu = x[0];
    double thetaee = x[1];
    
    double thetamue = (1-sqrt(1-thetamumu))*(1-sqrt(1-thetaee));
    
    return thetamue;


}

//Function with approximations. Takes into account that Umu4 << 1.
double FunctionTwo(double *x, double *par){

    double thetamumu = x[0];
    double thetaee = x[1];
    
    double thetamue = (thetamumu*thetaee)/4.;
    
    
    return thetamue;

}

void MixingAngles(){

    //Draws the function without approximations.
    TF2* angles = new TF2("angle",FunctionOne,0,0.5,0,0.5);
    
    angles->SetMinimum(0);
    angles->SetMaximum(0.09);
    
    angles->SetLineColor(kBlue);
    angles->GetXaxis()->SetTitle("sin^{2}(2#theta_{#mu#mu})");
    angles->GetYaxis()->SetTitle("sin^{2}(2#theta_{ee})");
    angles->GetZaxis()->SetTitle("sin^{2}(2#theta_{#mue})");
    angles->GetXaxis()->SetTitleOffset(2);
    angles->GetYaxis()->SetTitleOffset(2);
    angles->GetZaxis()->SetTitleOffset(1.4);
    angles->SetTitle("Green = With Approx., Blue = W/out Approx.");
    
    angles->Draw("SURF");
    
    //Draws the function with approximations in a different canvas.
    TF2* approx = new TF2("app",FunctionTwo,0,0.5,0,0.5);
    approx->SetMinimum(0);
    approx->SetMaximum(0.09);
    approx->SetLineColor(kGreen);
    
    approx->Draw("SURF,SAME");
    
    

}
