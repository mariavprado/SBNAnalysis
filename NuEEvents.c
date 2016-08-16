#include <stdio.h>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMinuit.h"
#include <errno.h>
#include <limits.h>
#include <TRandom3.h>
#include "TMath.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLegend.h"
#include <complex.h>
#include "TH2.h"
#include "TGraph.h"
#include <vector>


//void Rebinning(TH1D* hist, TH1D* hnew){
//
//    double xbin[12]={0.0,0.255,0.5,0.67,0.81,0.98,1.1,1.3,1.5,1.75,2.0,3.0};
//    int inc = 0;
//    
//    for (int i=0; i<12; i++) {
//        
//        
//        
//        double binwidth = hist->GetBinWidth(1);
//        
//        double NumInBin = TMath::Nint((xbin[i+1]-xbin[i])/binwidth);
//        
//        double EventPerBin = 0;
//        
//        for (int g=inc; g<(inc+NumInBin); g++) {
//            
//            EventPerBin = EventPerBin + hist->GetBinContent(g);
//            std::cout<<"      old center: "<<hist->GetBinCenter(g)<<std::endl;
//            
//        }
//        
//        //USE EDGES WITH IF STATEMENTS TO PLACE LIMITS 
//
//        //std::cout << "inc old: " << inc << std::endl;
//        std::cout << "xbin: " << xbin[i] << ", " << xbin[i+1] << std::endl;
//        //std::cout << "NumInBin: " << NumInBin << std::endl;
//
//        inc = inc + NumInBin;
//        hnew->SetBinContent(i,EventPerBin/(xbin[i+1]-xbin[i]));
//        //std::cout << "inc new: " << inc << std::endl;
//
//        //std::cout<<"new center: "<<hnew->GetBinCenter(i)<<std::endl;
//        double binHalfWidth = hnew->GetBinWidth(i)/2;
//        std::cout << "edges: " << hnew->GetBinCenter(i) - binHalfWidth << ", " << hnew->GetBinCenter(i) + binHalfWidth << std::endl;
//        
//        
//    }
//
//}

void Rescaling(TH1D* hist){

    for (int i=1; i<=hist->GetXaxis()->GetNbins(); i++) {
        double content = hist->GetBinContent(i);
        double binwidth = hist->GetBinWidth(i);
        
        hist->SetBinContent(i,content/binwidth);
        
        
    }

}

void ProbMuE(TH1D* hist, double DeltaM2, double L, double theta_mue){

    for (int i=1; i<=hist->GetXaxis()->GetNbins(); i++) {
        
        
        double EventOrig = hist->GetBinContent(i);
        
        
        double NuEnergy = hist->GetBinCenter(i);
        
        
        //theta_mue is really sin^2(2*theta_mue)
        double Prob = theta_mue*pow(sin(DeltaM2*L/(3*pow(10,8)*4*NuEnergy*pow(10,9)*6.582119514*pow(10,-16))),2);
        
        
        
        hist->SetBinContent(i,EventOrig*Prob);
        
        
        
    }
    

}

void ProbEE(TH1D* hist, double DeltaM2, double L, double theta_ee){
    
    
    for (int i=1; i<=hist->GetXaxis()->GetNbins(); i++) {
        
        double EventOrig = hist->GetBinContent(i);
        double NuEnergy = hist->GetBinCenter(i);

    
        double Prob = 1-theta_ee*pow(sin(DeltaM2*L/(3*pow(10,8)*4*NuEnergy*pow(10,9)*6.582119514*pow(10,-16))),2);
        
        
        
//        if (L==600) {
//            std::cout<<"bin center:   "<<NuEnergy<<std::endl;
//            std::cout<<"probEE:     "<<Prob<<std::endl;
//        }
        
        hist->SetBinContent(i,EventOrig*Prob);

        
    }


}


void ExpSignal(TH1D* hist1, TH1D* hist2, TH1D* hist3, TH1D* hist4, TH1D* hist5, TH1D* signalhist){

    for (int i=1; i<=hist1->GetXaxis()->GetNbins(); i++) {
        
        double event1 = hist1->GetBinContent(i);
        double event2 = hist2->GetBinContent(i);
        double event3 = hist3->GetBinContent(i);
        double event4 = hist4->GetBinContent(i);
        double event5 = hist5->GetBinContent(i);
        //double event6 = hist6->GetBinContent(i);
        //double event7 = hist7->GetBinContent(i);
        //double event8 = hist8->GetBinContent(i);
        
        double newcontent = event1+event2+event3+event4+event5;
        
        signalhist->SetBinContent(i,newcontent);
        
    }



}


double ChiSquare(TH1D* hist1, TH1D* hist2, TH1D* hist3, TH1D* hist4, TH1D* hist5, TH1D* hist6, TH1D* hist7, TH1D* hist8, TH1D* hist9, TH1D* sig1, TH1D* sig2, TH1D* sig3){
    
    std::cout<<"We have reached the chi square function"<<std::endl;

    std::vector <double> SBND;
    std::vector <double> MicroBooNE;
    std::vector <double> ICARUS;
    std::vector <double> Diff;
    std::vector <double> MCvar1;
    std::vector <double> MCvar2;
    std::vector <double> MCvar3;
    std::vector <double> MCAll;
    std::vector <std::vector <double> > Var;
    

    double NumHist = hist1->GetXaxis()->GetNbins();
    
    double handchi=0;

    


    
    for (int i=1; i<=NumHist; i++) {
        
        double MC1 = hist1->GetBinContent(i);
        double MC2 = hist2->GetBinContent(i);
        double MC3 = hist3->GetBinContent(i);
        double MC4 = hist4->GetBinContent(i);
        double MC5 = hist5->GetBinContent(i);
        double MC6 = hist6->GetBinContent(i);
        double MC7 = hist7->GetBinContent(i);
        double MC8 = hist8->GetBinContent(i);
        double MC9 = hist9->GetBinContent(i);
        
        double MCSBND = MC1+MC2+MC3;
        double MCMic = MC4+MC5+MC6;
        double MCICARUS = MC7+MC8+MC9;
        
//        std::cout<< "SBND green events for "<<i<<":   "<<MCSBND<<std::endl;
//        std::cout<< "MicroB real green events for "<<i<<":  "<<MCMic<<std::endl;
//        std::cout<< "ICARUS real green events for "<<i<<":  "<<MCICARUS<<std::endl;

        
        double DataSBND = sig1->GetBinContent(i);
        double DataMic = sig2->GetBinContent(i);
        double DataICARUS = sig3->GetBinContent(i);
        
//        std::cout<< "SBND data events for "<<i<<":   "<<DataSBND<<std::endl;
//        std::cout<< "MicroB data events for "<<i<<":  "<<DataMic<<std::endl;
//        std::cout<< "ICARUS data events for "<<i<<":  "<<DataICARUS<<std::endl;
//        
//        std::cout<< "SBND N-Nbar real events for "<<i<<":  "<<DataSBND-MCSBND<<std::endl;
//        std::cout<< "MicroB N-Nbar real events for "<<i<<":  "<<DataMic-MCMic<<std::endl;
//        std::cout<< "ICARUS N-Nbar real events for "<<i<<":  "<<DataICARUS-MCICARUS<<std::endl;
        
        
        
        SBND.push_back(DataSBND-MCSBND);
        MicroBooNE.push_back(DataMic-MCMic);
        ICARUS.push_back(DataICARUS-MCICARUS);
        
//        std::cout<<"SBND N-Nbar element "<<i<<":  "<<SBND[i-1]<<std::endl;
//        std::cout<<"MicroBooNE N-Nbar element "<<i<<":  "<<MicroBooNE[i-1]<<std::endl;
//        std::cout<<"ICARUS N-Nbar element "<<i<<":  "<<ICARUS[i-1]<<std::endl;


        
        MCvar1.push_back(MCSBND);
        MCvar2.push_back(MCMic);
        MCvar3.push_back(MCICARUS);
        
//        std::cout<<"SBND green events for "<<i<<":  "<<MCvar1[i-1]<<std::endl;
//        std::cout<<"MicroBooNE green events for "<<i<<":  "<<MCvar2[i-1]<<std::endl;
//        std::cout<<"ICARUS green events for "<<i<<":  "<<MCvar3[i-1]<<std::endl;
        
        
        
        handchi = handchi + (pow(DataSBND-MCSBND,2)/MCSBND)+(pow(DataMic-MCMic,2)/MCMic)+(pow(DataICARUS-MCICARUS,2)/MCICARUS);

        
    }
    
    
    Diff.insert(Diff.end(),SBND.begin(),SBND.end());
    Diff.insert(Diff.end(),MicroBooNE.begin(),MicroBooNE.end());
    Diff.insert(Diff.end(),ICARUS.begin(),ICARUS.end());
    
//    for (int i=0; i<NumHist*3; i++) {
//        std::cout<<"N-Nbar array:     "<<Diff[i]<<std::endl;
//
//    }
    
    
    MCAll.insert(MCAll.end(),MCvar1.begin(),MCvar1.end());
    MCAll.insert(MCAll.end(),MCvar2.begin(),MCvar2.end());
    MCAll.insert(MCAll.end(),MCvar3.begin(),MCvar3.end());

    
    Var.resize(NumHist*3);
    
    for (int i=0; i<NumHist*3; i++) {
        
        Var[i].resize(NumHist*3);
        
    }
    
    for (int i=0; i<NumHist*3; i++) {
        for (int g=0; g<NumHist*3; g++) {
            if(i==g){
                

                Var[i][g]=1/MCAll[i];
                
                
            } else {
                
                Var[i][g]=0;
                
            }
            
        }

    }
    
    
    
    
//    for (int i=0; i<NumHist*3; i++) {
//        for (int g=0; g<NumHist*3; g++) {
//            
//            std::cout<<Var[g][i]<<std::endl;
//
//        }
//        
//        std::cout<<"next column"<<std::endl;
//    
//    }
    
    
    
    
    
    
    double mid = 0;
    double chisq = 0;
    
    for (int i=0; i<NumHist*3; i++) {
        for (int j=0; j<NumHist*3; j++) {
            
            chisq += Diff[i]*Diff[j]*Var[j][i];
            
        }
        
    }
    
    
     double signif = sqrt(chisq);
    
    std::cout<<"hand chi square:    "<<handchi<<std::endl;
    std::cout<<"chi square:    "<<chisq<<std::endl;
    std::cout<<"significance:    "<<signif<<std::endl;

    
    
    return chisq;
    
}

void NuEEvents(){
    
    int NumBin = 12;
    double L1 = 110;
    double L2 = 470;
    double L3 = 600;
    //double m2 = 0.43;
    //double theta_mue = 0.013;
    double theta_ee = 0;
    

    TFile* file1 = new TFile("/Applications/ROOTfiles/ForDave/Nue_Plots/nominal_ntuples/hists/combined_ntuple_100m_nu_processed_nue.root_vePhot0.05_sg3_250w_0ecalo2_hists.root","READ");
    TFile* file2 = new TFile("/Applications/ROOTfiles/ForDave/Nue_Plots/nominal_ntuples/hists/combined_ntuple_470m_nu_processed_nue.root_vePhot0.05_sg3_250w_0ecalo2_hists.root","READ");
    TFile* file3 = new TFile("/Applications/ROOTfiles/ForDave/Nue_Plots/nominal_ntuples/hists/combined_ntuple_600m_onaxis_nu_processed_nue.root_vePhot0.05_sg3_250w_0ecalo2_hists.root","READ");
    TFile* file4 = new TFile("/Applications/ROOTfiles/ForDave/FinalNuMu_Plots/nominal_ntuples/LengthCut/combined_ntuple_100m_nu_processed_numu.root","READ");
    
    TFile* file5 = new TFile("/Applications/ROOTfiles/ForDave/FinalNuMu_Plots/nominal_ntuples/LengthCut/combined_ntuple_470m_nu_processed_numu.root","READ");
    
    TFile* file6 = new TFile("/Applications/ROOTfiles/ForDave/FinalNuMu_Plots/nominal_ntuples/LengthCut/combined_ntuple_600m_nu_processed_numu.root","READ");
    
    
    TH1D* hist1 = (TH1D*) file1->Get("kNueFromNueCC_muon1");
    TH1D* hist2 = (TH1D*) file1->Get("kNueFromNueCC_chargeKaon1");
    TH1D* hist3 = (TH1D*) file1->Get("kNueFromNueCC_neutKaon1");
    
    TH1D* hist1clone = (TH1D*) hist1->Clone();
    TH1D* hist2clone = (TH1D*) hist2->Clone();
    TH1D* hist3clone = (TH1D*) hist3->Clone();
    
    TH1D* hist1cloneCHI = (TH1D*) hist1->Clone();
    TH1D* hist2cloneCHI = (TH1D*) hist2->Clone();
    TH1D* hist3cloneCHI = (TH1D*) hist3->Clone();
    
    TH1D* hist4 = (TH1D*) file2->Get("kNueFromNueCC_muon2");
    TH1D* hist5 = (TH1D*) file2->Get("kNueFromNueCC_chargeKaon2");
    TH1D* hist6 = (TH1D*) file2->Get("kNueFromNueCC_neutKaon2");
    
    TH1D* hist4clone = (TH1D*) hist4->Clone();
    TH1D* hist5clone = (TH1D*) hist5->Clone();
    TH1D* hist6clone = (TH1D*) hist6->Clone();
    
    TH1D* hist4cloneCHI = (TH1D*) hist4->Clone();
    TH1D* hist5cloneCHI = (TH1D*) hist5->Clone();
    TH1D* hist6cloneCHI = (TH1D*) hist6->Clone();
    
    TH1D* hist7 = (TH1D*) file3->Get("kNueFromNueCC_muon3");
    TH1D* hist8 = (TH1D*) file3->Get("kNueFromNueCC_chargeKaon3");
    TH1D* hist9 = (TH1D*) file3->Get("kNueFromNueCC_neutKaon3");
    
    TH1D* hist7clone = (TH1D*) hist7->Clone();
    TH1D* hist8clone = (TH1D*) hist8->Clone();
    TH1D* hist9clone = (TH1D*) hist9->Clone();
    
    TH1D* hist7cloneCHI = (TH1D*) hist7->Clone();
    TH1D* hist8cloneCHI = (TH1D*) hist8->Clone();
    TH1D* hist9cloneCHI = (TH1D*) hist9->Clone();

    
    TH1D* signal1 = (TH1D*) file4->Get("NumuCC");
    TH1D* signal2 = (TH1D*) file5->Get("NumuCC");
    TH1D* signal3 = (TH1D*) file6->Get("NumuCC");
    
    
    //TH1D* back1 = (TH1D*) file1->Get("kDirt1");
    TH1D* back2 = (TH1D*) file1->Get("kNueFromNumuCC1");
    TH1D* back6 = (TH1D*) file2->Get("kNueFromNumuCC2");
    TH1D* back9 = (TH1D*) file3->Get("kNueFromNumuCC3");
    
    
    double xbin[13]={0.0,0.2,0.35,0.5,0.65,0.8,0.95,1.1,1.3,1.5,1.75,2.0,3.0};
    

     TCanvas* canvas = new TCanvas("canvas","Histogram1",200,20,1000,600);
    
    TH1D* clone1new = (TH1D*) hist1clone->Rebin(NumBin,"hist1old",xbin);
    TH1D* clone1CHI = (TH1D*) hist1cloneCHI->Rebin(NumBin,"CHI1",xbin);
    Rescaling(clone1new);
    clone1new->SetFillColor(32);
    
    TH1D* clone2new = (TH1D*) hist2clone->Rebin(NumBin,"hist2old",xbin);
    TH1D* clone2CHI = (TH1D*) hist2cloneCHI->Rebin(NumBin,"CHI2",xbin);
    Rescaling(clone2new);
    clone2new->SetFillColor(30);
    
    TH1D* clone3new = (TH1D*) hist3clone->Rebin(NumBin,"hist3old",xbin);
    TH1D* clone3CHI = (TH1D*) hist3cloneCHI->Rebin(NumBin,"CHI3",xbin);
    Rescaling(clone3new);
    clone3new->SetFillColor(kGreen);
    
    TH1D* clone4new = (TH1D*) hist4clone->Rebin(NumBin,"hist4old",xbin);
    TH1D* clone4CHI = (TH1D*) hist4cloneCHI->Rebin(NumBin,"CHI4",xbin);
    Rescaling(clone4new);
    clone4new->SetFillColor(32);
    
    TH1D* clone5new = (TH1D*) hist5clone->Rebin(NumBin,"hist5old",xbin);
    TH1D* clone5CHI = (TH1D*) hist5cloneCHI->Rebin(NumBin,"CHI5",xbin);
    Rescaling(clone5new);
    clone5new->SetFillColor(30);
    
    TH1D* clone6new = (TH1D*) hist6clone->Rebin(NumBin,"hist6old",xbin);
    TH1D* clone6CHI = (TH1D*) hist6cloneCHI->Rebin(NumBin,"CHI6",xbin);
    Rescaling(clone6new);
    clone6new->SetFillColor(kGreen);
    
    TH1D* clone7new = (TH1D*) hist7clone->Rebin(NumBin,"hist7old",xbin);
    TH1D* clone7CHI = (TH1D*) hist7cloneCHI->Rebin(NumBin,"CHI7",xbin);
    Rescaling(clone7new);
    clone7new->SetFillColor(32);
    
    TH1D* clone8new = (TH1D*) hist8clone->Rebin(NumBin,"hist8old",xbin);
    TH1D* clone8CHI = (TH1D*) hist8cloneCHI->Rebin(NumBin,"CHI8",xbin);
    Rescaling(clone8new);
    clone8new->SetFillColor(30);
    
    TH1D* clone9new = (TH1D*) hist9clone->Rebin(NumBin,"hist9old",xbin);
    TH1D* clone9CHI = (TH1D*) hist9cloneCHI->Rebin(NumBin,"CHI9",xbin);
    Rescaling(clone9new);
    clone9new->SetFillColor(kGreen);
    
    
    THStack *hs = new THStack("hs","SBND (110m)");
    
    hs->Add(clone1new);
    hs->Add(clone2new);
    hs->Add(clone3new);
    hs->Draw();
    
    hs->GetXaxis()->SetTitle("#nu_{e} Reconstructed Energy (GeV)");
    hs->GetXaxis()->SetTitleOffset(1);
    hs->GetYaxis()->SetTitle("Events per GeV");
    hs->GetYaxis()->SetTitleOffset(1.5);
    hs->SetMaximum(21000);
    
    TLegend* legendPostition = new TLegend(0.59, 0.58, 0.88, 0.89);
    legendPostition->AddEntry(clone1new, "#mu #rightarrow #nu_{e}");
    legendPostition->AddEntry(clone2new, "K^{+} #rightarrow #nu_{e}");
    legendPostition->AddEntry(clone3new, "K^{0} #rightarrow #nu_{e}");
    legendPostition->Draw();
    
    
    canvas->Modified();
    
    TCanvas* canvas1 = new TCanvas("canvas1","Histogram2",200,20,1000,600);
    
    THStack *hs1 = new THStack("hs1","MicroBooNE (470m)");
    
    hs1->Add(clone4new);
    hs1->Add(clone5new);
    hs1->Add(clone6new);
    hs1->Draw();
    
    hs1->GetXaxis()->SetTitle("#nu_{e} Reconstructed Energy (GeV)");
    hs1->GetXaxis()->SetTitleOffset(1);
    hs1->GetYaxis()->SetTitle("Events per GeV");
    hs1->GetYaxis()->SetTitleOffset(1.5);
    hs1->SetMaximum(1500);
    
    TLegend* legendPostition2 = new TLegend(0.59, 0.58, 0.88, 0.89);
    legendPostition2->AddEntry(clone4new, "#mu #rightarrow #nu_{e}");
    legendPostition2->AddEntry(clone5new, "K^{+} #rightarrow #nu_{e}");
    legendPostition2->AddEntry(clone6new, "K^{0} #rightarrow #nu_{e}");
    legendPostition2->Draw();
    
    canvas1->Modified();
    
    TCanvas* canvas2 = new TCanvas("canvas2","Histogram3",200,20,1000,600);
    
    THStack *hs2 = new THStack("hs2","ICARUS-T600 (600m)");
    
    hs2->Add(clone7new);
    hs2->Add(clone8new);
    hs2->Add(clone9new);
    hs2->Draw();
    
    hs2->GetXaxis()->SetTitle("#nu_{e} Reconstructed Energy (GeV)");
    hs2->GetXaxis()->SetTitleOffset(1);
    hs2->GetYaxis()->SetTitle("Events per GeV");
    hs2->GetYaxis()->SetTitleOffset(1.5);
    hs2->SetMaximum(3500);
    
    TLegend* legendPostition3 = new TLegend(0.59, 0.58, 0.88, 0.89);
    legendPostition3->AddEntry(clone7new, "#mu #rightarrow #nu_{e}");
    legendPostition3->AddEntry(clone8new, "K^{+} #rightarrow #nu_{e}");
    legendPostition3->AddEntry(clone9new, "K^{0} #rightarrow #nu_{e}");
    legendPostition3->Draw();
    
    canvas2->Modified();
    
    TCanvas* canvas3 = new TCanvas("canvas3","Histogram4",200,20,1000,600);

    
    TH1D* expsignal1 = new TH1D("signal","Data Signal SBND (110m)",NumBin,xbin);
    TH1D* expsignal2 = new TH1D("signal2","Data Signal MicroBooNE (470m)",NumBin,xbin);
    TH1D* expsignal3 = new TH1D("signal3","Data Signal ICARUS-T600 (600m)",NumBin,xbin);

    
    
    
    double M2min = pow(10,-2);
    double Thmin = pow(10,-4);
    double M2max = 100;
    double Thmax = 1;
    int numberStepsx = 50;
    int numberStepsy = 50;
    double thetaMUE[numberStepsx+1];
    double Chi[numberStepsx+1];



    
    TH2D* ChiSqSurface = new TH2D("chi","Chi Square Surface",numberStepsx,Thmin,Thmax,numberStepsy,M2min,M2max);
    
    TH1D* D1Chi = new TH1D("chi1","1 Dimensional Chi Square",numberStepsx,0,20000);
    
    
    //for (int i=0; i<=numberStepsy; i++) {
        for (int j=0; j<=numberStepsx; j++) {
            
            double m2 = 1;
        
            //log binning to span the entire chi square surface
           // double m2 = pow(10., (TMath::Log10(M2min)+ (i * (1./numberStepsy))*TMath::Log10(M2max/M2min)));
            double theta_mue = pow(10., (TMath::Log10(Thmin)+ (j * (1./numberStepsx))*TMath::Log10(Thmax/Thmin)));
            
            
            std::cout<<"m2:     "<<m2<<std::endl;
            std::cout<<"theta_mue:     "<<theta_mue<<std::endl;
    
            TH1D* hnew1 = (TH1D*) hist1->Clone();
            ProbEE(hnew1,m2,L1,theta_ee);
            hnew1->Rebin(NumBin,"hnew1",xbin);
            //Rescaling(hnew1);

            TH1D* hnew2 = (TH1D*) hist2->Clone();
            ProbEE(hnew2,m2,L1,theta_ee);
            hnew2->Rebin(NumBin,"hnew2",xbin);
            //Rescaling(hnew2);

            TH1D* hnew3 = (TH1D*) hist3->Clone();
            ProbEE(hnew3,m2,L1,theta_ee);
            hnew3->Rebin(NumBin,"hnew3",xbin);
            //Rescaling(hnew3);
    
            TH1D* newback2 = (TH1D*) back2->Clone();
            newback2->Rebin(NumBin,"back2",xbin);
            //Rescaling(newback2);

    
    
            
            TH1D* hnew4 = (TH1D*) hist4->Clone();
            ProbEE(hnew4,m2,L2,theta_ee);
            hnew4->Rebin(NumBin,"hnew4",xbin);
            //Rescaling(hnew4);

            TH1D* hnew5 = (TH1D*) hist5->Clone();
            ProbEE(hnew5,m2,L2,theta_ee);
            hnew5->Rebin(NumBin,"hnew5",xbin);
            //Rescaling(hnew5);
    
            TH1D* hnew6 = (TH1D*) hist6->Clone();
            ProbEE(hnew6,m2,L2,theta_ee);
            hnew6->Rebin(NumBin,"hnew6",xbin);
            //Rescaling(hnew6);
    
            TH1D* newback6 = (TH1D*) back6->Clone();
            newback6->Rebin(NumBin,"back6",xbin);
            //Rescaling(newback6);
    

    
            TH1D* hnew7 = (TH1D*) hist7->Clone();
            ProbEE(hnew7,m2,L3,theta_ee);
            hnew7->Rebin(NumBin,"hnew7",xbin);
            //Rescaling(hnew7);
    
            TH1D* hnew8 = (TH1D*) hist8->Clone();
            ProbEE(hnew8,m2,L3,theta_ee);
            hnew8->Rebin(NumBin,"hnew8",xbin);
            //Rescaling(hnew8);

            TH1D* hnew9 = (TH1D*) hist9->Clone();
            ProbEE(hnew9,m2,L3,theta_ee);
            hnew9->Rebin(NumBin,"hnew9",xbin);
            //Rescaling(hnew9);
    
            TH1D* newback9 = (TH1D*) back9->Clone();
            newback9->Rebin(NumBin,"back9",xbin);
            //Rescaling(newback9);

    

    
            TH1D* newsig1 = (TH1D*) signal1->Clone();
            ProbMuE(newsig1,m2,L1,theta_mue);
            newsig1->Rebin(NumBin,"signal1",xbin);
            //Rescaling(newsig1);

            TH1D* newsig2 = (TH1D*) signal2->Clone();
            ProbMuE(newsig2,m2,L2,theta_mue);
            newsig2->Rebin(NumBin,"signal2",xbin);
            //Rescaling(newsig2);
    
            TH1D* newsig3 = (TH1D*) signal3->Clone();
            ProbMuE(newsig3,m2,L3,theta_mue);
            newsig3->Rebin(NumBin,"signal3",xbin);
            //Rescaling(newsig3);

    
            ExpSignal(hnew1,hnew2,hnew3,newback2,newsig1,expsignal1);
            ExpSignal(hnew4,hnew5,hnew6,newback6,newsig2,expsignal2);
            ExpSignal(hnew7,hnew8,hnew9,newback9,newsig3,expsignal3);
    
    
    
            
            double chisq = ChiSquare(clone1CHI, clone2CHI, clone3CHI, clone4CHI, clone5CHI, clone6CHI, clone7CHI, clone8CHI, clone9CHI, expsignal1, expsignal2, expsignal3);
            

            std::cout<<"theta_mue after the chi sq function"<<theta_mue<<std::endl;
            std::cout<<"m2 after the chi sq function"<<m2<<std::endl;
            std::cout<<"chisq after the chi sq function"<<chisq<<std::endl;
            
            
            D1Chi->Fill(chisq);

            //ChiSqSurface->SetBinContent(j+1,i+1,chisq);
            
            
        
           thetaMUE[j] = theta_mue;
           Chi[j] = chisq;
            
            


    
    
    
    
    //CODE FOR IF YOU WANT TO PLOT THE SIGNAL PLOTS AS WELL. THEY ARE NOT GOING TO BE ORIGINALLY PLOTTED BECAUSE WITH THE LOOP THAT IS NEEDED TO MAKE THE CHI SQUARE SURFACE THERE WOULD BE AN INFINITE AMOUNT OF SIGNAL PLOTS FOR EVERY COMBINATION OF M^2 AND SIN^2(2*THETA_MUE)
    
//    TCanvas* canvasSig = new TCanvas("canvassig","theta_ee=0.9",200,20,1000,600);
//    expsignalCL1->GetXaxis()->SetTitle("#nu_{e} Reconstructed Energy (GeV)");
//    expsignalCL1->GetXaxis()->SetTitleOffset(1);
//    expsignalCL1->GetYaxis()->SetTitle("Events per GeV");
//    expsignalCL1->GetYaxis()->SetTitleOffset(1.5);
//    expsignalCL1->SetMaximum(22000);
//    expsignalCL1->Draw();
//    
//    TCanvas* canvasSig2 = new TCanvas("canvassig2","theta_ee=0.9",200,20,1000,600);
//    expsignalCL2->GetXaxis()->SetTitle("#nu_{e} Reconstructed Energy (GeV)");
//    expsignalCL2->GetXaxis()->SetTitleOffset(1);
//    expsignalCL2->GetYaxis()->SetTitle("Events per GeV");
//    expsignalCL2->GetYaxis()->SetTitleOffset(1.5);
//    expsignalCL2->SetMaximum(1500);
//    expsignalCL2->Draw();
//    
//    TCanvas* canvasSig3 = new TCanvas("canvassig3","theta_ee=0.9",200,20,1000,600);
//    expsignalCL3->GetXaxis()->SetTitle("#nu_{e} Reconstructed Energy (GeV)");
//    expsignalCL3->GetXaxis()->SetTitleOffset(1);
//    expsignalCL3->GetYaxis()->SetTitle("Events per GeV");
//    expsignalCL3->GetYaxis()->SetTitleOffset(1.5);
//    expsignalCL3->SetMaximum(3500);
//    expsignalCL3->Draw();
    

    
     //   }
    //}
    
    
    canvas3->SetLogx();
    canvas3->SetLogy();
    canvas3->SetLogz();
    ChiSqSurface->Draw("COLZ");
    
    //ChiSqSurface->Draw("SURF3");

    
    canvas3->Modified();

    
    TCanvas* canvas4 = new TCanvas("canvas4","Histogram5",200,20,1000,600);
    
    D1Chi->Draw();
    
    TCanvas* canvas5 = new TCanvas("canvas5","Scatter",200,20,1000,600);
    TGraph* g = new TGraph(numberStepsx+1,thetaMUE,Chi);
    g->SetLineColor(2);
    g->SetLineWidth(4);
    g->SetMarkerColor(4);
    g->SetMarkerSize(1.5);
    g->SetMarkerStyle(21);
    g->GetXaxis()->SetTitle("sin^{2}(2 #theta_{#mu e})");
    g->GetYaxis()->SetTitle("#chi^{2}");
    
    g->Draw("ap");


    


}
