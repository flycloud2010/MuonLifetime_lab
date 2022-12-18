#include <iostream>
#include <vector>
#include <string>

#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

#include "progressbar.hpp"

double t0 = 180;

double ToADC(double time, double offset, double slope)
{
    double adc = (time - offset) / slope;
    return(adc);
}

Double_t MyExp(Double_t* x, Double_t* par)
{
    Double_t fit = par[0]*TMath::Exp((x[0]-t0)*par[1]) + par[2];
    return fit;

}

Double_t DoubleExp(Double_t* x, Double_t* par)
{
    Double_t fit = par[0]*TMath::Exp((x[0]-t0)*par[1]) + par[2]*TMath::Exp((x[0]-t0)*par[3]) + par[4];
    return fit;

}


// parameters from TAC calibration
double a = 0.01588903364201166;
double b = 0.011793768586716;
double tau_short = ToADC(0.8466,a,b);
double tau_long = ToADC(2.1969811,a,b);
double R = 1.26;



void Fit(int N_sim=1, bool save = false, bool rebin = false, int weeks = 1, int n_tau_short = 4)
{   

    const int N_muons = weeks * 130000;    // number of events (mu+ and mu-) in one week
    
    //const int Np = weeks * 100000;       // mu +
    //const int Nm = N_muons - Np;         // mu-
    
    //double Ap = Np/tau_long;
    //double Am = Nm/tau_short;
    
    double Am = N_muons / (R*tau_long + tau_short);
    double Ap = R*Am;
    
    int Nm = Am*tau_short;
    int Np = Ap*tau_long;


    cout << "Nm " << Nm << endl;
    cout << "Np " << Np << endl;

    const double baseline = weeks * 4.;

    double Begin = 180.5;
    double End = 2048.5;

    const int N_base = baseline * (End - Begin);
    const int N_events = N_muons + N_base;

    double Begin_fit;
    double End_fit = 2048.5;

    double start_fit_expo = t0 + n_tau_short*tau_short;

    int nbins;
    float rebin_fact;


    if(rebin)
    {   
        nbins = 600;
        rebin_fact = (End - Begin) / 600;

    }
    else 
    {
        nbins = End - Begin;
        rebin_fact = 1;
    }

    cout << "Nbins: " << nbins << endl;
    cout << "Nentries " << N_events << endl;
    cout << "N_sim " << N_sim << endl;
    cout << "weeks: "  << weeks << endl;

    TString filename = Form("double_exp_%d_new_norebin.root", weeks);

    TFile* outfile = new TFile("root_files/" + filename,"UPDATE");
    TTree* outtree = new TTree(Form("exp_%d",n_tau_short),"Tree"); 


    double bs[10], bs_err[10];       
    double tau_plus[10], tau_plus_err[10];        
    double A_plus[10], A_plus_err[10];
    double A_minus[10], A_minus_err[10];
    double tau_minus[10], tau_minus_err[10];
    double delta_tau_minus[10], delta_tau_plus[10];
    double delta_bs[10];
    double delta_A_minus[10], delta_A_plus[10];
    double relerr_bs[10], relerr_tau_plus[10], relerr_tau_minus[10], relerr_A_plus[10], relerr_A_minus[10];


    outtree->Branch("bs",bs,"bs[10]/D");
    outtree->Branch("bs_err",bs_err,"bs_err[10]/D");
    outtree->Branch("delta_bs", delta_bs,"delta_bs[10]/D");
    outtree->Branch("relerr_bs", relerr_bs,"relerr_bs[10]/D");

    outtree->Branch("tau_plus",tau_plus,"tau_plus[10]/D");
    outtree->Branch("tau_plus_err",tau_plus_err,"tau_plus_err[10]/D");
    outtree->Branch("delta_tau_plus", delta_tau_plus,"delta_tau_plus[10]/D");
    outtree->Branch("relerr_tau_plus", relerr_tau_plus,"relerr_tau_plus[10]/D");

    outtree->Branch("tau_minus",tau_minus,"tau_minus[10]/D");
    outtree->Branch("tau_minus_err",tau_minus_err,"tau_minus_err[10]/D");
    outtree->Branch("delta_tau_minus", delta_tau_minus,"delta_tau_minus[10]/D");
    outtree->Branch("relerr_tau_minus", relerr_tau_minus,"relerr_tau_minus[10]/D");

    outtree->Branch("A_plus",A_plus,"A_plus[10]/D");
    outtree->Branch("A_plus_err",A_plus_err,"A_plus_err[10]/D");
    outtree->Branch("delta_A_plus", delta_A_plus,"delta_A_plus[10]/D");
    outtree->Branch("relerr_A_plus", relerr_A_plus,"relerr_A_plus[10]/D");

    outtree->Branch("A_minus",A_minus,"A_minus[10]/D");
    outtree->Branch("A_minus_err", A_minus_err,"A_minus_err[10]/D");
    outtree->Branch("delta_A_minus", delta_A_minus,"delta_A_minus[10]/D");
    outtree->Branch("relerr_A_minus", relerr_A_minus,"relerr_A_minus[10]/D");

    TRandom3* r = new TRandom3;

    progressbar bar(N_sim); 

    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");

    TH1D* h = nullptr;

    double y, t;
    float n_tau;
    double A_fixed;

    for (int j = 0; j < N_sim; j++)
        {
            bar.update();

            h = new TH1D("exp_base_histo","exp_base_histo",nbins,t0,End);

            r->SetSeed(j+1);

            
            for (int k = 0; k < N_base; k++ ) 
            {   
                h->Fill(r->Uniform(Begin,End));
            }

            for (int k = 0; k < Np; k++ ) 
            {   
                y = r->Uniform(0,1);
                t = - tau_long * log(1-y) + t0;
                h->Fill(t);
            }

            for (int k = 0; k < Nm; k++ ) 
            {   
                y = r->Uniform(0,1);
                t = - tau_short * log(1-y) + t0;
                h->Fill(t);
            }

            //h->Draw();
        

            for (int k = 0; k < 10; k++)
            
            {   

                TFitResultPtr basPtr = nullptr;
                TFitResultPtr tauPtr_long = nullptr;
                TFitResultPtr tauPtr_short = nullptr;
                TFitResultPtr totalPtr = nullptr;

                TF1* f = nullptr;
                TF1* f2 = nullptr;
                TF1* f1 = nullptr;
                TF1* f0 = nullptr;

                n_tau = 5 + 0.5*k;
    
                Begin_fit = t0 + n_tau * tau_long;

                f0 = new TF1("f0","pol0",Begin_fit,End_fit);

                // ---- 1) baseline ---- //

                basPtr = h->Fit(f0, "RLSNQ");
                
              
                // ---- 2) Long Exp with setpar for baseline ---- //
            
                f = new TF1("f",MyExp,start_fit_expo,End,3);

                //cout << "bin content " << h->GetBinContent(Begin_fit) << endl;
            
                //f->SetParLimits(0,Ap*rebin_fact*(1-0.01), Ap*rebin_fact*(1+0.01));
                //f->SetParLimits(1,-1/(tau_long*(1-0.0001)), 0);
                //cout << Begin_fit << endl;
                //f->SetParameter(0,Ap);
                f->SetParameter(1,-1/(tau_long*2));
                f->SetParameter(2,basPtr->Parameter(0));
                f->SetParLimits(2,basPtr->Parameter(0)-4*basPtr->ParError(0),basPtr->Parameter(0)+4*basPtr->ParError(0));

                tauPtr_long = h->Fit(f, "RLSNQ");

                f->SetLineColor(kRed);
                //f->Draw("same");


                // ---- 3) Short Exp prefit---- //
                f1 = new TF1("f1",DoubleExp,t0,End,5);
                
                
                //cout << "bin content " << h->GetBinContent(1) << endl;
                //cout << "max content " << h->GetMaximumBin() << endl;

                A_fixed = h->GetBinContent(1) - tauPtr_long->Parameter(0) - tauPtr_long->Parameter(2);

                //cout << "a minus: " << A_fixed << endl;
                //cout << "a minus th: " << Am << endl;
                


                f1->SetParameter(0,tauPtr_long->Parameter(0));
                f1->SetParameter(1,tauPtr_long->Parameter(1));
                f1->SetParameter(2,A_fixed);
                f1->SetParameter(3,-1/(tau_short*2));
                f1->SetParameter(4,tauPtr_long->Parameter(2));

                f1->SetParLimits(0, tauPtr_long->Parameter(0)-4*tauPtr_long->ParError(0),tauPtr_long->Parameter(0)+4*tauPtr_long->ParError(0));
                //f1->SetParLimits(3,-1/(tau_short*(1-0.2)), -1/(tau_short*(1+0.2)));
                //f1->SetParLimits(2,Am*rebin_fact*(1-0.2), Am*rebin_fact*(1+0.2));
                f1->SetParLimits(1, tauPtr_long->Parameter(1)-4*tauPtr_long->ParError(1),tauPtr_long->Parameter(1)+4*tauPtr_long->ParError(1));
                f1->SetParLimits(4, tauPtr_long->Parameter(2)-4*tauPtr_long->ParError(2), tauPtr_long->Parameter(2)+4*tauPtr_long->ParError(2));
                
                tauPtr_short = h->Fit(f1, "RLSNQ");

                //cout <<"taushort fit 2 " << -1/tauPtr_short->Parameter(3) << endl;

                f1->SetLineColor(kBlue);
                //f1->Draw("same");

                // ---- 4) Full fit---- //
                f2 = new TF1("f2",DoubleExp,t0,End,5);
                f2->SetParameters(tauPtr_short->Parameter(0), tauPtr_short->Parameter(1), tauPtr_short->Parameter(2), tauPtr_short->Parameter(3), tauPtr_short->Parameter(4));

                f2->SetParLimits(0,tauPtr_short->Parameter(0)-4*tauPtr_short->ParError(0), tauPtr_short->Parameter(0)+4*tauPtr_short->ParError(0));
                f2->SetParLimits(1,tauPtr_short->Parameter(1)-4*tauPtr_short->ParError(1), tauPtr_short->Parameter(1)+4*tauPtr_short->ParError(1));
                f2->SetParLimits(2,tauPtr_short->Parameter(2)-4*tauPtr_short->ParError(2), tauPtr_short->Parameter(2)+4*tauPtr_short->ParError(2));
                f2->SetParLimits(3,tauPtr_short->Parameter(3)-4*tauPtr_short->ParError(3), tauPtr_short->Parameter(3)+4*tauPtr_short->ParError(3));
                f2->SetParLimits(4,tauPtr_short->Parameter(4)-4*tauPtr_short->ParError(4), tauPtr_short->Parameter(4)+4*tauPtr_short->ParError(4));

                totalPtr = h->Fit(f2, "RLSNQ");
                //cout <<"taushort fit 3 " << -1/totalPtr->Parameter(3) << endl;

                A_plus[k] = totalPtr->Parameter(0);
                A_plus_err[k] = totalPtr->ParError(0);
                delta_A_plus[k] = abs(A_plus[k] - Ap) / A_plus_err[k];
                relerr_A_plus[k] = A_plus_err[k]/A_plus[k];

                tau_plus[k] = -1/totalPtr->Parameter(1);
                tau_plus_err[k] = totalPtr->ParError(1)/pow(totalPtr->Parameter(1),2);;
                delta_tau_plus[k] = abs(tau_plus[k] - tau_long) / tau_plus_err[k];
                relerr_tau_plus[k] = tau_plus_err[k]/tau_plus[k];

                A_minus[k] = totalPtr->Parameter(2);
                A_minus_err[k] = totalPtr->ParError(2);
                delta_A_minus[k] = abs(A_minus[k] - Am) / A_minus_err[k];
                relerr_A_minus[k] = A_minus_err[k]/A_minus[k];

                tau_minus[k] = -1/totalPtr->Parameter(3);
                tau_minus_err[k] = totalPtr->ParError(3)/pow(totalPtr->Parameter(3),2);;
                delta_tau_minus[k] = abs(tau_minus[k] - tau_short) / tau_minus_err[k];
                relerr_tau_minus[k] = tau_minus_err[k]/tau_minus[k];

                bs[k] = totalPtr->Parameter(4);
                bs_err[k] = totalPtr->ParError(4);
                delta_bs[k] = abs(bs[k] - baseline) / bs_err[k];
                relerr_bs[k] = bs_err[k]/bs[k];


                //cout <<"taushort delta " << delta_tau_minus[k] << endl;
                //cout <<"taulong delta " << delta_tau_plus[k] << endl;
                //cout <<"a- delta " << delta_A_minus[k] << endl;
                //cout <<"a+ delta " << delta_A_plus[k] << endl;
                //cout <<"bs delta " << delta_bs[k] << endl;
                //cout <<"a- rel " << relerr_A_minus[k] << endl;

                f2->SetLineColor(kGreen);
                //f2->Draw("same");

            }
            
            
            outtree->Fill();
            
            if( j == N_sim-1) 
            {   
                h->Write();
            }
//
            else h->Delete();
            
        }

    cout << endl;
    

    if(save) 
    {   
        cout << "Writing " << filename << endl;
        outtree->Write();  
        outfile->Close(); 
    }

}