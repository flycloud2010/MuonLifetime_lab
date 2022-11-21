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



double ToADC(double time, double offset, double slope)
{
    double adc = (time - offset) / slope;
    return(adc);
}

double MyExp(double* x, double* par)
{
    double fit;
    fit = par[0]*TMath::Exp(-x[0]/par[1]);
    return(fit);

}


// parameters from TAC calibration
double a = 0.02588903364201166;
double b = 0.011793768586716;
double tau_short = ToADC(0.8466,a,b);
double tau_long = ToADC(2.1969811,a,b);



void ExpBase(int N_sim, bool save = true, bool rebin = false, int weeks = 1, int n_tau_short = 4)
{   

    
    const int N_events = 61000; // number of events (mu+) in one week
    const double baseline = weeks * 4.;

    double Begin = 179.5;
    double End = 2048.5;

    const int N_base = baseline * (End - Begin);

    double Begin_fit;
    double End_fit = 2048.5;

    double start_fit_expo = n_tau_short*tau_short;

    int nbins;
    float rebin_fact;

    if(rebin)
    {   
        nbins = 310;
        rebin_fact = (End - Begin) / 310;

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

    TString filename = Form("exp_base_v2_%d.root", weeks);

    TFile* outfile = new TFile("root_files/" + filename,"UPDATE");
    TTree* outtree = new TTree(Form("baseline_exp_%d",n_tau_short),"Tree"); 

    double bs_LL[10], bs_LL_err[10];        // Log_Likelihood
    double bs_LL_relerr[10], bs_delta_L[10];
    double tau_LL[10], tau_LL_err[10];        
    double tau_LL_relerr[10], tau_delta_L[10];
    double tau_LL_noB[10], tau_LL_err_noB[10];        
    double tau_LL_relerr_noB[10], tau_delta_L_noB[10];

    outtree->Branch("bs_LL",bs_LL,"bs_LL[10]/D");
    outtree->Branch("bs_LL_err",bs_LL_err,"bs_LL_err[10]/D");
    outtree->Branch("bs_LL_relerr",bs_LL_relerr,"bs_LL_relerr[10]/D");
    outtree->Branch("bs_delta_L", bs_delta_L,"bs_delta_L[10]/D");

    outtree->Branch("tau_LL",tau_LL,"tau_LL[10]/D");
    outtree->Branch("tau_LL_err",tau_LL_err,"tau_LL_err[10]/D");
    outtree->Branch("tau_LL_relerr",tau_LL_relerr,"tau_LL_relerr[10]/D");
    outtree->Branch("tau_delta_L", tau_delta_L,"tau_delta_L[10]/D");

    outtree->Branch("tau_LL_noB",tau_LL_noB,"tau_LL_noB[10]/D");
    outtree->Branch("tau_LL_err_noB",tau_LL_err_noB,"tau_LL_err_noB[10]/D");
    outtree->Branch("tau_LL_relerr_noB",tau_LL_relerr_noB,"tau_LL_relerr_noB[10]/D");
    outtree->Branch("tau_delta_L_noB", tau_delta_L_noB,"tau_delta_L_noB[10]/D");

    TRandom3* r = new TRandom3;

    progressbar bar(N_sim); 

    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");

    TH1D* h = nullptr;

    double y, t;
    float n_tau;

    for (int j = 0; j < N_sim; j++)
        {
            bar.update();

            h = new TH1D("exp_base_histo","exp_base_histo",nbins,0,End);

            r->SetSeed(j+1);

            
            for (int k = 0; k < N_base; k++ ) 
            {   
                h->Fill(r->Uniform(Begin,End));
            }

            for (int k = 0; k < (N_events-N_base); k++ ) 
            {   
                y = r->Uniform(0,1);
                t = - tau_long * log(1-y) + 180.;
                h->Fill(t);
            }

            //h->Draw();

            
            for (int k = 0; k < 10; k++)
            
            {   

                TFitResultPtr basPtr = nullptr;
                TFitResultPtr tauPtr = nullptr;
                TFitResultPtr tauPtr1 = nullptr;

                TF1* f = nullptr;
                TF1* f1 = nullptr;
                TF1* f0 = nullptr;

                n_tau = 4 + 0.5*k;
                Begin_fit = 180. + n_tau * tau_long;
                f0 = new TF1("f0","pol0",Begin_fit,End_fit);

                // ---- baseline ---- //

                basPtr = h->Fit(f0, "RLSNQ");
                
                bs_LL[k] = basPtr->Parameter(0);
                bs_LL_err[k] = basPtr->ParError(0);
                bs_LL_relerr[k] = basPtr->ParError(0)/basPtr->Parameter(0);
                bs_delta_L[k] = abs(bs_LL[k] - baseline*rebin_fact) / bs_LL_err[k];

                //f0->Draw("same");
                // ---- exp with setpar for baseline ---- //
                f = new TF1("f","expo(0)+pol0(2)",start_fit_expo,End);
                f->SetParameter(bs_LL[k],2);

                tauPtr = h->Fit(f, "RLSNQ");
                tau_LL[k] = -1/tauPtr->Parameter(1);
                tau_LL_err[k] = tauPtr->ParError(1)/pow(tauPtr->Parameter(1),2);;
                tau_LL_relerr[k] = tauPtr->ParError(0)/tauPtr->Parameter(0);
                tau_delta_L[k] = abs(tau_LL[k] - tau_long) / tau_LL_err[k];

                // ---- complete fit exp+baseline ---- //
                f1 = new TF1("f","expo(0)+pol0(2)",start_fit_expo,End);
                tauPtr1 = h->Fit(f1, "RLSNQ");
                tau_LL_noB[k] = -1/tauPtr1->Parameter(1);
                tau_LL_err_noB[k] = tauPtr1->ParError(1)/pow(tauPtr1->Parameter(1),2);;
                tau_LL_relerr_noB[k] = tauPtr1->ParError(0)/tauPtr1->Parameter(0);
                tau_delta_L_noB[k] = abs(tau_LL_noB[k] - tau_long) / tau_LL_err_noB[k];

                //f->Draw("same");


            }
            
            
            outtree->Fill();
            
            if( j == N_sim-1) 
            {   
                h->Write();
            }

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