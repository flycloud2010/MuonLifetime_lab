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



void MC_CDF(int N_sim, bool save = true, bool rebin = false)
{   

    const int N_events = 140000; // number of events in one week

    double Begin = 179.5;
    double End = 2048.5;

    double Begin_fit;
    double End_fit = 2048.5;

    int nbins;

    (rebin) ? (nbins = 310) : (nbins = End - Begin );


    cout << "Nbins: " << nbins << endl;
    cout << "Nentries " << N_events << endl;
    cout << "N_sim " << N_sim << endl;

    TString filename = "taulong_invcdf.root";

    TFile* outfile = new TFile("root_files/" + filename,"RECREATE");
    TTree* invcdf = new TTree("expo_taulong","Tree"); 

    int k = 5;

    double tau_chi[5], tau_chi_err[5];      // chi square method
    double tau_chip[5], tau_chip_err[5];    // Pearson chi square
    double tau_LL[5], tau_LL_err[5];        // Log_Likelihood
    double delta_chi[5], delta_chip[5], delta_L[5];
    double entries;

    invcdf->Branch("tau_chi",tau_chi,"tau_chi[5]/D");
    invcdf->Branch("tau_chi_err",tau_chi_err,"tau_chi_err[5]/D");
    invcdf->Branch("tau_chip",tau_chip,"tau_chip[5]/D");
    invcdf->Branch("tau_chip_err",tau_chip_err,"tau_chip_err[5]/D");
    invcdf->Branch("tau_LL",tau_LL,"tau_LL[5]/D");
    invcdf->Branch("tau_LL_err",tau_LL_err,"tau_LL_err[5]/D");
    invcdf->Branch("delta_chi",delta_chi,"delta_chi[5]/D");
    invcdf->Branch("delta_chip",delta_chip,"delta_chip[5]/D");
    invcdf->Branch("delta_L", delta_L,"delta_L[5]/D");

    
    TRandom3* r = new TRandom3;

    progressbar bar(N_sim); 

    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");

    
    TF1* f = nullptr;
    TH1D* h = nullptr;

    double y, t;
    int n_tau;

    for (int j = 0; j < N_sim; j++)
        {
            bar.update();

            h = new TH1D("expMC","expMC",nbins,0,End);

            r->SetSeed(j+1);

            for (int k = 0; k < N_events; k++ ) 
            {   
                y = r->Uniform(0,1);
                t = - tau_long * log(1-y) + 180.;
                h->Fill(t);
            }

            //h->Draw();

            for (k = 0; k < 5; k++)
            
            {   
                n_tau = 5 + k;

                Begin_fit = n_tau * tau_short;
                f = new TF1("f","expo",Begin_fit,End);

                // ---- CHI SQUARE METHOD ---- //
            
                TFitResultPtr tauPtr1 = h->Fit(f, "RNSQ");;
                tau_chi[k] = -1/tauPtr1->Parameter(1);
                tau_chi_err[k] = tauPtr1->ParError(1)/pow(tauPtr1->Parameter(1),2);
                delta_chi[k] = abs(tau_chi[k] - tau_long) / tau_chi_err[k];


                 // ---- PEARSON CHI SQUARE METHOD ---- //

                TFitResultPtr tauPtr2 = h->Fit(f, "RPSNQ");;
                tau_chip[k] = -1/tauPtr2->Parameter(1);
                tau_chip_err[k] = tauPtr2->ParError(1)/pow(tauPtr2->Parameter(1),2);;
                delta_chip[k] = abs(tau_chip[k] - tau_long) / tau_chip_err[k];


                // ---- LOG LIKELIHOOD METHOD ---- //

                TFitResultPtr tauPtr3 = h->Fit(f, "RLSNQ");;
                tau_LL[k] = -1/tauPtr3->Parameter(1);
                tau_LL_err[k] = tauPtr3->ParError(1)/pow(tauPtr3->Parameter(1),2);;
                delta_L[k] = abs(tau_LL[k] - tau_long) / tau_LL_err[k];


            }

            invcdf->Fill();

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
        invcdf->Write();  
        outfile->Close(); 
    }
}