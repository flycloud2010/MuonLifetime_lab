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

void exponentialMC(int N_sim)
{ 
    const int N_events = 140000; // number of events in one week

    double Begin = 179.5;
    double End = 2048.5;

    double Begin_fit, End_fit;


    int nbins;

    nbins = (End - Begin);
   

    cout << "bins: " << nbins << endl;
    cout << "Nentries " << N_events << endl;

    TFile* outfile = new TFile("root_files/taushort_trandom.root","RECREATE");
    TTree* trandom = new TTree("expo_taushort","Tree"); 


    double tau_chi, tau_chi_err;      // chi square method
    double tau_chip, tau_chip_err;    // Pearson chi square
    double tau_LL, tau_LL_err;        // Log_Likelihood
    double delta_chi, delta_chip, delta_L;


    trandom->Branch("tau_chi",&tau_chi,"tau_chi/D");
    trandom->Branch("tau_chi_err",&tau_chi_err,"tau_chi_err/D");
    trandom->Branch("tau_chip",&tau_chip,"tau_chip/D");
    trandom->Branch("tau_chip_err",&tau_chip_err,"tau_chip_err/D");
    trandom->Branch("tau_LL",&tau_LL,"tau_LL/D");
    trandom->Branch("tau_LL_err",&tau_LL_err,"tau_LL_err/D");
    trandom->Branch("delta_chi",&delta_chi,"delta_chi/D");
    trandom->Branch("delta_chip",&delta_chip,"delta_chip/D");
    trandom->Branch("delta_L",&delta_L,"delta_L/D");


    TRandom3* r = new TRandom3;

    progressbar bar(N_sim); 

    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");

    //TCanvas* c = new TCanvas("c", "Simulations",1200,600);

    TF1* f = new TF1("f","expo",Begin,End);

    for (int j = 0; j < N_sim; j++)
        {
            bar.update();

            TH1D* h = new TH1D("exp","exp",nbins,Begin,End);

            r->SetSeed(j+1);

            for (int k = 0; k < N_events; k++ ) 
            {   
                h->Fill(r->Exp(tau_short));
            }

            //h->Draw();

            // ---- CHI SQUARE METHOD ---- //
            
            TFitResultPtr tauPtr1 = h->Fit(f, "RNSQ");;
            tau_chi = -1/tauPtr1->Parameter(1);
            tau_chi_err = tauPtr1->ParError(1)/pow(tauPtr1->Parameter(1),2);
            delta_chi = abs(tau_chi - tau_short) / tau_chi_err;
            

             // ---- PEARSON CHI SQUARE METHOD ---- //

            TFitResultPtr tauPtr2 = h->Fit(f, "RPSNQ");;
            tau_chip = -1/tauPtr2->Parameter(1);
            tau_chip_err = tauPtr2->ParError(1)/pow(tauPtr2->Parameter(1),2);;
            delta_chip = abs(tau_chip - tau_short) / tau_chip_err;


            // ---- LOG LIKELIHOOD METHOD ---- //

            TFitResultPtr tauPtr3 = h->Fit(f, "RLSNQ");;
            tau_LL = -1/tauPtr3->Parameter(1);
            tau_LL_err = tauPtr3->ParError(1)/pow(tauPtr3->Parameter(1),2);;
            delta_L = abs(tau_LL - tau_short) / tau_LL_err;

            trandom->Fill();

            //h->Draw();
            h->Delete();

        }

    cout << endl;

    trandom->Write();
    
    outfile->Close();



}



void MC_CDF(int N_sim, int nbins)
{   

    const int N_events = 140000; // number of events in one week

    double Begin = 179.5;
    double End = 2048.5;

    double Begin_fit, End_fit;


    //int nbins;

    //nbins = (End - Begin);
   

    cout << "bins: " << nbins << endl;
    cout << "Nentries " << N_events << endl;

    TFile* outfile = new TFile(Form("root_files/taushort_invcdf_%dbins.root",nbins),"RECREATE");
    TTree* invcdf = new TTree("expo_taushort","Tree"); 

    
    double tau_chi, tau_chi_err;      // chi square method
    double tau_chip, tau_chip_err;    // Pearson chi square
    double tau_LL, tau_LL_err;        // Log_Likelihood
    double delta_chi, delta_chip, delta_L;


    invcdf->Branch("tau_chi",&tau_chi,"tau_chi/D");
    invcdf->Branch("tau_chi_err",&tau_chi_err,"tau_chi_err/D");
    invcdf->Branch("tau_chip",&tau_chip,"tau_chip/D");
    invcdf->Branch("tau_chip_err",&tau_chip_err,"tau_chip_err/D");
    invcdf->Branch("tau_LL",&tau_LL,"tau_LL/D");
    invcdf->Branch("tau_LL_err",&tau_LL_err,"tau_LL_err/D");
    invcdf->Branch("delta_chi",&delta_chi,"delta_chi/D");
    invcdf->Branch("delta_chip",&delta_chip,"delta_chip/D");
    invcdf->Branch("delta_L",&delta_L,"delta_L/D");

    
    TRandom3* r = new TRandom3;

    progressbar bar(N_sim); 

    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");

    //TCanvas* c = new TCanvas("c", "Simulations",1200,600);

    TF1* f = new TF1("f","expo",Begin,End);

    double y, t;


    for (int j = 0; j < N_sim; j++)
        {
            bar.update();

            TH1D* h = new TH1D("exp","exp",nbins,Begin,End);

            r->SetSeed(j+1);

            for (int k = 0; k < N_events; k++ ) 
            {   
                y = r->Uniform(0,1);
                t = - tau_short * log(1-y);
                h->Fill(t);
            }

            //h->Draw();

            // ---- CHI SQUARE METHOD ---- //
            
            TFitResultPtr tauPtr1 = h->Fit(f, "RNSQ");;
            tau_chi = -1/tauPtr1->Parameter(1);
            tau_chi_err = tauPtr1->ParError(1)/pow(tauPtr1->Parameter(1),2);
            delta_chi = abs(tau_chi - tau_short) / tau_chi_err;
            

             // ---- PEARSON CHI SQUARE METHOD ---- //

            TFitResultPtr tauPtr2 = h->Fit(f, "RPSNQ");;
            tau_chip = -1/tauPtr2->Parameter(1);
            tau_chip_err = tauPtr2->ParError(1)/pow(tauPtr2->Parameter(1),2);;
            delta_chip = abs(tau_chip - tau_short) / tau_chip_err;


            // ---- LOG LIKELIHOOD METHOD ---- //

            TFitResultPtr tauPtr3 = h->Fit(f, "RLSNQ");;
            tau_LL = -1/tauPtr3->Parameter(1);
            tau_LL_err = tauPtr3->ParError(1)/pow(tauPtr3->Parameter(1),2);;
            delta_L = abs(tau_LL - tau_short) / tau_LL_err;

            invcdf->Fill();

        
            h->Delete();

        }

    cout << endl;

    invcdf->Write();   

}