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


#include "progressbar.hpp"


void baselineMC(int N_sim, int range = 1, int weeks = 1)
{ 
    const double baseline = weeks * 4.;

    double Begin = 179.5;
    double End = 2048.5;

    double Begin_fit, End_fit;

    int nbins, N_entries;

    if(range == 1) 
    {
        Begin_fit = 179.5;
        End_fit = 2048.5;
    }

    if(range == 2) 
    {
        Begin_fit = 1300.;
        End_fit = 2048.;
    }


    nbins = (End - Begin);
    N_entries = baseline * (End - Begin);

    cout << "bins: " << nbins << endl;
    cout << "Nentries " << N_entries << endl;

    TFile* outfile = new TFile(Form("root_files/baseline_%dw_range%d.root",weeks,range),"RECREATE");
    TTree* outtree = new TTree(Form("bs_%dweek",weeks),"Tree"); 

    cout << "writing file: " << Form("baseline_%dw_range%d.root",weeks,range) << endl;
    

    double bs_mean, bs_mean_err;    // mean value
    double bs_chi, bs_chi_err;      // chi square method
    double bs_chip, bs_chip_err;    // Pearson chi square
    double bs_LL, bs_LL_err;        // Log_Likelihood
    double delta_chi, delta_chip, delta_mean, delta_L;


    outtree->Branch("bs_mean",&bs_mean,"bs_mean/D");
    outtree->Branch("bs_mean_err",&bs_mean_err,"bs_mean_err/D");
    outtree->Branch("bs_chi",&bs_chi,"bs_chi/D");
    outtree->Branch("bs_chi_err",&bs_chi_err,"bs_chi_err/D");
    outtree->Branch("bs_chip",&bs_chip,"bs_chip/D");
    outtree->Branch("bs_chip_err",&bs_chip_err,"bs_chip_err/D");
    outtree->Branch("bs_LL",&bs_LL,"bs_LL/D");
    outtree->Branch("bs_LL_err",&bs_LL_err,"bs_LL_err/D");
    outtree->Branch("delta_mean",&delta_mean,"delta_mean/D");
    outtree->Branch("delta_chi",&delta_chi,"delta_chi/D");
    outtree->Branch("delta_chip",&delta_chip,"delta_chip/D");
    outtree->Branch("delta_L",&delta_L,"delta_L/D");
    
    //TCanvas* c = new TCanvas("c", "Simulations",1200,600);

    TRandom3* r = new TRandom3;


    progressbar bar(N_sim); 

    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");

    double sum, mean, errmean, sigma;   
    TF1* f0 = new TF1("f0","pol0",Begin_fit,End_fit);

    for (int j = 0; j < N_sim; j++)
        {

        bar.update();

        sigma = 0;
        sum = 0;

        TH1D* h = new TH1D("baseline","baseline",nbins,Begin,End);

        r->SetSeed(j+1);

        for (int k = 0; k < N_entries; k++ ) 
        {   
            h->Fill(r->Uniform(Begin,End));
        }


        if(range == 1) 
        {
            for(int i = 1; i <= nbins; i++)
            {
                sum += h->GetBinContent(i);
            }

            mean = sum / (nbins);
        }

        if(range == 2) 
        {
            for(int i = 1300; i <= nbins; i++)
            {
                sum += h->GetBinContent(i);
            }

            mean = sum / (nbins-1300);
        }

        

        if(range == 1) 
        {
            for(int i = 1; i <= nbins; i++)

            {
                sigma +=pow((h->GetBinContent(i) - mean),2);
            }

        }
        

        if(range == 2) 
        {
            for(int i = 1300; i <= nbins; i++)

            {
                sigma +=pow((h->GetBinContent(i) - mean),2);
            }

        }

        errmean = sqrt(sigma/sqrt(N_sim-1)) / sqrt(N_sim);
        
        bs_mean = mean;
        bs_mean_err = errmean;
        delta_mean = abs(mean - baseline) / errmean;


        // ---- CHI SQUARE METHOD ---- //

        TFitResultPtr basePtr1 = h->Fit(f0, "RNSQ");;
        bs_chi = basePtr1->Parameter(0);
        bs_chi_err = basePtr1->ParError(0);
        delta_chi = abs(bs_chi - baseline) / bs_chi_err;

        // ---- PEARSON CHI SQUARE METHOD ---- //

        TFitResultPtr basePtr2 = h->Fit(f0, "RPSNQ");;
        bs_chip = basePtr2->Parameter(0);
        bs_chip_err = basePtr2->ParError(0);
        delta_chip = abs(bs_chip - baseline) / bs_chip_err;

        // ---- LOG LIKELIHOOD METHOD ---- //

        TFitResultPtr basePtr3 = h->Fit(f0, "RLSNQ");;
        bs_LL = basePtr3->Parameter(0);
        bs_LL_err = basePtr3->ParError(0);
        delta_L = abs(bs_LL - baseline) / bs_LL_err;




        
        outtree->Fill();

        //h->Draw();
        h->Delete();
        }


    cout << endl;

    outtree->Write();
    outfile->Close();

}

