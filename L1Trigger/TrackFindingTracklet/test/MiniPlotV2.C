// ----------------------------------------------------------------------------------------------------------------
// Basic example ROOT script for making quick plots of track distributions using the ntuples produced by L1TrackNtupleMaker.cc
//
// e.g. in ROOT do: .L MiniPlot.C++, MiniPlot("TTbar_PU200_D49")
//
// By Louise Skinnari
// ----------------------------------------------------------------------------------------------------------------

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include <TError.h>
#include <TSystem.h>
#include <TArrayI.h>
#include <TClassTable.h>
#include "TLine.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

float PI = 3.14159265358979323846;

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x, Double_t y, Double_t tsize, Color_t color, char* text);

// ----------------------------------------------------------------------------------------------------------------
// comparing tracks from two different root files, overlaid by seed
// ----------------------------------------------------------------------------------------------------------------

void trk_overlay(TString file1, TString param, int seedNo, double *xRange, int const bins, TString file2 = "0"){

  // file1 = first filename w/o .root
  // file2 = optional, second filename w/o .root

  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;

  SetPlotStyle();
  gStyle->SetOptTitle(1);
  gStyle->SetTitleY(0.97); //y2
  gStyle->SetTitleX(0.60); //x2

  // ----------------------------------------------------------------------------------------------------------------
  // read ntuples
  TString dir = "../src/CMSSW_12_0_0_pre4/src/L1Trigger/TrackFindingTracklet/test/";
  TString filename = dir + file1 + ".root";

  if( gROOT->FindObject(filename.Data()) )
    ((TFile *)gROOT->FindObject(filename.Data()))->Close(); //file location

  TFile *file = TFile::Open(filename.Data());

  file->cd("L1TrackNtuple"); //go into folder

  //add file check

  TTree *eventTree = (TTree *) gROOT->FindObject("eventTree"); //finds trk_phi TTree

    if( gROOT->FindObject("c1") )
      ((TCanvas *)gROOT->FindObject("c1"))->Close();

    TCanvas *c1 = new TCanvas("c1","c1",500,500);
    c1->SetLeftMargin(0.25);
    c1->SetBottomMargin(0.15);
    c1->SetTopMargin(0.1);

    //first track
    if(seedNo == -1){ //show all seeds
      eventTree->Draw(Form("eventTree.%s>>htmp1(%d, 0., 0.)", param.Data(), bins)); // xRange[0], xRange[1]));
    }
    else if(seedNo > 0){ //show particular seed
      // eventTree->Draw(Form("eventTree.%s>>htmp1(%d, 0., 0.)", param.Data(), bins), Form("eventTree.trk_seed == %i", seedNo));
      eventTree->Draw(Form("eventTree.%s>>htmp1(%d, 0., 0.)", param.Data(), bins), Form("eventTree.trk_seed == %i && abs(eventTree.trk_d0) > 4", seedNo));
    }
    TH1F *htmp1 = ((TH1F *)gROOT->FindObject("htmp1"));
    //file->Close();

    //track overlay
    if(file2 != "0"){
      dir = "../src/CMSSW_12_0_0_pre4/src/L1Trigger/TrackFindingTracklet/test/";
      filename = dir + file2 + ".root";

      if( gROOT->FindObject(filename.Data()) )
        ((TFile *)gROOT->FindObject(filename.Data()))->Close(); //file location

      TFile *file = TFile::Open(filename.Data());

      file->cd("L1TrackNtuple"); //go into folder

      TTree *eventTree1 = (TTree *) gROOT->FindObject("eventTree"); //finds trk_phi TTree

      //first track
      if(seedNo == -1){ //show all seeds
        eventTree1->Draw(Form("eventTree.%s>>htmp2(%d, 0., 0.)", param.Data(), bins)); // xRange[0], xRange[1]));
      }
      else if(seedNo > 0){ //show particular seed
        // eventTree1->Draw(Form("eventTree.%s>>htmp2(%d, 0., 0.)", param.Data(), bins), Form("eventTree.trk_seed == %i", seedNo));
        eventTree1->Draw(Form("eventTree.%s>>htmp2(%d, 0., 0.)", param.Data(), bins), Form("eventTree.trk_seed == %i && abs(eventTree.trk_d0) > 4", seedNo));
      }
  }


    TH1F *htmp2 = ((TH1F *)gROOT->FindObject("htmp2"));

    htmp1->SetLineColor(1);
    htmp1->SetMarkerColor(1);
    htmp1->SetMarkerStyle(8);
    htmp1->GetYaxis()->SetRangeUser(0, 180);
    htmp1->Draw("");
    htmp1->SetTitle(Form("seed %i, |d_0| > 4", seedNo));
    htmp1->SetXTitle(param);
    htmp1->SetYTitle(Form("# of Tracks / %f", (xRange[1]-xRange[0])/bins));
    htmp1->GetYaxis()->SetTitleOffset(1.7);

    htmp2->SetLineColor(2);
    htmp2->SetMarkerColor(2);
    htmp2->SetMarkerStyle(24);
    htmp2->Draw(",same");
    htmp2->SetTitle(Form("seed %i, |d_0| > 4", seedNo));
    htmp2->SetXTitle(param);
    htmp2->SetYTitle(Form("# of Tracks / %f", (xRange[1]-xRange[0])/bins));
    htmp2->GetYaxis()->SetTitleOffset(1.7);

    TLegend* l;
    if (param.Contains("trk")) {
      l = new TLegend(0.55,0.22,0.75,0.32);
      mySmallText(0.25,0.22,0.025, 1, "Disp muon, PU=0");
    }
    l->SetFillColor(0);
    l->SetLineColor(0);
    l->SetTextSize(0.02);
    l->AddEntry(htmp1,"approx","ep");
    l->AddEntry(htmp2,"no approx","ep");
    l->SetTextFont(42);
    l->Draw();

    c1->Draw();
}

// *** visualizing the projected stubs found in layers/disks per seed ***
void stub_plot(TString file1, TString file2, int seedNo, TString param){

  // gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;
  SetPlotStyle();
  gStyle->SetOptTitle(1);

    double xRange[2] = {0, 11};
    int bins = 11;

    gSystem->cd("/Users/justinvega/Documents/CERN/src/CMSSW_12_0_0_pre4/src/L1Trigger/TrackFindingTracklet/test/");

    TH1F* h_trk_tracklet_hits_pos = new TH1F("trklt_hits_pos", Form("Stubs in Seed %i, #eta > 0, |d0| > 4;layer/disk; # of Tracks", seedNo), bins, xRange[0], xRange[1]);
    TH1F* h_trk_tracklet_hits_pos1 = new TH1F("trklt_hits_pos1", Form("Stubs in Seed %i, #eta > 0, |d0| > 4;layer/disk; # of Tracks", seedNo), bins, xRange[0], xRange[1]);

    TH1F* h_trk_tracklet_hits_neg = new TH1F("trklt_hits_neg", Form("Stubs in Seed %i, #eta < 0, |d0| > 4;layer/disk; # of Tracks", seedNo), bins, xRange[0], xRange[1]);
    TH1F* h_trk_tracklet_hits_neg1 = new TH1F("trklt_hits_neg1", Form("Stubs in Seed %i, #eta < 0, |d0| > 4;layer/disk; # of Tracks", seedNo), bins, xRange[0], xRange[1]);

    TString filelist[2] = {file1, file2};
    std::string pos_histo[2] = {"trklt_hits_pos", "trklt_hits_pos1"};
    std::string neg_histo[2] = {"trklt_hits_neg", "trklt_hits_neg1"};

    for (int fileidx = 0; fileidx < sizeof(filelist)/sizeof(filelist[0]); fileidx++){

      // ----------------------------------------------------------------------------------------------------------------
      // read ntuples
      //TString dir = "../src/CMSSW_12_0_0_pre4/src/L1Trigger/TrackFindingTracklet/test/";
      TString filename = filelist[fileidx] + ".root";

      if( gROOT->FindObject(filename.Data()) )
        ((TFile *)gROOT->FindObject(filename.Data()))->Close(); //file location

      TFile *file = TFile::Open(filename.Data());
      file->cd("L1TrackNtuple"); //go into folder

      //add file check

      TTree *tree = (TTree *) gROOT->FindObject("eventTree"); //finds trk_phi TTree

        if( gROOT->FindObject("c1") )
          ((TCanvas *)gROOT->FindObject("c1"))->Close();

        //TCanvas *c1 = new TCanvas("c1","c1",600,600);
        vector<float>* trk_eta;
        vector<float>* trk_phi;
        vector<float>* trk_invr;
        vector<int>* trk_lhits;
        vector<int>* trk_dhits;
        vector<int>* trk_seed;
        vector<int>* trk_d0;

        trk_eta = 0;
        trk_phi = 0;
        trk_lhits = 0;
        trk_dhits = 0;
        trk_seed = 0;
        trk_d0 = 0;
        trk_invr = 0;

        tree->SetBranchAddress("trk_eta", &trk_eta); //name of leaf, initial value, branch
        tree->SetBranchAddress("trk_phi", &trk_phi);
        tree->SetBranchAddress("trk_d0", &trk_d0);
        tree->SetBranchAddress("trk_lhits", &trk_lhits);
        tree->SetBranchAddress("trk_dhits", &trk_dhits);
        tree->SetBranchAddress("trk_seed", &trk_seed);
        tree->SetBranchAddress("trk_invr", &trk_invr);

          int nevt = tree->GetEntries();
          // cout << "number of events = " << nevt << endl;

          for (int i = 0; i < nevt; i++) {
            tree->GetEntry(i, 0);

            // ----------------------------------------------------------------------------------------------------------------
            // loop over *all* found found L1 tracks
            for (int it = 0; it < (int)trk_lhits->size(); it++) {

              // here you can place some selection criteria and fill histograms for the L1 tracks

              // create an 11-bit long iterable from lhits and dhits
              int num_layers = 6;
              int num_discs = 5;
              int lhits = trk_lhits->at(it);
              int dhits = trk_dhits->at(it);
              int eta = trk_eta->at(it);
              int d0 = trk_d0->at(it);
              int seed = trk_seed->at(it);
              // float rInv = trk_invr->at(it);
              std::vector<int> layers = {};

              for (int layer_index = 0; layer_index < num_layers + num_discs; layer_index++) {
                if (layer_index < num_layers) {
                  layers.push_back(lhits % 10);
                  lhits /= 10;
                } else {
                  layers.push_back(dhits % 10);
                  dhits /= 10;
                }
              }

              for (unsigned int layer = 0; layer < layers.size(); layer++){
                if (layers.at(layer)) {
                  if(eta<0){
                    if(abs(d0)>4){
                      //if(rInv > 0){
                        if(seed == seedNo){
                          TH1F* htmp = (TH1F *) gROOT->FindObjectAny(neg_histo[fileidx].c_str());  // if there was a hit at this layer...
                          htmp->Fill(layer);  // ...fill this bin with the layer of the track.
                        }
                      //}
                    }
                  }
                  else if(eta>0){
                    if(abs(d0)>4){
                      if(seed == seedNo){
                        TH1F* htmp = (TH1F *) gROOT->FindObjectAny(pos_histo[fileidx].c_str());
                        htmp->Fill(layer);
                      }
                    }
                  }
                }
              }

            }
          }

        }
        TCanvas *c1 = new TCanvas("c1","c1",500,500);
        c1->SetLeftMargin(0.25);
        c1->SetBottomMargin(0.15);
        c1->SetTopMargin(0.1);
        gStyle->SetTitleY(0.97); //y2
        gStyle->SetTitleX(0.60); //x2
        // gStyle->SetTitleOffset(1);

        // determine max bin for y-range
        int negmaxbin = (h_trk_tracklet_hits_neg > h_trk_tracklet_hits_neg1) ? h_trk_tracklet_hits_neg->GetBinContent(h_trk_tracklet_hits_neg->GetMaximumBin()) : h_trk_tracklet_hits_neg1->GetBinContent(h_trk_tracklet_hits_neg1->GetMaximumBin());
        int posmaxbin = (h_trk_tracklet_hits_pos > h_trk_tracklet_hits_pos1) ? h_trk_tracklet_hits_pos->GetBinContent(h_trk_tracklet_hits_pos->GetMaximumBin()) : h_trk_tracklet_hits_pos1->GetBinContent(h_trk_tracklet_hits_pos1->GetMaximumBin());
        int negentry = h_trk_tracklet_hits_neg->GetEntries();
        int neg1entry = h_trk_tracklet_hits_neg1->GetEntries();
        int posentry = h_trk_tracklet_hits_pos->GetEntries();
        int pos1entry = h_trk_tracklet_hits_pos1->GetEntries();

        cout << Form("%i, %i, %i, %i", negentry, neg1entry, posentry, pos1entry) << endl;

        int pad = 0;
        if(seedNo == 9){
          pad = 50;
        }
        else {
          pad = 400;
        }

        cout << Form("%i, %i", negmaxbin, posmaxbin) << endl;

        if (param == "neg"){
          h_trk_tracklet_hits_neg1->SetLineColor(1);
          h_trk_tracklet_hits_neg1->SetMarkerColor(1);
          h_trk_tracklet_hits_neg1->SetMarkerStyle(24);
          h_trk_tracklet_hits_neg1->GetXaxis()->SetTitleOffset(1.3);
          h_trk_tracklet_hits_neg1->GetXaxis()->CenterTitle();
          h_trk_tracklet_hits_neg1->GetYaxis()->SetTitleOffset(1.7);
          h_trk_tracklet_hits_neg1->GetYaxis()->CenterTitle();
          h_trk_tracklet_hits_neg1->GetYaxis()->SetRangeUser(0, negmaxbin + pad);
          h_trk_tracklet_hits_neg1->Draw();
          c1->Update();
          h_trk_tracklet_hits_neg->SetLineColor(2);
          h_trk_tracklet_hits_neg->SetMarkerColor(2);
          h_trk_tracklet_hits_neg->SetMarkerStyle(8);
          h_trk_tracklet_hits_neg->GetYaxis()->SetRangeUser(0, negmaxbin + pad);
          h_trk_tracklet_hits_neg->GetXaxis()->SetTitleOffset(1.3);
          h_trk_tracklet_hits_neg->GetXaxis()->CenterTitle();
          h_trk_tracklet_hits_neg->GetYaxis()->SetTitleOffset(1.7);
          h_trk_tracklet_hits_neg->GetYaxis()->CenterTitle();
          h_trk_tracklet_hits_neg->Draw(",same");

          TLegend* l;
            l = new TLegend(0.75,0.80,0.95,0.90);
            mySmallText(0.55,0.82,0.025, 1, "Disp muon, PU=0");
          l->SetFillColor(0);
          l->SetLineColor(1);
          l->SetTextSize(0.02);
          l->AddEntry(h_trk_tracklet_hits_neg,"approx","l");
          l->AddEntry(h_trk_tracklet_hits_neg1,"no approx","l");
          l->SetTextFont(42);
          l->Draw();
            }

          else if (param == "pos"){
            h_trk_tracklet_hits_pos1->SetTitleOffset(1);
            h_trk_tracklet_hits_pos1->SetLineColor(1);
            h_trk_tracklet_hits_pos1->SetMarkerColor(1);
            h_trk_tracklet_hits_pos1->SetMarkerStyle(24);
            h_trk_tracklet_hits_pos1->GetXaxis()->SetTitleOffset(1.3);
            h_trk_tracklet_hits_pos1->GetXaxis()->CenterTitle();
            h_trk_tracklet_hits_pos1->GetYaxis()->SetTitleOffset(1.7);
            h_trk_tracklet_hits_pos1->GetYaxis()->CenterTitle();
            h_trk_tracklet_hits_pos1->GetYaxis()->SetRangeUser(0, posmaxbin + pad);
            h_trk_tracklet_hits_pos1->Draw();
            // c1->Update();
            // h_trk_tracklet_hits_pos->SetLineColor(2);
            // h_trk_tracklet_hits_pos->SetMarkerColor(2);
            // h_trk_tracklet_hits_pos->SetMarkerStyle(8);
            // h_trk_tracklet_hits_pos->GetYaxis()->SetRangeUser(0, posmaxbin + pad);
            // h_trk_tracklet_hits_pos->GetXaxis()->SetTitleOffset(1.3);
            // h_trk_tracklet_hits_pos->GetXaxis()->CenterTitle();
            // h_trk_tracklet_hits_pos->GetYaxis()->SetTitleOffset(1.7);
            // h_trk_tracklet_hits_pos->GetYaxis()->CenterTitle();
            // h_trk_tracklet_hits_pos->Draw(",same");
            //
            // TLegend* l;
            //   l = new TLegend(0.75,0.80,0.95,0.90);
            //   mySmallText(0.55,0.84,0.025, 1, "Disp muon, PU=0");
            // l->SetFillColor(0);
            // l->SetLineColor(1);
            // l->SetTextSize(0.02);
            // l->AddEntry(h_trk_tracklet_hits_pos,"approx","l");
            // l->AddEntry(h_trk_tracklet_hits_pos1,"no approx","l");
            // l->SetTextFont(42);
            // l->Draw();
          }

          // bin labels
            std::vector<string> SEED = {"L1", "L2", "L3", "L4", "L5", "L6", "D1", "D2", "D3", "D4", "D5"};

            for(int r=0; r< SEED.size(); r++){
              h_trk_tracklet_hits_neg->GetXaxis()->SetBinLabel(r+1, SEED[r].c_str());
              h_trk_tracklet_hits_pos->GetXaxis()->SetBinLabel(r+1, SEED[r].c_str());
              h_trk_tracklet_hits_neg1->GetXaxis()->SetBinLabel(r+1, SEED[r].c_str());
              h_trk_tracklet_hits_pos1->GetXaxis()->SetBinLabel(r+1, SEED[r].c_str());
            }

            h_trk_tracklet_hits_neg->GetXaxis()->SetLabelSize(12);
            h_trk_tracklet_hits_neg->GetYaxis()->SetLabelSize(16);
            h_trk_tracklet_hits_pos->GetXaxis()->SetLabelSize(12);
            h_trk_tracklet_hits_pos->GetYaxis()->SetLabelSize(16);
            h_trk_tracklet_hits_neg1->GetXaxis()->SetLabelSize(12);
            h_trk_tracklet_hits_neg1->GetYaxis()->SetLabelSize(16);
            h_trk_tracklet_hits_pos1->GetXaxis()->SetLabelSize(12);
            h_trk_tracklet_hits_pos1->GetYaxis()->SetLabelSize(16);

            c1->Draw();
            // c1->SaveAs("seedplots/stubs_layerdisk_seed8.pdf");

}
// *** with charge separation/filtering ***
void stub_plotV2(TString file1, TString file2, int seedNo){

  // gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;
  SetPlotStyle();
  gStyle->SetOptTitle(1);

    double xRange[2] = {0, 11};
    int bins = 11;

    gSystem->cd("/Users/justinvega/Documents/CERN/src/CMSSW_12_0_0_pre4/src/L1Trigger/TrackFindingTracklet/test/");

    TH1F* h_trk_tracklet_hits_poseta[5]; //0 = (+, file1), 1 = (-, file1), etc.
    TH1F* h_trk_tracklet_hits_negeta[5]; //0 = (+, file1), 1 = (-, file1), etc.

    for(int j = 0; j < 4; j++){
      h_trk_tracklet_hits_poseta[j] = new TH1F(Form("trklt_hits_poseta%i", j), Form("Stubs in Seed %i, #eta > 0, |d0| > 4;layer/disk; # of Tracks", seedNo), bins, xRange[0], xRange[1]);
      h_trk_tracklet_hits_negeta[j] = new TH1F(Form("trklt_hits_negeta%i", j), Form("Stubs in Seed %i, #eta < 0, |d0| > 4;layer/disk; # of Tracks", seedNo), bins, xRange[0], xRange[1]);
    }

    TString filelist[2] = {file1, file2};
    int qpos_idx[2] = {0, 2};
    int qneg_idx[2] = {1, 3};

    for(int fileidx = 0; fileidx < sizeof(filelist)/sizeof(filelist[0]); fileidx++){
      // ----------------------------------------------------------------------------------------------------------------
      // read ntuples
      //TString dir = "../src/CMSSW_12_0_0_pre4/src/L1Trigger/TrackFindingTracklet/test/";
      TString filename = filelist[fileidx] + ".root";

      if( gROOT->FindObject(filename.Data()) )
        ((TFile *)gROOT->FindObject(filename.Data()))->Close(); //file location

      TFile *file = TFile::Open(filename.Data());
      file->cd("L1TrackNtuple"); //go into folder

      //add file check?

      TTree *tree = (TTree *) gROOT->FindObject("eventTree"); //finds trk_phi TTree

        if( gROOT->FindObject("c1") )
          ((TCanvas *)gROOT->FindObject("c1"))->Close();

        // initializing branches
        vector<float>* trk_eta;
        vector<float>* trk_phi;
        vector<float>* trk_invr;
        vector<int>* trk_lhits;
        vector<int>* trk_dhits;
        vector<int>* trk_seed;
        vector<int>* trk_d0;

        trk_eta = 0;
        trk_phi = 0;
        trk_lhits = 0;
        trk_dhits = 0;
        trk_seed = 0;
        trk_d0 = 0;
        trk_invr = 0;

        tree->SetBranchAddress("trk_eta", &trk_eta); //name of leaf, initial value, branch
        tree->SetBranchAddress("trk_phi", &trk_phi);
        tree->SetBranchAddress("trk_d0", &trk_d0);
        tree->SetBranchAddress("trk_lhits", &trk_lhits);
        tree->SetBranchAddress("trk_dhits", &trk_dhits);
        tree->SetBranchAddress("trk_seed", &trk_seed);
        tree->SetBranchAddress("trk_invr", &trk_invr);

          int nevt = tree->GetEntries();
          // cout << "number of events = " << nevt << endl;

          for (int i = 0; i < nevt; i++) { //for each event
            tree->GetEntry(i, 0);

            // ----------------------------------------------------------------------------------------------------------------
            // loop over *all* found found L1 tracks
            for (int it = 0; it < (int)trk_lhits->size(); it++) {

              // here you can place some selection criteria and fill histograms for the L1 tracks

              // create an 11-bit long iterable from lhits and dhits
              int num_layers = 6;
              int num_discs = 5;
              int lhits = trk_lhits->at(it);
              int dhits = trk_dhits->at(it);
              int eta = trk_eta->at(it);
              int d0 = trk_d0->at(it);
              int seed = trk_seed->at(it);
              float rInv = trk_invr->at(it);
              std::vector<int> layers = {};

              for (int layer_index = 0; layer_index < num_layers + num_discs; layer_index++) {
                if (layer_index < num_layers) {
                  layers.push_back(lhits % 10);
                  lhits /= 10;
                }
                else {
                  layers.push_back(dhits % 10);
                  dhits /= 10;
                }
              }

              for (unsigned int layer = 0; layer < layers.size(); layer++){
                if (layers.at(layer)) {
                  if(eta<0){ //negative eta
                    if(abs(d0)>4){
                      if(seed == seedNo){
                        if(rInv > 0){ //positive charge (0, 2)
                          //find + charge, negative eta (h_trk_tracklet_hits_negeta[0] or [2])
                          TString histoName= Form("trklt_hits_negeta%i", qpos_idx[fileidx]);
                          TH1F* htmp = (TH1F *) gROOT->FindObjectAny(histoName.Data());  // if there was a hit at this layer...
                          htmp->Fill(layer);  // ...fill this bin with the layer of the track.
                        }
                        else if (rInv < 0){ //negative charge (1, 3)
                          //find - charge, negative eta (h_trk_tracklet_hits_negeta[1] or [3])
                          TString histoName = Form("trklt_hits_negeta%i", qneg_idx[fileidx]);
                          TH1F* htmp = (TH1F *) gROOT->FindObjectAny(histoName.Data());  // if there was a hit at this layer...
                          htmp->Fill(layer);  // ...fill this bin with the layer of the track.
                        } //"trklt_hits_poseta%i", trklt_hits_negeta%i
                      }
                    }
                  }
                  else if(eta>0){ //positive eta
                    if(abs(d0)>4){
                      if(seed == seedNo){
                        if(rInv > 0){ //positive charge (0, 2)
                          //find + charge, positive eta (h_trk_tracklet_hits_poseta[0] or [2])
                          TString histoName = Form("trklt_hits_poseta%i", qpos_idx[fileidx]);
                          TH1F* htmp = (TH1F *) gROOT->FindObjectAny(histoName.Data());  // if there was a hit at this layer...
                          htmp->Fill(layer);  // ...fill this bin with the layer of the track.
                        }
                        else if (rInv < 0){ //negative charge (1, 3)
                          //find - charge, positive eta (h_trk_tracklet_hits_poseta[1] or [3])
                          TString histoName = Form("trklt_hits_poseta%i", qneg_idx[fileidx]);
                          TH1F* htmp = (TH1F *) gROOT->FindObjectAny(histoName.Data());  // if there was a hit at this layer...
                          htmp->Fill(layer);  // ...fill this bin with the layer of the track.
                        }
                      }

                    }
                  }
                }
              }

            }
          }
        }

        //canvas set-up
        TCanvas *c2 = new TCanvas("c2","c2",500,500);
        c2->SetLeftMargin(0.25);
        c2->SetBottomMargin(0.15);
        c2->SetTopMargin(0.1);
        gStyle->SetTitleY(0.97); //y2
        gStyle->SetTitleX(0.60); //x2

          for(int chargeidx=0; chargeidx < 2; chargeidx++){
            //set conditions for approxtrue/file1 (0, 1)
            h_trk_tracklet_hits_negeta[chargeidx]->SetLineColor(2);
            h_trk_tracklet_hits_negeta[chargeidx]->SetMarkerColor(2);
            h_trk_tracklet_hits_negeta[chargeidx]->SetMarkerStyle(24);
            h_trk_tracklet_hits_poseta[chargeidx]->SetLineColor(2);
            h_trk_tracklet_hits_poseta[chargeidx]->SetMarkerColor(2);
            h_trk_tracklet_hits_poseta[chargeidx]->SetMarkerStyle(24);
          }

          for(int charge1idx=2; charge1idx < 4; charge1idx++){
            //set conditions for approxfalse/file2 (2, 3)
            h_trk_tracklet_hits_negeta[charge1idx]->SetLineColor(1);
            h_trk_tracklet_hits_negeta[charge1idx]->SetMarkerColor(1);
            h_trk_tracklet_hits_negeta[charge1idx]->SetMarkerStyle(8);
            // h_trk_tracklet_hits_negeta[charge1idx]->GetYaxis()->SetRangeUser(0, negmaxbin + pad);
            h_trk_tracklet_hits_poseta[charge1idx]->SetLineColor(1);
            h_trk_tracklet_hits_poseta[charge1idx]->SetMarkerColor(1);
            h_trk_tracklet_hits_poseta[charge1idx]->SetMarkerStyle(8);
            // h_trk_tracklet_hits_poseta[charge1idx]->GetYaxis()->SetRangeUser(0, posmaxbin + pad);
          }

          // bin labels
          std::vector<string> SEED = {"L1", "L2", "L3", "L4", "L5", "L6", "D1", "D2", "D3", "D4", "D5"};

          for(int allplots=0; allplots < 4; allplots++){
            for(int r=0; r< SEED.size(); r++){
              h_trk_tracklet_hits_negeta[allplots]->GetXaxis()->SetBinLabel(r+1, SEED[r].c_str());
              h_trk_tracklet_hits_negeta[allplots]->GetXaxis()->SetLabelSize(12);
              //h_trk_tracklet_hits_neg->GetYaxis()->SetLabelSize(16);
              h_trk_tracklet_hits_poseta[allplots]->GetXaxis()->SetBinLabel(r+1, SEED[r].c_str());
              h_trk_tracklet_hits_poseta[allplots]->GetXaxis()->SetLabelSize(12);
            }
            if(seedNo==11){
              h_trk_tracklet_hits_poseta[allplots]->GetYaxis()->SetRangeUser(0., 600.);
              h_trk_tracklet_hits_negeta[allplots]->GetYaxis()->SetRangeUser(0., 600.);
            }
          }


          //for negative eta, positive charge
          THStack *hs = new THStack("hs",Form("Stubs in Seed %i, #eta < 0, |d0| > 4, + charge;layer/disk; # of Tracks", seedNo));
          hs->Add(h_trk_tracklet_hits_negeta[2]);
          hs->Add(h_trk_tracklet_hits_negeta[0]);
          hs->Draw("nostack");

          TLegend* l;
            l = new TLegend(0.75,0.80,0.95,0.90);
            mySmallText(0.35,0.05,0.025, 1, "Disp muon, PU=0");
          l->SetFillColor(0);
          l->SetLineColor(1);
          l->SetTextSize(0.02);
          l->SetTextFont(42);
          l->AddEntry(h_trk_tracklet_hits_negeta[0],"approx","l");
          l->AddEntry(h_trk_tracklet_hits_negeta[2],"no approx","l");
          l->Draw();
          c2->Draw();
          c2->SaveAs(Form("seedplots/stubs_layerdisk_seed%i_negeta_qpos.pdf", seedNo));
          l->Clear();

          //for negative eta, negative charge
          THStack *hs1 = new THStack("hs",Form("Stubs in Seed %i, #eta < 0, |d0| > 4, - charge;layer/disk; # of Tracks", seedNo));
          hs1->Add(h_trk_tracklet_hits_negeta[3]);
          hs1->Add(h_trk_tracklet_hits_negeta[1]);
          hs1->Draw("nostack");

          l->AddEntry(h_trk_tracklet_hits_negeta[1],"approx","l");
          l->AddEntry(h_trk_tracklet_hits_negeta[3],"no approx","l");
          mySmallText(0.35,0.05,0.025, 1, "Disp muon, PU=0");
          l->Draw();

          c2->SaveAs(Form("seedplots/stubs_layerdisk_seed%i_negeta_qneg.pdf", seedNo));
          l->Clear();

          //for positive eta, positive charge
          THStack *hs2 = new THStack("hs",Form("Stubs in Seed %i, #eta > 0, |d0| > 4, + charge;layer/disk; # of Tracks", seedNo));
          hs2->Add(h_trk_tracklet_hits_poseta[2]);
          hs2->Add(h_trk_tracklet_hits_poseta[0]);
          hs2->Draw("nostack");

          l->AddEntry(h_trk_tracklet_hits_poseta[0],"approx","l");
          l->AddEntry(h_trk_tracklet_hits_poseta[2],"no approx","l");
          mySmallText(0.35,0.05,0.025, 1, "Disp muon, PU=0");
          l->Draw();

          c2->SaveAs(Form("seedplots/stubs_layerdisk_seed%i_poseta_qpos.pdf", seedNo));
          l->Clear();

          //for positive eta, negative charge
          THStack *hs3 = new THStack("hs",Form("Stubs in Seed %i, #eta > 0, |d0| > 4, - charge;layer/disk; # of Tracks", seedNo));
          hs3->Add(h_trk_tracklet_hits_poseta[3]);
          hs3->Add(h_trk_tracklet_hits_poseta[1]);
          hs3->Draw("nostack");

          l->AddEntry(h_trk_tracklet_hits_poseta[1],"approx","l");
          l->AddEntry(h_trk_tracklet_hits_poseta[3],"no approx","l");
          mySmallText(0.35,0.05,0.025, 1, "Disp muon, PU=0");
          l->Draw();

          c2->SaveAs(Form("seedplots/stubs_layerdisk_seed%i_poseta_qneg.pdf", seedNo));

}

void SetPlotStyle() {

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20, 26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.3);
  gStyle->SetTitleYOffset(1.7);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42, "x");
  gStyle->SetTitleFont(42, "x");
  gStyle->SetLabelFont(42, "y");
  gStyle->SetTitleFont(42, "y");
  gStyle->SetLabelFont(42, "z");
  gStyle->SetTitleFont(42, "z");
  gStyle->SetLabelSize(0.03, "x");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetLabelSize(0.04, "y");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetLabelSize(0.05, "z");
  gStyle->SetTitleSize(0.05, "z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2, "[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

void mySmallText(Double_t x, Double_t y, Double_t tsize, Color_t color, char* text) {
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}
