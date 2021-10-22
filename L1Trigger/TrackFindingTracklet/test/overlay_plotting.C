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

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text);


// ----------------------------------------------------------------------------------------------------------------
// Main script for overlaid histograms, from Louise
void overlay(TString what, TString file1, TString file2) { //'what' is a specific histogram?

  SetPlotStyle();

  if( gROOT->FindObject(file1.Data()) )
    ((TFile *)gROOT->FindObject(file1.Data()))->Close();

  if( gROOT->FindObject(file2.Data()) )
    ((TFile *)gROOT->FindObject(file2.Data()))->Close();

    TString dir = "../src/CMSSW_12_0_0_pre4/src/L1Trigger/TrackFindingTracklet/test/";

  // files to pull from
  TFile* tree1 = TFile::Open(dir + file1.Data());
  TFile* tree2 = TFile::Open(dir + file2.Data());

  TH1F* h1 = (TH1F*) tree1->Get(what);
  TH1F* h2 = (TH1F*) tree2->Get(what);
  //TH1F* h3 = (TH1F*) tree3->Get(what);
  //TH1F* h4 = (TH1F*) tree4->Get(what);

  if (what.Contains("ntrk")) {
    h1->Rebin(10);
    h2->Rebin(10);
    //h3->Rebin(10);
  }

  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  h1->SetLineColor(1);
  h1->SetMarkerColor(1);
  h1->SetMarkerStyle(8);
  h1->Draw(""); //l draws a line, p draws a marker
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->SetMarkerStyle(24);
  h2->Draw(",same");
  // h3->SetLineColor(4);
  // h3->SetMarkerColor(4);
  // h3->SetMarkerStyle(22);
  // h3->Draw("lp,same");

  //h4->SetLineColor(8);
  //h4->SetMarkerColor(8);
  //h4->SetMarkerStyle(23);
  //h4->Draw("lp,same");

  TLegend* l;
  if (what.Contains("eff")) {
    l = new TLegend(0.55,0.22,0.75,0.40);
    mySmallText(0.30,0.22,1, "Disp muon, PU=0");
  }
  else {
    l = new TLegend(0.2,0.72,0.5,0.9);
    mySmallText(0.6,0.82,1, "Disp muon, PU=0");
  }
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h1,"approx","ep");
  l->AddEntry(h2,"no approx","ep");
  // l->AddEntry(h3,"Tracklet sim","lep");
  //l->AddEntry(h4,"Extended","lep");
  l->SetTextFont(42);
  l->Draw();

  gPad->SetGridy();

  c1->Modified();
  c1->Update();
  c1->Draw();
  //TString dir = "/Users/justinvega/Documents/CERN/src/CMSSW_12_0_0_pre4/src/L1Trigger/TrackFindingTracklet/test/";
  c1->SaveAs(dir + "Overlay/disp_muon_overlay_"+what+"_approx.pdf");

}


void SetPlotStyle() {

  // from ATLAS plot style macro

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
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42,"x");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]");

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
// *************
// Justin's script for plotting overlaid ntuples from separate .root files (this one is for tracks specifically)
// *************

void overlay_trk(TString what, TString file1, TString file2, double *xRange, int bins){

  TChain* tree = new TChain("L1TrackNtuple/eventTree");
  TChain* tree1 = new TChain("L1TrackNtuple/eventTree");
  TString dir = "../src/CMSSW_12_0_0_pre4/src/L1Trigger/TrackFindingTracklet/test/";
  tree->Add(dir + file1 + ".root");
  tree1->Add(dir + file2 + ".root");

  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }

  if (tree1->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }

  vector<float>* trk_z0;
  vector<float>* trk_d0;
  vector<float>* trk_pt;
  vector<float>* trk_eta;
  vector<float>* trk_phi;
  vector<float>* trk_z01;
  vector<float>* trk_d01;
  vector<float>* trk_pt1;
  vector<float>* trk_eta1;
  vector<float>* trk_phi1;

  TBranch* b_trk_z0;
  TBranch* b_trk_d0;
  TBranch* b_trk_pt;
  TBranch* b_trk_eta;
  TBranch* b_trk_phi;
  TBranch* b_trk_z01;
  TBranch* b_trk_d01;
  TBranch* b_trk_pt1;
  TBranch* b_trk_eta1;
  TBranch* b_trk_phi1;

  trk_z0 = 0;
  trk_d0 = 0;
  trk_pt = 0;
  trk_eta = 0;
  trk_phi = 0;
  trk_z01 = 0;
  trk_d01 = 0;
  trk_pt1 = 0;
  trk_eta1 = 0;
  trk_phi1 = 0;

  tree->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
  tree->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
  tree->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
  tree->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
  tree->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
  tree1->SetBranchAddress("trk_z0", &trk_z01, &b_trk_z01);
  tree1->SetBranchAddress("trk_d0", &trk_d01, &b_trk_d01);
  tree1->SetBranchAddress("trk_pt", &trk_pt1, &b_trk_pt1);
  tree1->SetBranchAddress("trk_eta", &trk_eta1, &b_trk_eta1);
  tree1->SetBranchAddress("trk_phi", &trk_phi1, &b_trk_phi1);


  TH1F* h_trk_pt = new TH1F("trk_pt", Form(";Track p_{T} (GeV); Tracks / %f GeV", (xRange[1]-xRange[0])/bins), bins, xRange[0], xRange[1]);
  TH1F* h_trk_pt1 = new TH1F("trk_pt1", Form(";Track p_{T} (GeV); Tracks / %f GeV", (xRange[1]-xRange[0])/bins), bins, xRange[0], xRange[1]);

  TH1F* h_trk_eta = new TH1F("trk_eta", Form(";Track eta; Tracks / %f", (xRange[1]-xRange[0])/bins), bins, xRange[0], xRange[1]);
  TH1F* h_trk_eta1 = new TH1F("trk_eta1", Form(";Track eta; Tracks / %f", (xRange[1]-xRange[0])/bins), bins, xRange[0], xRange[1]);

  TH1F* h_trk_phi = new TH1F("trk_phi", Form(";Track phi (rad); Tracks / %f rad", (xRange[1]-xRange[0])/bins), bins, xRange[0], xRange[1]);
  TH1F* h_trk_phi1 = new TH1F("trk_phi1", Form(";Track phi (rad); Tracks / %f rad", (xRange[1]-xRange[0])/bins), bins, xRange[0], xRange[1]);

  TH1F* h_trk_d0 = new TH1F("trk_d0", Form(";Track d0 (cm); Tracks / %f cm", (xRange[1]-xRange[0])/bins), bins, xRange[0], xRange[1]);
  TH1F* h_trk_d01 = new TH1F("trk_d01", Form(";Track d0 (cm); Tracks / %f cm", (xRange[1]-xRange[0])/bins), bins, xRange[0], xRange[1]);

  TH1F* h_trk_z0 = new TH1F("trk_z0", Form(";Track z0 (cm); Tracks / %f cm", (xRange[1]-xRange[0])/bins), bins, xRange[0], xRange[1]);
  TH1F* h_trk_z01 = new TH1F("trk_z01", Form(";Track z0 (cm); Tracks / %f cm", (xRange[1]-xRange[0])/bins), bins, xRange[0], xRange[1]);


  int nevt = tree->GetEntries();
  cout << "number of events = " << nevt << endl;


  // ----------------------------------------------------------------------------------------------------------------
  // event loop, for each event
  for (int i = 0; i < nevt; i++) {
    tree->GetEntry(i, 0);
    tree1->GetEntry(i, 0);

    // ----------------------------------------------------------------------------------------------------------------
    // loop over *all* found found L1 tracks
    for (int it = 0; it < (int)trk_pt->size(); it++) {

      // here you can place some selection criteria and fill histograms for the L1 tracks

      h_trk_pt->Fill(trk_pt->at(it)); // fill histogram "h_trk_pt" with the track pt for each track in the event
      h_trk_eta->Fill(trk_eta->at(it));
      h_trk_phi->Fill(trk_phi->at(it));
      h_trk_d0->Fill(trk_d0->at(it));
      h_trk_z0->Fill(trk_z0->at(it));

    } // end of track loop

    for (int it = 0; it < (int)trk_pt1->size(); it++) {
      h_trk_pt1->Fill(trk_pt1->at(it));
      h_trk_eta1->Fill(trk_eta1->at(it));
      h_trk_phi1->Fill(trk_phi1->at(it));
      h_trk_d01->Fill(trk_d01->at(it));
      h_trk_z01->Fill(trk_z01->at(it));
    }
  }

  TH1F *h1 = new TH1F("h1", "", 100,0,100.);
  TH1F *h2 = new TH1F("h2", "", 100,0,100.);

  if(what=="trk_eta"){
    h1 = ((TH1F *)gROOT->FindObject(what));
    h2 = ((TH1F *)gROOT->FindObject(what+"1"));
  }
    else if(what=="trk_phi"){
      h1 = ((TH1F *)gROOT->FindObject(what));
      h2 = ((TH1F *)gROOT->FindObject(what+"1"));
    }
    else if(what=="trk_d0"){
      h1 = ((TH1F *)gROOT->FindObject(what));
      h2 = ((TH1F *)gROOT->FindObject(what+"1"));
    }
    else if(what=="trk_z0"){
      h1 = ((TH1F *)gROOT->FindObject(what));
      h2 = ((TH1F *)gROOT->FindObject(what+"1"));
    }
    else if(what=="trk_pt"){
      h1 = ((TH1F *)gROOT->FindObject(what));
      h2 = ((TH1F *)gROOT->FindObject(what+"1"));
    }

  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  h1->SetLineColor(1);
  h1->SetMarkerColor(1);
  h1->SetMarkerStyle(8);
  h1->Draw(""); //l draws a line, p draws a marker
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->SetMarkerStyle(24);
  h2->Draw(",same");

  TLegend* l;
  if (what.Contains("trk")) {
    l = new TLegend(0.55,0.22,0.85,0.40);
    mySmallText(0.25,0.22,1, "Disp muon, PU=0");
  }
  else {
    l = new TLegend(0.2,0.72,0.5,0.9);
    mySmallText(0.6,0.82,1, "Disp muon, PU=0");
  }
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h1,"no rebin","ep");
  l->AddEntry(h2,"no rebin + no approx","ep");
  l->SetTextFont(42);
  l->Draw();

  gPad->SetGridy();

  c1->Draw();

}


void mySmallText(Double_t x,Double_t y,Color_t color,char *text) {
  Double_t tsize=0.044;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
