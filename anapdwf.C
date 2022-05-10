#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TChain.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

//===================Function filter===============================================================

void TV1D_denoise(vector<float>& input, vector<float>& output, const int width, const float lambda){
  if (width>0) {                /*to avoid invalid memory access to input[0]*/
    int k=0, k0=0;            /*k: current sample location, k0: beginning of current segment*/
    float umin=lambda, umax=-lambda;    /*u is the dual variable*/
    float vmin=input[0]-lambda, vmax=input[0]+lambda;    /*bounds for the segment's value*/
    int kplus=0, kminus=0;     /*last positions where umax=-lambda, umin=lambda, respectively*/
    const float twolambda=2.0*lambda;    /*auxiliary variable*/
    const float minlambda=-lambda;        /*auxiliary variable*/
    for (;;) {                /*simple loop, the exit test is inside*/
      //cout<<"k = "<<k<<endl;
      if (k>=width) return;
      while (k==width-1) {    /*we use the right boundary condition*/
        if (umin<0.0) {            /*vmin is too high -> negative jump necessary*/
          do{
            //cout<<"k0 = "<<k0<<" kminus = "<<kminus<<" vmin = "<<vmin<<endl;
            output[k0++]=vmin;
          }while (k0<=kminus);
          umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
        } else if (umax>0.0) {    /*vmax is too low -> positive jump necessary*/
          do{
            //cout<<"k0 = "<<k0<<" kplus = "<<kplus<<" vmax = "<<vmax<<endl;            
            output[k0++]=vmax;
          }while (k0<=kplus);
          umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
        } else {
          vmin+=umin/(k-k0+1);
          do{
            //cout<<"k0 = "<<k0<<" k = "<<k<<" vmin = "<<vmin<<endl;                        
            output[k0++]=vmin;
          }while(k0<=k);
          return;
        }
      }
      if ((umin+=input[k+1]-vmin)<minlambda) {        /*negative jump necessary*/
        do{
          //cout<<"negative jump, k0 = "<<k0<<" vmin = "<<vmin<<endl;
          output[k0++]=vmin;
        }while (k0<=kminus);
        vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
        umin=lambda; umax=minlambda;
      } else if ((umax+=input[k+1]-vmax)>lambda) {    /*positive jump necessary*/
        do{
          //cout<<"positive jump, k0 = "<<k0<<" vmax = "<<vmax<<endl;
          output[k0++]=vmax;
        }while (k0<=kplus);
        vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
        umin=lambda; umax=minlambda;
      } else {     /*no jump necessary, we continue*/
        k++;
        if (umin>=lambda) {        /*update of vmin*/
          vmin+=(umin-lambda)/((kminus=k)-k0+1);
          umin=lambda;
        }
        if (umax<=minlambda) {    /*update of vmax*/
          vmax+=(umax+lambda)/((kplus=k)-k0+1);
          umax=minlambda;
        }
      }
    }
  }
}


void anapdwf(){

  gStyle->SetOptStat(0);
  
  TFile f1("hist.root","recreate");
  TH1F *wfor = new TH1F("wfor","original waveform", 2000, 0, 2000);
  TH1F *wfma = new TH1F("wfma","moving average", 2000, 0, 2000);
  TH1F *wfden = new TH1F("wfden","denoised", 2000, 0, 2000);

  const int nch = 288;
  const int nwf = 1000;
  vector<vector<TH1F*>> wfnoise(nch);
  vector<vector<TH1F*>> wfsignal(nch);
  for (size_t i = 0; i<wfnoise.size(); ++i){
    wfnoise[i].resize(nwf);
    for (size_t j = 0; j<wfnoise[i].size(); ++j){
      wfnoise[i][j] = new TH1F(Form("wfnoise_%zu_%zu",i,j), Form("wfnoise_%zu_%zu",i,j), 2000, 0, 2000);
    }
    wfsignal[i].resize(nwf);
    for (size_t j = 0; j<wfsignal[i].size(); ++j){
      wfsignal[i][j] = new TH1F(Form("wfsignal_%zu_%zu",i,j), Form("wfsignal_%zu_%zu",i,j), 2000, 0, 2000);
    }
  }
  
  TH1F *wfsum[nch];
  for (int i = 0; i<nch; ++i){
    wfsum[i] = new TH1F(Form("wfsum%d",i),Form("ch %d",i), 1100, -1000, 10000);
  }
  
  map<int, int> signalcounter;
  map<int, int> noisecounter;
  
  TChain* pdtree = new TChain("pdtree");
  std::ifstream in;
  in.open("pd_6849.xroot");
  //in.open("pd_6852.xroot");
  char line[1024];
  
  while(1){
    in.getline(line,1024);
    if (!in.good()) break;
    pdtree->Add(Form("%s/pdwaveform/pdtree", line));
  }
  in.close();
  in.clear();
//  TFile f("pdwaveform.root");
//  TTree *pdtree = (TTree*)f.Get("pdwaveform/pdtree");
  int run, event, daqch;
  vector<short> *onda = 0;
  pdtree->SetBranchAddress("run", &run);
  pdtree->SetBranchAddress("event", &event);
  pdtree->SetBranchAddress("daqch", &daqch);
  pdtree->SetBranchAddress("onda", &onda);

  std::vector<float> waveform(2000);
  std::vector<float> ondama(2000);
  std::vector<float> ondaden(2000);
  std::vector<float> ondaor(2000);
  std::vector<float> input(2000);
  std::vector<float> output(2000);

  unsigned int mobileAVG = 4;
  int ngraphs = 10;
  int igraph = 0;
  TCanvas *can = new TCanvas("can", "can", 800, 400);
  can->Print("pdwf.ps[");

  int nentries = pdtree->GetEntries();
  for (int i = 0; i<nentries; ++i){
    //for (int i = 0; i<10; ++i){
    if (i%1000000==0) cout<<i<<"/"<<nentries<<endl;
    //if (i!=3125312) continue;
    pdtree->GetEntry(i);
    if (daqch!=143) continue;
    //cout<<onda->size()<<endl;
    wfor->Reset();
    wfma->Reset();
    wfden->Reset();
    wfor->SetTitle(Form("Run %d Event %d Ch %d;time tick;ADC",run,event,daqch));
    fill(waveform.begin(), waveform.end(), 0);
    fill(ondama.begin(), ondama.end(), 0);
    fill(ondaor.begin(), ondaor.end(), 0);
    fill(ondaden.begin(), ondaden.end(), 0);
    fill(input.begin(), input.end(), 0);
    fill(output.begin(), output.end(), 0);
    for (size_t j = 0; j<onda->size(); ++j){
      waveform[j] = (*onda)[j];
    }
    //========= Moving average  ============================================================
    float avg=0.;
    int c0=0, c1=0,  c2=0;
    for(size_t n=0; n<onda->size() ; ++n){
      if(n>=mobileAVG && n<onda->size()-mobileAVG){
        for(size_t k=n-mobileAVG; k<=n+mobileAVG; ++k){
          avg=avg+waveform[k];
          ++c0;
        }
        avg=avg/c0;
        ondama[n]=avg;
        avg = 0;
        c0 = 0;
      }
      else{
        if(n<mobileAVG){
          for(size_t k=0; k<=n+mobileAVG; ++k){
            avg=avg+waveform[k];
            c1=c1+1;
          }
          avg=avg/c1;
          ondama[n]=avg;
          avg = 0;
          c1 = 0;
        }
        else if(n>=onda->size()-mobileAVG){
          for(size_t k=n-mobileAVG; k<onda->size(); ++k){
            avg=avg+waveform[k];
            ++c2;
          }
          avg=avg/c2;
          ondama[n]=avg;
          avg = 0;
          c2 = 0;
        }
      }
    }

    //========= Denoising ============================================================
    
    for(size_t j=0; j<onda->size(); ++j){
      input[j]=ondama[j];
      output[j]=input[j];
    }
    TV1D_denoise(input,output,2000,10);
    for(size_t j=0; j<onda->size(); ++j){
      ondaden[j]=output[j];
    }

    //===========================BASELINE Histo========================================
    float base=0;
    TH1F *basehelp= new TH1F("basehelp","basehelp",2000, 1300,1800);
    for(size_t j=0; j<onda->size(); j++){
      if(j<1000){
        basehelp->Fill(ondaden[j]);
      }
    }
    int basebinmax = basehelp->GetMaximumBin();
    base = basehelp->GetXaxis()->GetBinCenter(basebinmax);
    basehelp->Delete();
    float sum = 0;
    for(size_t j = 0; j < onda->size(); ++j){
      ondaden[j] = ondaden[j]-base;
      ondaor[j] = waveform[j]-base;
      ondama[j] = ondama[j]-base;
      wfor->SetBinContent(j+1, ondaor[j]);
      wfma->SetBinContent(j+1, ondama[j]);
      wfden->SetBinContent(j+1, ondaden[j]);
//      f1.cd();
//      wfor->Write(Form("wfor_%d_%d_%d",run,event,daqch));
//      wfma->Write(Form("wfma_%d_%d_%d",run,event,daqch));
//      wfden->Write(Form("wfden_%d_%d_%d",run,event,daqch));
      if (j>=1000 && j<1500) sum += ondaden[j];
    }
    if (igraph<ngraphs){
      wfor->Draw();
      wfma->SetLineColor(2);
      wfma->Draw("same");
      wfden->SetLineColor(3);
      wfden->Draw("same");
      TLegend *leg = new TLegend(0.15,0.75,0.35,0.9);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(wfor,"Original","l");
      leg->AddEntry(wfma,"After moving averaging","l");
      leg->AddEntry(wfden,"After denoising","l");
      leg->Draw();
      can->Print("pdwf.ps");
      ++igraph;
    }
    wfsum[daqch]->Fill(sum);
    if (sum > -300 && sum < 300){
      if (noisecounter[daqch]<nwf){
        for (size_t j = 0; j<ondaor.size(); ++j){
          wfnoise[daqch][noisecounter[daqch]]->SetBinContent(j+1, ondaor[j]);
        }
        ++noisecounter[daqch];
      }
    }
    else if (sum > 600 && sum < 1000){
      if (signalcounter[daqch]<nwf){
        for (size_t j = 0; j<ondaor.size(); ++j){
          wfsignal[daqch][signalcounter[daqch]]->SetBinContent(j+1, ondaor[j]);
        }
        ++signalcounter[daqch];
      }
    }
    if (noisecounter[daqch]==nwf && signalcounter[daqch]==nwf) break;
  }
  can->Print("pdwf.ps]");
  f1.cd();
  for (int i = 0; i<nch; ++i){
    wfsum[i]->Write();
    for (int j = 0; j<nwf; ++j){
      wfsignal[i][j]->Write();
      wfnoise[i][j]->Write();
    }
  }
  f1.Close();
}
