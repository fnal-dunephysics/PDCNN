double fitf(double *x, double *par){

  if (x[0]<par[3]) return 0;
  else if (x[0]<par[4]) return par[0]*(1-exp(-(x[0]-par[3])/par[1]));
  else return par[0]*(1-exp(-(par[4]-par[3])/par[1]))*exp(-(x[0]-par[4])/par[2]);
}

void SPE(){

  TFile f("fft.root");

  TH1F *avgsignal = (TH1F*)f.Get("avgsignal");

  TF1 *func = new TF1("fitf", "fitf", 0, 500, 5);
  func->SetParameter(0, 9);
  func->SetParameter(1, 8);
  func->SetParameter(2, 70);
  func->SetParameter(3, 72);
  func->SetParameter(4, 90);
  //func->FixParameter(4, 90);
  avgsignal->Fit(func);

  TCanvas *c1 = new TCanvas("c1","c1");
  avgsignal->DrawCopy();

  for (int i = 0; i<5; ++i){
    cout<<"par["<<i<<"]="<<func->GetParameter(i)<<endl;
  } 

  TCanvas *c2 = new TCanvas("c2","c2");
  
}

  
