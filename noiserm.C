double fitf(double *x, double *par){

  if (x[0]<par[3]) return 0;
  else if (x[0]<par[4]) return par[0]*(1-exp(-(x[0]-par[3])/par[1]));
  else return par[0]*(1-exp(-(par[4]-par[3])/par[1]))*exp(-(x[0]-par[4])/par[2]);
}

void noiserm(){

  TFile f("hist_6852.root");

  const int n = 1000;
  
  TH1D *waveform[n];
  TH1D *wfsum = new TH1D("wfsum","summed ADC", 1100, -100, 10000);
  TH1D *wfsumnf = new TH1D("wfsumnf","summed ADC, noise filtered", 1100, -100, 10000);

  int ntbin = 2000;
  double samplerate = 6.67; // [ns]
  double binwidth = 1./(ntbin*samplerate*1e-3); // [MHz]

  double *real_full = new double[ntbin];
  double *img_full = new double[ntbin];
  double *real_spe = new double[ntbin];
  double *img_spe = new double[ntbin];

  TF1 *filtfun = new TF1("filtfun","[0]*exp(-0.5*pow(x/[1],2))", 0, 250);
  filtfun->SetParameters(1./ntbin, 2.01642);

  TF1 *func = new TF1("fitf", "fitf", 0, 10000, 5);
  func->SetParameter(0, 10.2314);
  func->SetParameter(1, 8.10114);
  func->SetParameter(2, 66.4333);
  func->SetParameter(3, 72.9817);
  func->SetParameter(4, 94.6331);
  TH1D *spe = new TH1D("spe", "spe", ntbin, 0, ntbin);
  for (int i = 1; i<=spe->GetNbinsX(); ++i){
    spe->SetBinContent(i, func->Eval(spe->GetBinCenter(i)));
  }
  TH1 *re_spe = 0;
  re_spe = spe->FFT(re_spe, "RE");
  
  TH1 *im_spe = 0;
  im_spe = spe->FFT(im_spe, "IM");
  
  TCanvas *c1 = new TCanvas("c1","c1",500,900);
  c1->Divide(1,3);
  c1->Print("noiserm.ps[");
  for (int i = 0; i<n; ++i){
    double sum = 0;
    waveform[i] = (TH1D*)f.Get(Form("wfsignal_143_%d",i));
    for (int j = 1000; j<1500; ++j) sum+= waveform[i]->GetBinContent(j+1);
    wfsum->Fill(sum);
    if (!waveform[i]){
      cout<<"Cannot get histogram "<<Form("wfsignal_143_%d",i)<<endl;
      exit(1);
    }

    c1->cd(1);
    waveform[i]->Draw("hist");

    //Filter noise
    c1->cd(2);    
    TH1 *re_fft = 0;
    re_fft = waveform[i]->FFT(re_fft, "RE");

    TH1 *im_fft = 0;
    im_fft = waveform[i]->FFT(im_fft, "IM");
    
    //cout<<re_fft->GetNbinsX()<<" "<<im_fft->GetNbinsX()<<endl;

    for (int j = 0; j<ntbin; ++j){
      real_full[j] = re_fft->GetBinContent(j+1);
      img_full[j] = im_fft->GetBinContent(j+1);
      // If the following is commented out, the original waveform will be retrieved. 
      double frequency = (j+0.5)*binwidth;
      double filt = filtfun->Eval(frequency);
      real_full[j] *= filt;
      img_full[j] *= filt;
//      if (j==0){
//        real_full[j] = 0;
//        img_full[j] = 0;
//      }
    }
    TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &ntbin, "C2R M K");
    fft_back->SetPointsComplex(&real_full[0], &img_full[0]);
    fft_back->Transform();
    TH1 *h1_fft_back = (TH1*)waveform[i]->Clone("h1_fft_back");
    h1_fft_back->Reset();
    h1_fft_back = TH1::TransformHisto(fft_back,h1_fft_back,"Re");
    h1_fft_back->DrawCopy("hist");
    sum = 0;
    for (int j = 1000; j<1500; ++j) sum+= h1_fft_back->GetBinContent(j+1);
    wfsumnf->Fill(sum);

    c1->cd(3);
    for (int j = 0; j<ntbin; ++j){
      real_full[j] = re_fft->GetBinContent(j+1);
      img_full[j] = im_fft->GetBinContent(j+1);
      real_spe[j] = re_spe->GetBinContent(j+1);
      img_spe[j] = im_spe->GetBinContent(j+1);
      TComplex cn1(real_full[j], img_full[j]);
      TComplex cn2(real_spe[j], img_spe[j]);
      cn1/=cn2;
      real_full[j] = cn1.Re();
      img_full[j] = cn1.Im();
      double frequency = (j+0.5)*binwidth;
      double filt = filtfun->Eval(frequency);
      real_full[j] *= filt;
      img_full[j] *= filt;
    }
    TVirtualFFT *fft_back1 = TVirtualFFT::FFT(1, &ntbin, "C2R M K");
    fft_back1->SetPointsComplex(&real_full[0], &img_full[0]);
    fft_back1->Transform();
    TH1 *h1_fft_back1 = (TH1*)waveform[i]->Clone("h1_fft_back1");
    h1_fft_back1->Reset();
    h1_fft_back1 = TH1::TransformHisto(fft_back1,h1_fft_back1,"Re");
    h1_fft_back1->DrawCopy("hist");
    
    c1->Print("noiserm.ps");
    delete re_fft;
    delete im_fft;
    delete h1_fft_back;
  }
  c1->Print("noiserm.ps]");

//  c1->cd(1);
//  wfsum->Draw();
//  c1->cd(2);
//  wfsumnf->Draw();
  
}
