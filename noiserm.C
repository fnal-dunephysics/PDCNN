{

  TFile f("hist_6852.root");

  const int n = 100;
  
  TH1D *waveform[n];

  int ntbin = 2000;
  double samplerate = 6.67; // [ns]
  double binwidth = 1./(ntbin*samplerate*1e-3); // [MHz]

  double *real_full = new double[ntbin];
  double *img_full = new double[ntbin];

  TF1 *filtfun = new TF1("filtfun","[0]*exp(-0.5*pow(x/[1],2))", 0, 250);
  filtfun->SetParameters(1./ntbin, 2.01642);

  TCanvas *c1 = new TCanvas("c1","c1",500,800);
  c1->Divide(1,2);
  c1->Print("noiserm.ps[");
  for (int i = 0; i<n; ++i){

    waveform[i] = (TH1D*)f.Get(Form("wfsignal_143_%d",i));
    if (!waveform[i]){
      cout<<"Cannot get histogram "<<Form("wfsignal_143_%d",i)<<endl;
      exit(1);
    }

    c1->cd(1);
    waveform[i]->Draw("hist");
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
    c1->cd(2);
    h1_fft_back->DrawCopy("hist");
    c1->Print("noiserm.ps");
    delete re_fft;
    delete im_fft;
    delete h1_fft_back;
  }
  c1->Print("noiserm.ps]");

}
