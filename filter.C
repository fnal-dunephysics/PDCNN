{

  TFile f("fft.root");

  TProfile *freq_sig = (TProfile*)f.Get("freq_signal");
  TProfile *freq_noise = (TProfile*)f.Get("freq_noise");
  TH1D *hfreq_sig = freq_sig->ProjectionX("hfreq_sig");
  TH1D *hfreq_noise = freq_noise->ProjectionX("hfreq_noise");
  
  TCanvas *c1 = new TCanvas("c1","c1");
  freq_sig->SetTitle("");
  freq_sig->Draw();
  //gPad->SetLogy();
  //freq_sig->GetYaxis()->SetRangeUser(0,0.5);
  freq_noise->SetLineColor(2);
  freq_noise->SetMarkerColor(2);
  freq_noise->Draw("same");
  TLegend *leg1 = new TLegend(0.3,0.6,0.6,0.9);
  leg1->SetFillStyle(0);
  leg1->AddEntry(freq_sig, "1PE", "p");
  leg1->AddEntry(freq_noise, "Pedestal", "p");
  leg1->Draw();
  c1->Print("freq.png");
  c1->Print("freq.pdf");
  
  TH1D *filt = (TH1D*)hfreq_sig->Clone("filt");
  filt->Add(hfreq_noise, -1);
  TH1D *dem = (TH1D*)hfreq_sig->Clone("dem");
  dem->Add(hfreq_noise);
  filt->Divide(dem);
  TCanvas *c2 = new TCanvas("c2","c2");
  filt->SetTitle("");
  filt->GetYaxis()->SetTitle("(1PE-Pedestal)/(1PE+Pedestal)");
  filt->Draw();

  TF1 *fun = new TF1("fun","[0]*exp(-0.5*pow(x/[1],2))", 0, 10000);
  fun->SetParameters(0.7, 2);
  //fun->Draw("same");
  filt->Fit("fun","R","",0.,40);
  cout<<fun->GetParameter(0)<<" "<<fun->GetParameter(1)<<endl;
  c2->Print("filter.pdf");
  c2->Print("filter.png");
  
  
  // Convert filter to time domain
  int ntbin = 2000;
  double samplerate = 6.67; // [ns]
  double binwidth = 1./(ntbin*samplerate*1e-3); // [MHz]
  double *real_full = new double[ntbin];
  double *img_full = new double[ntbin];

  fun->SetParameter(0, 1./ntbin);

  for (int j = 0; j<ntbin; ++j){
    double frequency = (j+0.5)*binwidth;
    double filt = fun->Eval(frequency);
    real_full[j] = filt;
    img_full[j] = 0;
//    if (j==0){
//      real_full[j] = 0;
//      img_full[j] = 0;
//    }
  }
  
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &ntbin, "C2R M K");
  fft_back->SetPointsComplex(&real_full[0], &img_full[0]);
  fft_back->Transform();
  TH1 *h1_fft_back = new TH1F("h1_fft_back","filter",ntbin, 0, ntbin);
  h1_fft_back = TH1::TransformHisto(fft_back,h1_fft_back,"Re");
  TCanvas *c3 = new TCanvas("c3","c3");
  h1_fft_back->Draw("hist");
  cout<<h1_fft_back->Integral()<<endl;
}
