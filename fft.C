void fft(){

  TFile f("hist_6852.root");

  //const int nch = 288;
  const int nwf = 1000;

  TH1F *wfnoise[nwf];
  TH1F *wfsignal[nwf];

  for (int i = 0; i<nwf; ++i){
    wfnoise[i] = (TH1F*)f.Get(Form("wfnoise_143_%d",i));
    wfsignal[i] = (TH1F*)f.Get(Form("wfsignal_143_%d",i));
    if (!wfnoise[i] || !wfsignal[i]){
      cout<<"histogram empty"<<endl;
      exit(1);
    }
  }

  int ntbin = 500;
  double samplerate = 6.67; //ns
  double binwidth = 1./(ntbin*samplerate*1e-3); // [MHz]
  
  TH1D *htemp = new TH1D("htemp","htemp",ntbin,0,ntbin);

  TFile fout("fft.root","recreate");
  
  TProfile *freq_signal = new TProfile("freq_signal","Signal;Frequency (MHz);Amplitude",ntbin,0,ntbin*binwidth);
  TProfile *freq_noise = new TProfile("freq_noise","Noise;Frequency (MHz);Amplitude",ntbin,0,ntbin*binwidth);

  TH1F *avgnoise = new TH1F("avgnoise","Average noise spectrum", ntbin, 0, ntbin);
  TH1F *avgsignal = new TH1F("avgsignal","Average signal spectrum", ntbin, 0, ntbin);

  for (int i = 0; i<nwf; ++i){
    htemp->Reset();
    for (int j = 1; j<=ntbin; ++j){
      htemp->SetBinContent(j, wfnoise[i]->GetBinContent(j+1000));
      avgnoise->SetBinContent(j, avgnoise->GetBinContent(j) + wfnoise[i]->GetBinContent(j+1000));
    }
    TH1 *h1_fft = 0;
    h1_fft = htemp->FFT(h1_fft, "MAG");
    for (int j = 0; j<h1_fft->GetNbinsX(); ++j){
      double frequency = (j+0.5)*binwidth;
      double amp = h1_fft->GetBinContent(j+1);
      freq_noise->Fill(frequency, amp);
    }
    if (h1_fft) delete h1_fft;
  }

  for (int i = 0; i<nwf; ++i){
    htemp->Reset();
    for (int j = 1; j<=ntbin; ++j){
      htemp->SetBinContent(j, wfsignal[i]->GetBinContent(j+1000));
      avgsignal->SetBinContent(j, avgsignal->GetBinContent(j) + wfsignal[i]->GetBinContent(j+1000));
    }
    TH1 *h1_fft = 0;
    h1_fft = htemp->FFT(h1_fft, "MAG");
    for (int j = 0; j<h1_fft->GetNbinsX(); ++j){
      double frequency = (j+0.5)*binwidth;
      double amp = h1_fft->GetBinContent(j+1);
      freq_signal->Fill(frequency, amp);
    }
    if (h1_fft) delete h1_fft;
  }

  avgnoise->Scale(1./nwf);
  avgsignal->Scale(1./nwf);
  
  fout.Write();
  fout.Close();
  
}
