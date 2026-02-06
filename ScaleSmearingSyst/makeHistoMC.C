Double_t DoubleSidedCB2(double x, double mu, double width, double a1, double p1, double a2, double p2)
{
  double u   = (x-mu)/width;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(1);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}


double DoubleSidedCB(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]));
}

void FitHisto(TH1D* _hist, TString name, string period)
{
  TF1 *dcb = new TF1("dcb","DoubleSidedCB",80, 105,7);
  dcb->SetParLimits(1,80.,105.);
  dcb->SetParLimits(2,0.,5.);
  dcb->SetParLimits(3,0.,2.);
  dcb->SetParLimits(4,0.,10.);
  dcb->SetParLimits(5,0.,2.);
  dcb->SetParLimits(6,0.,10.);

  dcb->SetParameters(_hist->Integral()/5.,91.,2.5,1.,1.,1.,1.);

  dcb->SetParName(0,"C");
  dcb->SetParName(1,"Mean");
  dcb->SetParName(2,"Sigma");
  dcb->SetParName(3,"#alpha_{1}");
  dcb->SetParName(4,"n_{1}");
  dcb->SetParName(5,"#alpha_{2}");
  dcb->SetParName(6,"n_{2}");

  _hist->Fit("dcb");
  Double_t mass = dcb->GetParameter(1);
  Double_t width = dcb->GetParameter(2);
  Double_t mass_err = dcb->GetParError(1);
  Double_t width_err = dcb->GetParError(2);

  TCanvas * canvas = new TCanvas("canvas","canvas",800,800);
  _hist->Draw();
  dcb->Draw("SAME");
  canvas->SaveAs((TString("plots_") + period.c_str() + "/" + name + ".png").Data());
  canvas->SaveAs((TString("plots_") + period.c_str() + "/" + name + ".pdf").Data());


  std::ofstream outfile;
  outfile.open("FitResults/"+(TString)name+".txt", std::ios_base::trunc);
  outfile << Form("%.3f \n", mass);
  outfile << Form("%.3f \n", width);
  outfile << Form("%.3f \n", mass_err);
  outfile << Form("%.4f \n", width_err);

}

float makemass(vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi) {
    TLorentzVector lep_1;
        TLorentzVector lep_2;
        lep_1.SetPtEtaPhiM(LepPt[0], LepEta[0], LepPhi[0], 0.00001);
        lep_2.SetPtEtaPhiM(LepPt[1], LepEta[1], LepPhi[1], 0.00001);
        return (lep_1+lep_2).M();
}

void makeHistoMC(){
  gROOT->SetBatch(kTRUE);

  string fpath;
  string proc = "DY_2023postBPix";

   double lumi = 0.0;

    // Luminosity map
    std::map<std::string, double> lumi_map = {
        {"2022",  7.98},
        {"2022EE",  26.67},
        {"2023preBPix", 18.06},
        {"2023postBPix", 9.69},
        {"2024",108.82}
    };

    // Extract dataset period from proc string
    string period = proc.substr(proc.find("_") + 1);

    // Lookup luminosity
    if (lumi_map.count(period))
        lumi = lumi_map[period];
    else {
        cerr << "ERROR: Unknown lumi period in proc: " << proc << endl;
        return;
    }

    cout << "Using lumi = " << lumi << " fb^-1" << endl;

  fpath = "/eos/user/m/mmanoni/ScaleSmearingSyst_skimmed_rootfiles/" 
        + period 
        + "DY_ZLCand_ScaleSmearingSyst.root";
  
  //"/eos/user/m/mmanoni/ScaleSmearingSyst_skimmed_rootfiles/2022EEDY_ZLCand_ScaleSmearingSyst.root";
  //fpath = "/eos/user/m/mmanoni/ScaleSmearingSyst_skimmed_rootfiles/2022EEDY_ZLCand_ScaleSmearingSyst.root";


  TFile * f          = new TFile((TString)fpath,"READ");
  TH1F  * t1          = (TH1F*)f->Get("CRZLTree/Counters");
  double gen_sumWeights = t1->GetBinContent(40);

  ROOT::RDataFrame rdf("CRZLTree/candTree", fpath);

  auto ext_df = (rdf.Define("ZMass_ScaleUp", makemass, {"LepPt_ScaleUp", "LepEta", "LepPhi"})
                    .Define("ZMass_ScaleDn", makemass, {"LepPt_ScaleDn", "LepEta", "LepPhi"})
                    .Define("ZMass_SmearUp", makemass, {"LepPt_SmearUp", "LepEta", "LepPhi"})
                    .Define("ZMass_SmearDn", makemass, {"LepPt_SmearDn", "LepEta", "LepPhi"})
                );

  vector<double> eta_cut {0, -2.6, -1.95, -1.3, -0.65, 0, 0.65, 1.3, 1.95, 2.6};
  double low_eta = -10;
  double high_eta = 10;
  string s_eta_cut = "";

  TFile * histoMC = new TFile("FitResults/histoMC_"+(TString)proc+".root", "RECREATE");
  histoMC->cd();

  for(int i = 2; i < eta_cut.size() ; i++){
      low_eta = eta_cut.at(i-1);
      high_eta = eta_cut.at(i);
      std::cout << low_eta << " < eta < " << high_eta << std::endl;

      s_eta_cut = Form("etaBin%i", i-1);

      std::cout << s_eta_cut << std::endl;
      
      auto hist_ele = rdf.Filter(Form("Z1Mass>=80 && Z1Mass<=105 && Z1Flav==-121 && LepEta[0] < %f && LepEta[0] > %f && LepEta[1] < %f && LepEta[1] > %f", high_eta, low_eta, high_eta, low_eta)).Define("weight", Form("1000 * overallEventWeight * %f / %f", lumi, gen_sumWeights)).Histo1D({"Z1Mass", "Z1Mass", 20, 80., 105.},"Z1Mass", "weight");

      auto hist_ele_smearUp = ext_df.Filter(Form("ZMass_SmearUp>=80 && ZMass_SmearUp<=105 && Z1Flav==-121 && LepEta[0] < %f && LepEta[0] > %f && LepEta[1] < %f && LepEta[1] > %f", high_eta, low_eta, high_eta, low_eta)).Define("weight", Form("1000 * overallEventWeight * %f / %f", lumi, gen_sumWeights)).Histo1D({"ZMass_SmearUp", "ZMass_SmearUp", 20, 80., 105.},"ZMass_SmearUp", "weight");

      auto hist_ele_smearDn = ext_df.Filter(Form("ZMass_SmearDn>=80 && ZMass_SmearDn<=105 && Z1Flav==-121 && LepEta[0] < %f && LepEta[0] > %f && LepEta[1] < %f && LepEta[1] > %f", high_eta, low_eta, high_eta, low_eta)).Define("weight", Form("1000 * overallEventWeight * %f / %f", lumi, gen_sumWeights)).Histo1D({"ZMass_SmearDn", "ZMass_SmearDn", 20, 80., 105.},"ZMass_SmearDn", "weight");

      auto hist_ele_scaleUp = ext_df.Filter(Form("ZMass_ScaleUp>=80 && ZMass_ScaleUp<=105 && Z1Flav==-121 && LepEta[0] < %f && LepEta[0] > %f && LepEta[1] < %f && LepEta[1] > %f", high_eta, low_eta, high_eta, low_eta)).Define("weight", Form("1000 * overallEventWeight * %f / %f", lumi, gen_sumWeights)).Histo1D({"ZMass_ScaleUp", "ZMass_ScaleUp", 20, 80., 105.},"ZMass_ScaleUp", "weight");

      auto hist_ele_scaleDn = ext_df.Filter(Form("ZMass_ScaleDn>=80 && ZMass_ScaleDn<=105 && Z1Flav==-121 && LepEta[0] < %f && LepEta[0] > %f && LepEta[1] < %f && LepEta[1] > %f", high_eta, low_eta, high_eta, low_eta)).Define("weight", Form("1000 * overallEventWeight * %f / %f", lumi, gen_sumWeights)).Histo1D({"ZMass_ScaleDn", "ZMass_ScaleDn", 20, 80., 105.},"ZMass_ScaleDn", "weight");

      hist_ele->SetTitle(Form("Nominal (%s)", s_eta_cut.c_str()));
      hist_ele_smearUp->SetTitle(Form("Smear Up (%s)", s_eta_cut.c_str()));
      hist_ele_smearDn->SetTitle(Form("Smear Dn (%s)", s_eta_cut.c_str()));
      hist_ele_scaleUp->SetTitle(Form("Scale Up (%s)", s_eta_cut.c_str()));
      hist_ele_scaleDn->SetTitle(Form("Scale Dn (%s)", s_eta_cut.c_str()));

      hist_ele->Write();
      hist_ele_smearUp->Write();
      hist_ele_smearDn->Write();
      hist_ele_scaleUp->Write();
      hist_ele_scaleDn->Write();

      FitHisto(hist_ele.GetPtr(), proc+"_ele_"+s_eta_cut.c_str(), period);
      FitHisto(hist_ele_smearUp.GetPtr(), proc+"_ele_smearUp_"+s_eta_cut.c_str(), period);
      FitHisto(hist_ele_smearDn.GetPtr(), proc+"_ele_smearDn_"+s_eta_cut.c_str(), period);
      FitHisto(hist_ele_scaleUp.GetPtr(), proc+"_ele_scaleUp_"+s_eta_cut.c_str(), period);
      FitHisto(hist_ele_scaleDn.GetPtr(), proc+"_ele_scaleDn_"+s_eta_cut.c_str(), period);
  }
  histoMC->Close();

}
