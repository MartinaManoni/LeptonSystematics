# Study of Lepton Scale Factor (SF) uncertainties in H→ZZ→4l analysis
# This script:
#  - Loads ROOT ntuples with H→ZZ→4l MC samples
#  - Computes yields for different final states (4e, 4μ, 2e2μ)
#  - Propagates lepton SF uncertainties (trigger/reco/selection)
#  - Separates correlated vs uncorrelated effects across leptons
#  - Outputs relative yield variations (up/down) as percentages
################################################################################

import math
import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *
import numpy as np


import uproot
import matplotlib.pyplot as plt
import numpy as np
import os


gRandom.SetSeed(101)

################################################################################
# ANALYSIS CONFIGURATION
################################################################################

year = "2022"

if (year=="2022"):
   lumi= 7.98
   folder = '/eos/user/m/mmanoni/SKIMMED_syst_final/'
   file_name = '/ZZ4lAnalysis_SKIMMED_idx_new.root'

   List = ['2022']

elif (year == "2022EE"):
   lumi= 26.671

   folder = '/eos/user/m/mmanoni/SKIMMED_syst_final/'
   file_name = '/ZZ4lAnalysis_SKIMMED_idx_new.root'

   List = ['2022EE']

elif (year == "2023preBPix"):
   lumi= 18.062

   folder = '/eos/user/m/mmanoni/SKIMMED_syst_final/'
   file_name = '/ZZ4lAnalysis_SKIMMED_idx_new.root'
   List = ['2023preBPix']

elif (year == "2023postBPix"):
   lumi= 9.693

   folder = '/eos/user/m/mmanoni/SKIMMED_syst_final/'
   file_name = '/ZZ4lAnalysis_SKIMMED_idx_new.root'
   List = ['2023postBPix']

elif (year == "2024"):
   lumi= 108.822

   folder = '/eos/user/m/mmanoni/SKIMMED_syst_final/'
   file_name = '/final_2024.root'
   List = ['2024']

else:
   raise ValueError(f"Unknown year: {year}")

# Flags to control uncertainty treatment
correlated_leptons = True     # treat lepton SF uncertainties as correlated across all 4 leptons
uncorrelated_leptons = True   # treat lepton SF uncertainties as independent (quadrature sum)
corr_factor = 0               # correlation coefficient (ρ = 0 means no correlation)

print("List of samples", List)


################################################################################
# HELPER FUNCTION: combine uncertainties across leptons with correlation ρ
################################################################################


def sigma_event(rho, SF1, SF2, SF3, SF4, sigma1, sigma2, sigma3, sigma4):
   #print("rho, SF1, SF2, SF3, SF4, sigma1, sigma2, sigma3, sigma4")
   #print(rho, SF1, SF2, SF3, SF4, sigma1, sigma2, sigma3, sigma4)
   """
   Computes total fractional uncertainty for 4 leptons
   given individual SF errors σ_i and correlation ρ.
   Formula: Var_total = sum_i (σ_i/SF_i)^2 + 2ρ * sum_{i<j}(σ_iσ_j / SF_iSF_j)
   """
   rez = (sigma1/SF1)**2 + (sigma2/SF2)**2 + (sigma3/SF3)**2 + (sigma4/SF4)**2 + 2*rho*( sigma1*sigma2/SF1/SF2 + sigma1*sigma3/SF1/SF3 + sigma1*sigma4/SF1/SF4 + sigma2*sigma3/SF2/SF3 + sigma2*sigma4/SF2/SF4 + sigma3*sigma4/SF3/SF4)
   if rez < 0.000001:
      return 0
   else:
      return rez

################################################################################
# LOOP OVER SAMPLES
################################################################################

for Type in List:
   # Open the root file and get tree
   mc=TFile.Open(folder+Type+file_name)
   tree = mc.Get("ZZTree/candTree")
   counters = mc.Get("Counters")
   print(type(counters))
   NGen = counters.GetBinContent(40)
   print("NGen", NGen)
   
   # Set mass range and resolution
   lo=105
   hi=140
   nbins=(hi-lo)*1 # 1 GeV resolution


   # Calculate nominal value of yield
   h_nom = TH1F("h_nom", "h_nom", nbins, lo, hi)

   # Compute nominal event yields for each channel (4e, 4μ, 2e2μ)
   # Each Draw applies a selection and fills the histogram
   tree.Draw("ZZMass >> h_nom" , "(abs(LepLepId[0]) == 11 && abs(LepLepId[3]) == 11)*overallEventWeight*dataMCWeight*1000*"+str(lumi)+"/"+str(NGen))
   nom_yield_4e = h_nom.Integral()
   print("nom_yield_4e", nom_yield_4e)
   
   tree.Draw("ZZMass >> h_nom" , "(abs(LepLepId[0]) == 13 && abs(LepLepId[3]) == 13)*overallEventWeight*dataMCWeight*1000*"+str(lumi)+"/"+str(NGen))
   nom_yield_4mu = h_nom.Integral()
   print("nom_yield_4mu", nom_yield_4mu)
   
   tree.Draw("ZZMass >> h_nom" , "(abs(abs(LepLepId[0]) - (LepLepId[3])) == 2)*overallEventWeight*dataMCWeight*1000*"+str(lumi)+"/"+str(NGen))
   nom_yield_2e2mu = h_nom.Integral()
   print("nom_yield_2e2mu", nom_yield_2e2mu)

   #############################################################################
   # INITIALIZATION: counters, arrays, labels
   #############################################################################
   br_data = 0     # event counter
   br_fs4e = 0     # number of 4e events
   br_fs4mu = 0    # number of 4μ events
   br_fs2e2mu = 0  # number of 2e2μ events
   
   # Store yield variations for trigger and reco/selection (2 types of variations)
   variation_name = ["TRIGGER", "RECO_SEL"]

   # Initialize all arrays (each with 2 entries → one per variation type)
   yield_4e_up = [0., 0.]
   yield_4mu_up = [0., 0.]
   yield_2e2mu_e_up = [0., 0.]
   yield_2e2mu_mu_up = [0., 0.]
   yield_4e_dn = [0., 0.]
   yield_4mu_dn = [0., 0.]
   yield_2e2mu_e_dn = [0., 0.]
   yield_2e2mu_mu_dn = [0., 0.]

   # Same for uncorrelated case
   uncor_4e_up = [0., 0.]
   uncor_4e_dn = [0., 0.]
   uncor_4mu_up = [0., 0.]
   uncor_4mu_dn = [0., 0.]
   uncor_2e2mu_e_up = [0., 0.]
   uncor_2e2mu_e_dn = [0., 0.]
   uncor_2e2mu_mu_up = [0., 0.]
   uncor_2e2mu_mu_dn = [0., 0.]

   
   #############################################################################
   # LOOP OVER EVENTS
   #############################################################################
   
   SF_sigma_4e_0 = []
   SF_sigma_4e_1 = []
   SF_sigma_4e_2 = []
   SF_sigma_4e_3 = []

   SF_sigma_4mu_0 = []
   SF_sigma_4mu_1 = []
   SF_sigma_4mu_2 = []
   SF_sigma_4mu_3 = []

   plot_Pt = []
   plot_Eta = []
   plot_sf_ratio2 = []

   plot_Pt_mu = []
   plot_Eta_mu = []
   plot_sf_ratio2_mu = []

   nEvents = tree.GetEntries()

   for event in range(nEvents):
      tree.GetEntry(event)

      # RESET EVERY EVENT — no exceptions
      SF_lep       = [0., 0., 0., 0.]
      err_lep_up   = [0., 0., 0., 0.]
      err_lep_dn   = [0., 0., 0., 0.]

      SF_lep_trig       = [1., 1., 1., 1.]
      err_lep_trig_up   = [0., 0., 0., 0.]
      err_lep_trig_dn   = [0., 0., 0., 0.]

      br_data+=1
      
      mass4l = tree.ZZMass
      if( mass4l < 105. or mass4l >140.): continue # Skip events that are not in the mass window

      # Identify final state based on lepton IDs
      idL1 = abs(tree.LepLepId[0])
      idL3 = abs(tree.LepLepId[3])

      # Increment counters for each final state
      if (idL1==11 and idL3==11):
         br_fs4e +=1
      elif (idL1==13 and idL3==13):
         br_fs4mu += 1
      elif (abs(idL1-idL3)==2):
         br_fs2e2mu += 1
      
      # Nominal value of total SF, product of 4 lepton nominal SF
      SF_tot_nom = tree.dataMCWeight

      # Calculate nominal weigh using central value of SF
      weight_nom = tree.overallEventWeight * 1000 * SF_tot_nom * lumi / NGen #xsec is already inside overallEventWeight --- what about L1prefiringWeight?
 
      for i in range (0,4):
         # Hard-coded trigger SF and unc (Since Run2)
         err_lep_trig_dn = [0.01,0,0,0]
         SF_lep_trig[i] = 1.
 
         if(idL1==11 and idL3==11 and tree.LepPt[3] < 12):
               err_lep_trig_up = [0.02,0,0,0]
               err_lep_trig_dn = [0.11,0,0,0]
         if(idL1==11 and idL3==11 and tree.LepPt[3] >= 12):
               err_lep_trig_up = [0.005,0,0,0]
               err_lep_trig_dn = [0.01,0,0,0]
         if(idL1==13 and idL3==13 and tree.LepPt[3] < 7):
               err_lep_trig_up = [0.001,0,0,0]
               err_lep_trig_dn = [0.032,0,0,0]
         if(idL1==13 and idL3==13 and tree.LepPt[3] >= 7 and tree.LepPt[3] < 12):
               err_lep_trig_up = [0.001,0,0,0]
               err_lep_trig_dn = [0.015,0,0,0]
         if(idL1==13 and idL3==13 and tree.LepPt[3] >= 12):
               err_lep_trig_up = [0.001,0,0,0]
               err_lep_trig_dn = [0.015,0,0,0]
                  
         if(abs(idL1-idL3)==2 and tree.LepPt[3] < 7):
               err_lep_trig_up = [0.005,0,0,0]
               err_lep_trig_dn = [0.08,0,0,0]
         if(abs(idL1-idL3)==2 and tree.LepPt[3] >= 7 and tree.LepPt[3] < 12):
               err_lep_trig_up = [0.005,0,0,0]
               err_lep_trig_dn = [0.032,0,0,0]
         if(abs(idL1-idL3)==2 and tree.LepPt[3] >= 12):
               err_lep_trig_up = [0.001,0,0,0]
               err_lep_trig_dn = [0.01,0,0,0]
      
         if abs(tree.LepLepId[i]) == 11:

            SF_lep[i] = tree.LepSF[i]
            err_lep_up[i] = tree.LepSFUnc[i]
            err_lep_dn[i] = tree.LepSFUnc[i]

         elif abs(tree.LepLepId[i]) == 13:
            SF_lep[i] = tree.LepSF[i]
            err_lep_up[i] = tree.LepSFUnc[i]
            err_lep_dn[i] = tree.LepSFUnc[i]
      
      for i in range(4):

         # Compute selections using existing variables
         eta_ok = abs(tree.LepEta[i]) <= 2.5
         pt_ok = (
            (abs(tree.LepLepId[i]) == 13 and tree.LepPt[i] >= 5  and tree.LepPt[i] <= 120) or
            (abs(tree.LepLepId[i]) == 11 and tree.LepPt[i] >= 7  and tree.LepPt[i] <= 500)
         )

      for k in range (0,2):
         # Calculate each variation independently
         if(k==0):
            TRIG = 1
            RECO_SEL = 0
         elif(k==1):
            TRIG = 0
            RECO_SEL = 1


         SF_var_up = 1.
         SF_var_dn = 1.
         SF_var_e_up = 1.
         SF_var_e_dn = 1.
         SF_var_mu_up = 1.
         SF_var_mu_dn = 1.
         
         if (idL1==11 and idL3==11):
            # Vary SF of each lepton up and down
            for i in range (0,4):
               if (correlated_leptons ):
                  SF_var_up *= (SF_lep_trig[i] + TRIG*err_lep_trig_up[i]) * (SF_lep[i] + RECO_SEL*err_lep_up[i])
                  SF_var_dn *= (SF_lep_trig[i] - TRIG*err_lep_trig_dn[i]) * (SF_lep[i] - RECO_SEL*err_lep_dn[i])
           
            if (uncorrelated_leptons):
               #print("here")
               uncor_4e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_up[0],err_lep_up[1],err_lep_up[2],err_lep_up[3]))
               SF_sigma_4e_0.append((err_lep_up[0] / SF_lep[0])**2)
               SF_sigma_4e_1.append((err_lep_up[1] / SF_lep[1])**2)
               SF_sigma_4e_2.append((err_lep_up[2] / SF_lep[2])**2)
               SF_sigma_4e_3.append((err_lep_up[3] / SF_lep[3])**2)

               # Loop over the 4 leptons to check SF_sigma_4e > 1
               for i in range(4):
                  sf_ratio2 = (err_lep_up[i] / SF_lep[i])**2
                  plot_Pt.append(tree.LepPt[i])
                  plot_Eta.append(tree.LepEta[i])
                  plot_sf_ratio2.append(sf_ratio2)
                  

               uncor_4e_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],err_lep_trig_dn[2],err_lep_trig_dn[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_dn[0],err_lep_dn[1],err_lep_dn[2],err_lep_dn[3])) 
               
            yield_4e_up[k] += weight_nom/SF_tot_nom * SF_var_up
            yield_4e_dn[k] += weight_nom/SF_tot_nom * SF_var_dn
      

         elif (idL1==13 and idL3==13):
            for i in range (0,4):
               if ( correlated_leptons ):
                  SF_var_up *= (SF_lep_trig[i] + TRIG*err_lep_trig_up[i]) * (SF_lep[i] + RECO_SEL*err_lep_up[i])
                  SF_var_dn *= (SF_lep_trig[i] - TRIG*err_lep_trig_dn[i]) * (SF_lep[i] - RECO_SEL*err_lep_dn[i])
            if (uncorrelated_leptons):

               uncor_4mu_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_up[0],err_lep_up[1],err_lep_up[2],err_lep_up[3]))
               uncor_4mu_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],err_lep_trig_dn[2],err_lep_trig_dn[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_dn[0],err_lep_dn[1],err_lep_dn[2],err_lep_dn[3])) 

               for i in range(4):
                  sf_ratio2 = (err_lep_up[i] / SF_lep[i])**2
                  if sf_ratio2 == 0.25:
                     plot_Pt_mu.append(tree.LepPt[i])
                     plot_Eta_mu.append(tree.LepEta[i])
                     plot_sf_ratio2_mu.append(sf_ratio2)

            yield_4mu_up[k] += weight_nom/SF_tot_nom * SF_var_up
            yield_4mu_dn[k] += weight_nom/SF_tot_nom * SF_var_dn

         elif (abs(idL1-idL3)==2):
            # Vary electron SF while fixing muon and vice-versa
            if ( idL1 == 11):
               if (correlated_leptons):
                  SF_var_e_up  = (SF_lep_trig[0] + TRIG*err_lep_trig_up[0]) * (SF_lep[0] + RECO_SEL*err_lep_up[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_up[1]) * (SF_lep[1] + RECO_SEL*err_lep_up[1]) * SF_lep_trig[2] * SF_lep[2] * SF_lep_trig[3] * SF_lep[3]
                  SF_var_mu_up = (SF_lep_trig[2] + TRIG*err_lep_trig_up[2]) * (SF_lep[2] + RECO_SEL*err_lep_up[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_up[3]) * (SF_lep[3] + RECO_SEL*err_lep_up[3]) * SF_lep_trig[0] * SF_lep[0] * SF_lep_trig[1] * SF_lep[1]

                  SF_var_e_dn  = (SF_lep_trig[0] + TRIG*err_lep_trig_dn[0]) * (SF_lep[0] + RECO_SEL*err_lep_dn[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_dn[1]) * (SF_lep[1] + RECO_SEL*err_lep_dn[1]) * SF_lep_trig[2] * SF_lep[2] * SF_lep_trig[3] * SF_lep[3]
                  SF_var_mu_dn = (SF_lep_trig[2] + TRIG*err_lep_trig_dn[2]) * (SF_lep[2] + RECO_SEL*err_lep_dn[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_dn[3]) * (SF_lep[3] + RECO_SEL*err_lep_dn[3]) * SF_lep_trig[0] * SF_lep[0] * SF_lep_trig[1] * SF_lep[1]
               
               if (uncorrelated_leptons):

                  uncor_2e2mu_e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],0,0)) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_up[0],err_lep_up[1],0,0))
                  uncor_2e2mu_e_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],0,0)) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_dn[0],err_lep_dn[1],0,0))
                  uncor_2e2mu_mu_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],0,0,err_lep_up[2],err_lep_up[3]))
                  uncor_2e2mu_mu_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_dn[2],err_lep_trig_dn[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],0,0,err_lep_dn[2],err_lep_dn[3]))
            
            elif ( idL1 == 13):
               if ( correlated_leptons ):
                  SF_var_mu_up  = (SF_lep_trig[0] + TRIG*err_lep_trig_up[0]) * (SF_lep[0] + RECO_SEL*err_lep_up[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_up[1]) * (SF_lep[1] + RECO_SEL*err_lep_up[1]) * SF_lep_trig[2] * SF_lep[2] * SF_lep_trig[3] * SF_lep[3]
                  SF_var_e_up = (SF_lep_trig[2] + TRIG*err_lep_trig_up[2]) * (SF_lep[2] + RECO_SEL*err_lep_up[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_up[3]) * (SF_lep[3] + RECO_SEL*err_lep_up[3]) * SF_lep_trig[0] * SF_lep[0] * SF_lep_trig[1] * SF_lep[1]

                  SF_var_mu_dn  = (SF_lep_trig[0] + TRIG*err_lep_trig_dn[0]) * (SF_lep[0] + RECO_SEL*err_lep_dn[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_dn[1]) * (SF_lep[1] + RECO_SEL*err_lep_dn[1]) * SF_lep_trig[2] * SF_lep[2] * SF_lep_trig[3] * SF_lep[3]
                  SF_var_e_dn = (SF_lep_trig[2] + TRIG*err_lep_trig_dn[2]) * (SF_lep[2] + RECO_SEL*err_lep_dn[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_dn[3]) * (SF_lep[3] + RECO_SEL*err_lep_dn[3]) * SF_lep_trig[0] * SF_lep[0] * SF_lep_trig[1] * SF_lep[1]
                     
               if (uncorrelated_leptons):
                  uncor_2e2mu_mu_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],0,0)) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_up[0],err_lep_up[1],0,0))
                  uncor_2e2mu_mu_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],0,0)) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_dn[0],err_lep_dn[1],0,0))
                  uncor_2e2mu_e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],0,0,err_lep_up[2],err_lep_up[3]))
                  uncor_2e2mu_e_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_dn[2],err_lep_trig_dn[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],0,0,err_lep_dn[2],err_lep_dn[3]))

            yield_2e2mu_e_up[k]  += weight_nom/SF_tot_nom * SF_var_e_up
            yield_2e2mu_mu_up[k] += weight_nom/SF_tot_nom * SF_var_mu_up
               
            yield_2e2mu_e_dn[k]  += weight_nom/SF_tot_nom * SF_var_e_dn
            yield_2e2mu_mu_dn[k] += weight_nom/SF_tot_nom * SF_var_mu_dn
      


   comb_4e_up = 0.
   comb_4mu_up = 0.
   comb_2e2mu_e_up = 0.
   comb_2e2mu_mu_up = 0.

   comb_4e_dn = 0.
   comb_4mu_dn = 0.
   comb_2e2mu_e_dn = 0.
   comb_2e2mu_mu_dn = 0.
   
   comb_uncor_4e_up = 0.
   comb_uncor_4mu_up = 0.
   comb_uncor_2e2mu_e_up = 0.
   comb_uncor_2e2mu_mu_up = 0.
   
   comb_uncor_4e_dn = 0.
   comb_uncor_4mu_dn = 0.
   comb_uncor_2e2mu_e_dn = 0.
   comb_uncor_2e2mu_mu_dn = 0.

   # Sum relative uncertainties in quadrature
   for k in range (0,2):
      if ( correlated_leptons ):
         comb_4e_up += ((yield_4e_up[k]-nom_yield_4e)/nom_yield_4e*100)**2
         comb_4mu_up += ((yield_4mu_up[k]-nom_yield_4mu)/nom_yield_4mu*100)**2
         comb_2e2mu_e_up += ((yield_2e2mu_e_up[k]-nom_yield_2e2mu)/nom_yield_2e2mu*100)**2
         comb_2e2mu_mu_up += ((yield_2e2mu_mu_up[k]-nom_yield_2e2mu)/nom_yield_2e2mu*100)**2
         
         comb_4e_dn += ((-yield_4e_dn[k]+nom_yield_4e)/nom_yield_4e*100)**2
         comb_4mu_dn += ((-yield_4mu_dn[k]+nom_yield_4mu)/nom_yield_4mu*100)**2
         comb_2e2mu_e_dn += ((-yield_2e2mu_e_dn[k]+nom_yield_2e2mu)/nom_yield_2e2mu*100)**2
         comb_2e2mu_mu_dn += ((-yield_2e2mu_mu_dn[k]+nom_yield_2e2mu)/nom_yield_2e2mu*100)**2
         
      if (uncorrelated_leptons):
         comb_uncor_4e_up += (uncor_4e_up[k]/br_fs4e)
         comb_uncor_4mu_up += (uncor_4mu_up[k]/br_fs4mu)

         #print("uncor_4mu_up[k]", uncor_4mu_up[k])
         #print('br_fs4mu',br_fs4mu)
         #print("comb_uncor_4mu_up", comb_uncor_4mu_up)

         comb_uncor_2e2mu_e_up += (uncor_2e2mu_e_up[k]/br_fs2e2mu)
         comb_uncor_2e2mu_mu_up += (uncor_2e2mu_mu_up[k]/br_fs2e2mu)
         
         comb_uncor_4e_dn += (uncor_4e_dn[k]/br_fs4e)
         comb_uncor_4mu_dn += (uncor_4mu_dn[k]/br_fs4mu)
         comb_uncor_2e2mu_e_dn += (uncor_2e2mu_e_dn[k]/br_fs2e2mu)
         comb_uncor_2e2mu_mu_dn += (uncor_2e2mu_mu_dn[k]/br_fs2e2mu)
      
      # Print the results
      if ( correlated_leptons ):
         print("==================================================================")
         print("| CORRELATED LEPTON UNCERTAINTIES FOR {} IN {} SAMPLE |".format(variation_name[k],Type))
         print("==================================================================") 
         print("4e = +{:.1f} % -{:.1f}%  4mu = +{:.1f}% -{:.1f} %".format(((yield_4e_up[k]-nom_yield_4e)/nom_yield_4e*100),((-yield_4e_dn[k]+nom_yield_4e)/nom_yield_4e*100),((yield_4mu_up[k]-nom_yield_4mu)/nom_yield_4mu*100),((-yield_4mu_dn[k]+nom_yield_4mu)/nom_yield_4mu*100)))
         print("2e2mu/ele = +{:.1f} % -{:.1f} % 2e2mu/mu = +{:.1f}% -{:.1f}%".format(((yield_2e2mu_e_up[k]-nom_yield_2e2mu)/nom_yield_2e2mu*100),((-yield_2e2mu_e_dn[k]+nom_yield_2e2mu)/nom_yield_2e2mu*100),((yield_2e2mu_mu_up[k]-nom_yield_2e2mu)/nom_yield_2e2mu*100),((-yield_2e2mu_mu_dn[k]+nom_yield_2e2mu)/nom_yield_2e2mu*100)))
   
      if (uncorrelated_leptons):
         print("==================================================================")
         print("| UNCORRELATED LEPTON UNCERTAINTIES FOR {} IN {} SAMPLE |".format(variation_name[k],Type))
         print("==================================================================")
         print("4e = +{:.1f}% -{:.1f}%  4mu = +{:.1f}% -{:.1f}%".format(100*(uncor_4e_up[k]/br_fs4e)**0.5,100*(uncor_4e_dn[k]/br_fs4e)**0.5,100*(uncor_4mu_up[k]/br_fs4mu)**0.5,100*(uncor_4mu_dn[k]/br_fs4mu)**0.5))
         print("2e2mu/ele = +{:.1f}% -{:.1f}% 2e2mu/mu = +{:.1f}% -{:.1f}%".format(100*(uncor_2e2mu_e_up[k]/br_fs2e2mu)**0.5,100*(uncor_2e2mu_e_dn[k]/br_fs2e2mu)**0.5,100*(uncor_2e2mu_mu_up[k]/br_fs2e2mu)**0.5,100*(uncor_2e2mu_mu_dn[k]/br_fs2e2mu)**0.5))

   if ( correlated_leptons ):
      print("===============================")
      print("| CORRELATED COMBINATION |")
      print("===============================")
      print("4e = +{:.1f}% -{:.1f}%  4mu = +{:.1f}% -{:.1f}%".format(comb_4e_up**0.5,comb_4e_dn**0.5,comb_4mu_up**0.5,comb_4mu_dn**0.5))
      print("2e2mu/ele = +{:.1f}% -{:.1f}% 2e2mu/mu = +{:.1f}% -{:.1f}%".format(comb_2e2mu_e_up**0.5,comb_2e2mu_e_dn**0.5,comb_2e2mu_mu_up**0.5,comb_2e2mu_mu_dn**0.5))
      print("==============================================================================")

   if (uncorrelated_leptons):
      print("comb_uncor_4mu_up", comb_uncor_4mu_up)
      print("comb_uncor_4mu_up", comb_uncor_4mu_up**0.5)
      print("================================")
      print("| UNCORRELATED COMBINATION |")
      print("================================")
      print("4e = +{:.1f}% -{:.1f}%  4mu = +{:.1f}% -{:.1f}%".format(100*(comb_uncor_4e_up**0.5),100*(comb_uncor_4e_dn**0.5),100*(comb_uncor_4mu_up**0.5),100*(comb_uncor_4mu_dn**0.5)))
      print("2e2mu/ele = +{:.1f}% -{:.1f}% 2e2mu/mu = +{:.1f}% -{:.1f}%".format(100*(comb_uncor_2e2mu_e_up**0.5),100*(comb_uncor_2e2mu_e_dn**0.5),100*(comb_uncor_2e2mu_mu_up**0.5),100*(comb_uncor_2e2mu_mu_dn**0.5)))
      print("==============================================================================")