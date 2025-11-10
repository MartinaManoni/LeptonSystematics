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

gRandom.SetSeed(101)

################################################################################
# ANALYSIS CONFIGURATION
################################################################################

year = "2022"

if (year=="2022"):
   lumi= 7.98
   folder = '/eos/user/m/mmanoni/HZZ_prod_031125_22_23/2022_MC/PROD_samplesNano_2022_MC_a3d06bb6/'
   #'/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-24-013/RunIII_byZ1Z2/240820/2022/'#
   #"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2022_MC/" # Define input folder name
#file_name = '/ZZ4lAnalysis_SKIMMED.root' # Define input file name
file_name = '/ZZ4lAnalysis_LepSyst.root'
#file_name = '/ZZ4lAnalysis_SKIMMED.root'
#ZZ4lAnalysis_SKIMMED.root

# Flags to control uncertainty treatment
correlated_leptons = True     # treat lepton SF uncertainties as correlated across all 4 leptons
uncorrelated_leptons = False   # treat lepton SF uncertainties as independent (quadrature sum)
corr_factor = 0               # correlation coefficient (ρ = 0 means no correlation)

# List of samples to run on
List = [
#'ggH125'
'ggH125_LepSyst',
#'ZZTo4l',
# 'ggTo4l'
]

print("List of samples", List)


################################################################################
# HELPER FUNCTION: combine uncertainties across leptons with correlation ρ
################################################################################

def sigma_event(rho, SF1, SF2, SF3, SF4, sigma1, sigma2, sigma3, sigma4):
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
   tree.Draw("ZZMass >> h_nom" , "(abs(LepLepId[0]) == 11 && abs(LepLepId[3]) == 11)*overallEventWeight*1000*"+str(lumi)+"/"+str(NGen))
   nom_yield_4e = h_nom.Integral()
   print("nom_yield_4e", nom_yield_4e)
   
   tree.Draw("ZZMass >> h_nom" , "(abs(LepLepId[0]) == 13 && abs(LepLepId[3]) == 13)*overallEventWeight*1000*"+str(lumi)+"/"+str(NGen))
   nom_yield_4mu = h_nom.Integral()
   print("nom_yield_4mu", nom_yield_4mu)
   
   tree.Draw("ZZMass >> h_nom" , "(abs(abs(LepLepId[0]) - (LepLepId[3])) == 2)*overallEventWeight*1000*"+str(lumi)+"/"+str(NGen))
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

   for event in tree:# Loop over all events in tree
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
      
      # Calculate nominal weigh using central value of SF
      weight_nom = event.overallEventWeight * 1000 * lumi / NGen #xsec is already inside overallEventWeight --- what about L1prefiringWeight?
      
      # Nominal value of total SF, product of 4 lepton nominal SF
      SF_tot_nom = event.dataMCWeight
      
      SF_lep_trig = []
      err_lep_trig_up = []
      err_lep_trig_dn = []

      SF_lep = []
      err_lep_up = []
      err_lep_dn = []
      
      for i in range (0,4):
         # print(err_lep_trig_dn)
         SF_lep_trig.append(0.)
         # err_lep_trig_up.append(0.)
         # err_lep_trig_dn.append(0.)
     
         SF_lep.append(0.)
         err_lep_up.append(0.)
         err_lep_dn.append(0.)
         
      for i in range (0,4):
         # Hard-coded trigger SF and unc (Since Run2)
         err_lep_trig_dn = [0.01,0,0,0]
         SF_lep_trig[i] = 1.
 
         if(idL1==11 and idL3==11 and event.LepPt[3] < 12):
               err_lep_trig_up = [0.02,0,0,0]
               err_lep_trig_dn = [0.11,0,0,0]
         if(idL1==11 and idL3==11 and event.LepPt[3] >= 12):
               err_lep_trig_up = [0.005,0,0,0]
               err_lep_trig_dn = [0.01,0,0,0]
         if(idL1==13 and idL3==13 and event.LepPt[3] < 7):
               err_lep_trig_up = [0.001,0,0,0]
               err_lep_trig_dn = [0.032,0,0,0]
         if(idL1==13 and idL3==13 and event.LepPt[3] >= 7 and event.LepPt[3] < 12):
               err_lep_trig_up = [0.001,0,0,0]
               err_lep_trig_dn = [0.015,0,0,0]
         if(idL1==13 and idL3==13 and event.LepPt[3] >= 12):
               err_lep_trig_up = [0.001,0,0,0]
               err_lep_trig_dn = [0.015,0,0,0]
                  
         if(abs(idL1-idL3)==2 and event.LepPt[3] < 7):
               err_lep_trig_up = [0.005,0,0,0]
               err_lep_trig_dn = [0.08,0,0,0]
         if(abs(idL1-idL3)==2 and event.LepPt[3] >= 7 and event.LepPt[3] < 12):
               err_lep_trig_up = [0.005,0,0,0]
               err_lep_trig_dn = [0.032,0,0,0]
         if(abs(idL1-idL3)==2 and event.LepPt[3] >= 12):
               err_lep_trig_up = [0.001,0,0,0]
               err_lep_trig_dn = [0.01,0,0,0]
       
         # Load reco and selection SF and unc
         # for electrons, reconstrucion SF unc and HZZ selection SF unc are combined and provided in the branch stored in our NTuple;
         # print(err_lep_trig_up)
         if (abs(tree.LepLepId[i]) == 11):

            SF_lep[i] = event.Electron_SF[i]
            err_lep_up[i] = event.Electron_SFUnc[i]
            err_lep_dn[i] = event.Electron_SFUnc[i]

            #print("ELE SF: ", SF_lep, " ##########")
            #print("ELE: err_lep_up", err_lep_up)
            #print("ELE: err_lep_dn", err_lep_dn)

         elif (abs(tree.LepLepId[i]) == 13):
            SF_lep[i] = event.Muon_SF[i]
            err_lep_up[i] = event.Muon_SFUnc[i]
            err_lep_dn[i] = event.Muon_SFUnc[i]

            #print("MUON SF: ", SF_lep, " ##########")
            #print("MUON: err_lep_up", err_lep_up)
            #print("MUON: err_lep_dn", err_lep_dn)
            
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
               if ( correlated_leptons ):
                  SF_var_up *= (SF_lep_trig[i] + TRIG*err_lep_trig_up[i]) * (SF_lep[i] + RECO_SEL*err_lep_up[i])
                  SF_var_dn *= (SF_lep_trig[i] - TRIG*err_lep_trig_dn[i]) * (SF_lep[i] - RECO_SEL*err_lep_dn[i])
           
            if (uncorrelated_leptons):
               uncor_4e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_up[0],err_lep_up[1],err_lep_up[2],err_lep_up[3]))
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

            yield_4mu_up[k] += weight_nom/SF_tot_nom * SF_var_up
            yield_4mu_dn[k] += weight_nom/SF_tot_nom * SF_var_dn

         elif (abs(idL1-idL3)==2):
            # Vary electron SF while fixing muon and vice-versa
            if ( idL1 == 11):
               if ( correlated_leptons ):
                  SF_var_e_up  = (SF_lep_trig[0] + TRIG*err_lep_trig_up[0]) * (SF_lep[0] + RECO_SEL*err_lep_up[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_up[1]) * (SF_lep[1] + RECO_SEL*err_lep_up[1]) * SF_lep_trig[2] * SF_lep[2] * SF_lep_trig[3] * SF_lep[3]
                  SF_var_mu_up = (SF_lep_trig[2] + TRIG*err_lep_trig_up[2]) * (SF_lep[2] + RECO_SEL*err_lep_up[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_up[3]) * (SF_lep[3] + RECO_SEL*err_lep_up[3]) * SF_lep_trig[0] * SF_lep[0] * SF_lep_trig[1] * SF_lep[1]

                  SF_var_e_dn  = (SF_lep_trig[0] + TRIG*err_lep_trig_dn[0]) * (SF_lep[0] + RECO_SEL*err_lep_dn[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_dn[1]) * (SF_lep[1] + RECO_SEL*err_lep_dn[1]) * SF_lep_trig[2] * SF_lep[2] * SF_lep_trig[3] * SF_lep[3]
                  SF_var_mu_dn = (SF_lep_trig[2] + TRIG*err_lep_trig_dn[2]) * (SF_lep[2] + RECO_SEL*err_lep_dn[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_dn[3]) * (SF_lep[3] + RECO_SEL*err_lep_dn[3]) * SF_lep_trig[0] * SF_lep[0] * SF_lep_trig[1] * SF_lep[1]
               
               if (uncorrelated_leptons):
                  uncor_2e2mu_e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],0,0)) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_up[0],err_lep_up[1],0,0))
                  uncor_2e2mu_e_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],0,0)) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_dn[0],err_lep_dn[1],0,0))
                  uncor_2e2mu_mu_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],0,0,err_lep_up[2],err_lep_dn[3]))
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
                  uncor_2e2mu_e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],0,0,err_lep_up[2],err_lep_dn[3]))
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
         print("uncor_4e_up[k]", uncor_4e_up[k])
         print('br_fs4e',br_fs4e)
         comb_uncor_4mu_up += (uncor_4mu_up[k]/br_fs4mu)
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
      print("================================")
      print("| UNCORRELATED COMBINATION |")
      print("================================")
      print("4e = +{:.1f}% -{:.1f}%  4mu = +{:.1f}% -{:.1f}%".format(100*(comb_uncor_4e_up**0.5),100*(comb_uncor_4e_dn**0.5),100*(comb_uncor_4mu_up**0.5),100*(comb_uncor_4mu_dn**0.5)))
      print("2e2mu/ele = +{:.1f}% -{:.1f}% 2e2mu/mu = +{:.1f}% -{:.1f}%".format(100*(comb_uncor_2e2mu_e_up**0.5),100*(comb_uncor_2e2mu_e_dn**0.5),100*(comb_uncor_2e2mu_mu_up**0.5),100*(comb_uncor_2e2mu_mu_dn**0.5)))
      print("==============================================================================")
