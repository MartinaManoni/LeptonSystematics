import sys, getopt
import os
import glob
import ROOT
from math import sqrt
import time
from pathlib import Path
from tqdm import trange, tqdm
import numpy as np
from array import *
from collections import Counter
from decimal import *
import uproot
from functools import reduce
import itertools
import argparse
import awkward as ak 
#from root_numpy import array2tree, tree2array
from AnalysisTools.data import gConstants as gConstants
from AnalysisTools.data import cConstants as cConstants
from AnalysisTools.Utils import Config as Config
from AnalysisTools.Utils import OnShell_Category as OnShell_Category
from AnalysisTools.Utils import Discriminants_Vec as Discriminants
from AnalysisTools.Utils import OnShell_Help as OnShell_Help
from AnalysisTools.Utils.Variable_Names_For_Tagging import * 
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
import vector 

basePath = "/afs/cern.ch/work/g/gritsan/public_write/nipinto/CMSSW_14_1_6/src/ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/"
f_eleID          = basePath+"SF2022eleID_postEE.root"
f_eleReco_highPt = "/eos/user/n/nipinto/old_CMSSW_13_3_3/src/On_Shell_Analysis_Run3/AnalysisTools/data/LepSF/egammaEffi.txt_EGM2D.root"
f_eleReco_midPt  = "/eos/user/n/nipinto/old_CMSSW_13_3_3/src/On_Shell_Analysis_Run3/AnalysisTools/data/LepSF/egammaEffi.txt_EGM2D.root"
f_eleReco_lowPt  = "/eos/user/n/nipinto/old_CMSSW_13_3_3/src/On_Shell_Analysis_Run3/AnalysisTools/data/LepSF/egammaEffi.txt_EGM2D.root"


root_file = ROOT.TFile.Open(f_eleID, "READ")
h_Ele_ID = root_file.Get("EGamma_SF2D").Clone("h_Ele_ID")
h_Ele_ID.SetDirectory(0)
root_file.Close()


root_file = ROOT.TFile.Open(f_eleReco_highPt,"READ")
h_Ele_Reco_highPt = root_file.Get("EGamma_SF2D").Clone("h_Ele_Reco_highPt")
h_Ele_Reco_highPt.SetDirectory(0)
root_file.Close()

root_file = ROOT.TFile.Open(f_eleReco_lowPt,"READ")
h_Ele_Reco_lowPt = root_file.Get("EGamma_SF2D").Clone("h_Ele_Reco_lowPt")
h_Ele_Reco_lowPt.SetDirectory(0)
root_file.Close()

root_file = ROOT.TFile.Open(f_eleReco_midPt,"READ")
h_Ele_Reco_midPt = root_file.Get("EGamma_SF2D").Clone("h_Ele_Reco_midPt")
h_Ele_Reco_midPt.SetDirectory(0)
root_file.Close()


#MUONS
f_mu = basePath+"final_HZZ_SF_2023D_RMS_mupogsysts.root"
root_file = ROOT.TFile.Open(f_mu,"READ")
h_Mu_SF  = root_file.Get("FINAL").Clone("h_Mu_SF")
h_Mu_Unc = root_file.Get("ERROR").Clone("h_Mu_Unc")
h_Mu_SF.SetDirectory(0)
h_Mu_Unc.SetDirectory(0)
root_file.Close()

### PDF UNCERTAINTY CALCULATION: 
def GetLHEPDFUncertainty(tree): 
    LHEPDFWeights = tree.arrays("LHEPdfWeight")["LHEPdfWeight"]
    Sig0 = LHEPDFWeights[:,0]
    diff = []
    diff = ak.Array([(LHEPDFWeights[:,i] - Sig0)**2 for i in np.arange(101)])
    newDiff = np.transpose(diff)
    PdfUnc = np.sqrt(ak.sum(newDiff, axis = 1))
    AsDn = LHEPDFWeights[:,101]
    AsUp = LHEPDFWeights[:,102]
    return PdfUnc, AsDn, AsUp

def GetLHEScaleVariations(tree): 
    LHEScaleWeights = tree.arrays("LHEScaleWeight")["LHEScaleWeight"]
    muF_0p5_muR_0p5 = LHEScaleWeights[:,0]
    muF_1_muR_0p5 = LHEScaleWeights[:,1]
    muF_0p5_muR_1 = LHEScaleWeights[:,3]
    muF_1_muR_1 = LHEScaleWeights[:,4]
    muF_2_muR_1 = LHEScaleWeights[:,5]
    muF_1_muR_2 = LHEScaleWeights[:,7]
    muF_2_muR_2 = LHEScaleWeights[:,8]
    return muF_0p5_muR_0p5, muF_1_muR_0p5, muF_0p5_muR_1, muF_1_muR_1, muF_2_muR_1, muF_1_muR_2, muF_2_muR_2

### GGZZ QCD KFactor Determination: 
def evalSpline(sp, xval):
        xmin = sp.GetXmin()
        xmax = sp.GetXmax()
        res = 0.
        if xval<xmin :
            res=sp.Eval(xmin)
            deriv=sp.Derivative(xmin)
            res += deriv*(xval-xmin)
        elif (xval>xmax) :
            res=sp.Eval(xmax)
            deriv=sp.Derivative(xmax)
            res += deriv*(xval-xmax)
        else :
            res=sp.Eval(xval)
        return res

# GGZZ QCD SF Setup: 
strZZGGKFVar = ["Nominal", "PDFScaleDn", "PDFScaleUp", "QCDScaleDn", "QCDScaleUp", "AsDn", "AsUp", "PDFReplicaDn", "PDFReplicaUp"]
# NNLO
ggZZKFactorFile = ROOT.TFile.Open('/eos/user/n/nipinto/old_CMSSW_13_3_3/src/On_Shell_Analysis_Run3/AnalysisTools/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root')
spkfactor_ggzz_nnlo = [None]*9
spkfactor_ggzz_nlo = [None]*9
for i in range(0,9):
    for i in range(0,9):
        spkfactor_ggzz_nnlo[i] = (ggZZKFactorFile.Get('sp_kfactor_%s' % strZZGGKFVar[i])).Clone('sp_kfactor_%s_NNLO' % strZZGGKFVar[i])
ggZZKFactorFile.Close()
# NLO
ggZZKFactorFile = ROOT.TFile.Open('/eos/user/n/nipinto/old_CMSSW_13_3_3/src/On_Shell_Analysis_Run3/AnalysisTools/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NLO_NNPDF_NarrowWidth_13TeV.root')
for i in range(0,9):
    for i in range(0,9):
        spkfactor_ggzz_nlo[i] = (ggZZKFactorFile.Get('sp_kfactor_%s' % strZZGGKFVar[i])).Clone('sp_kfactor_%s_NNLO' % strZZGGKFVar[i])
ggZZKFactorFile.Close()

def GetGGZZ_QCD_KFactor(tree, mass_string="GenZZ_mass", variation=0):
    print("**GetGGZZ_QCD_KFactor: Mass string is: ", mass_string)
    masses = tree.arrays(mass_string)[mass_string]
    ### SHOULD USE GEN LEVEL ZZ MASS to evaluate KFACTOR!!!
    KFactor_QCD_ggZZ_NNLO_over_LO = np.zeros(len(masses))
    KFactor_QCD_ggZZ_NLO_over_LO = np.zeros(len(masses))
    KFactor_QCD_ggZZ_NNLO_over_NLO = np.zeros(len(masses))
    for i in range(len(masses)):
        KFactor_QCD_ggZZ_NNLO_over_LO[i] = evalSpline(spkfactor_ggzz_nnlo[variation], masses[i]) #1: NNLO/LO
        #KFactor_QCD_ggZZ_NLO_over_LO[i] = evalSpline(spkfactor_ggzz_nlo[variation], masses[i]) #2: NLO/LO
        #KFactor_QCD_ggZZ_NNLO_over_NLO[i] = KFactor_QCD_ggZZ_NNLO_over_LO[i] / KFactor_QCD_ggZZ_NLO_over_LO[i] #3: NNLO/LO 
    return KFactor_QCD_ggZZ_NNLO_over_LO #, KFactor_QCD_ggZZ_NLO_over_LO, KFactor_QCD_ggZZ_NNLO_over_NLO

def GetQQZZ_EW_KFactor_Unc(tree, mass_string="GenZZ_mass"): 
    print("**GetQQZZ_EW_KFactor_Unc")
    Gen_pdgIds = tree.arrays("GenPart_pdgId")["GenPart_pdgId"]
    Gen_pts = tree.arrays("GenPart_pt")["GenPart_pt"]
    Gen_etas = tree.arrays("GenPart_eta")["GenPart_eta"]
    Gen_phi = tree.arrays("GenPart_phi")["GenPart_phi"]
    Gen_masses = tree.arrays("GenPart_mass")["GenPart_mass"]

    Z1_L1_Idxs = tree.arrays("GenZZ_Z1l1Idx")["GenZZ_Z1l1Idx"]
    Z1_L2_Idxs = tree.arrays("GenZZ_Z1l2Idx")["GenZZ_Z1l2Idx"]
    Z2_L1_Idxs = tree.arrays("GenZZ_Z2l1Idx")["GenZZ_Z2l1Idx"]
    Z2_L2_Idxs = tree.arrays("GenZZ_Z2l2Idx")["GenZZ_Z2l2Idx"]

    Z1_L1_pt = Gen_pts[ak.local_index(Gen_pts, axis = 0), Z1_L1_Idxs]
    Z1_L2_pt = Gen_pts[ak.local_index(Gen_pts, axis = 0), Z1_L2_Idxs]
    Z2_L1_pt = Gen_pts[ak.local_index(Gen_pts, axis = 0), Z2_L1_Idxs]
    Z2_L2_pt = Gen_pts[ak.local_index(Gen_pts, axis = 0), Z2_L2_Idxs]

    Z1_L1_eta = Gen_etas[ak.local_index(Gen_etas, axis = 0), Z1_L1_Idxs]
    Z1_L2_eta = Gen_etas[ak.local_index(Gen_etas, axis = 0), Z1_L2_Idxs]
    Z2_L1_eta = Gen_etas[ak.local_index(Gen_etas, axis = 0), Z2_L1_Idxs]
    Z2_L2_eta = Gen_etas[ak.local_index(Gen_etas, axis = 0), Z2_L2_Idxs]

    Z1_L1_pt = Gen_pts[ak.local_index(Gen_pts, axis = 0), Z1_L1_Idxs]
    Z1_L2_pt = Gen_pts[ak.local_index(Gen_pts, axis = 0), Z1_L2_Idxs]
    Z2_L1_pt = Gen_pts[ak.local_index(Gen_pts, axis = 0), Z2_L1_Idxs]
    Z2_L2_pt = Gen_pts[ak.local_index(Gen_pts, axis = 0), Z2_L2_Idxs]

def getSF( flav,  pt,  eta,  SCeta,  isCrack): 
    


    RecoSF = 1.0
    SelSF = 1.0
    SF = 1.0

    RecoSF_Unc = 0.0
    SelSF_Unc = 0.0
    SFError = 0.0
   

   

    if(abs(flav) == 11):
        if(pt < 20.):
            RecoSF     = h_Ele_Reco_lowPt.GetBinContent(h_Ele_Reco_lowPt.GetXaxis().FindBin(SCeta),h_Ele_Reco_lowPt.GetYaxis().FindBin(15.))
            RecoSF_Unc = h_Ele_Reco_lowPt.GetBinError  (h_Ele_Reco_lowPt.GetXaxis().FindBin(SCeta),h_Ele_Reco_lowPt.GetYaxis().FindBin(15.))
        elif (pt < 75. and h_Ele_Reco_midPt!= 0):
            RecoSF     = h_Ele_Reco_midPt.GetBinContent(h_Ele_Reco_midPt.GetXaxis().FindBin(SCeta),h_Ele_Reco_midPt.GetYaxis().FindBin(min(pt,75)))
            RecoSF_Unc = h_Ele_Reco_midPt.GetBinError  (h_Ele_Reco_midPt.GetXaxis().FindBin(SCeta),h_Ele_Reco_midPt.GetYaxis().FindBin(min(pt,75)))
        
        
        else:
            RecoSF     = h_Ele_Reco_highPt.GetBinContent(h_Ele_Reco_highPt.GetXaxis().FindBin(SCeta),h_Ele_Reco_highPt.GetYaxis().FindBin(min(pt,499)))
            RecoSF_Unc = h_Ele_Reco_highPt.GetBinError  (h_Ele_Reco_highPt.GetXaxis().FindBin(SCeta),h_Ele_Reco_highPt.GetYaxis().FindBin(min(pt,499)))
        
        
        

        SF = RecoSF*SelSF   
        SFError = sqrt( RecoSF_Unc*RecoSF_Unc/(RecoSF*RecoSF) + SelSF_Unc*SelSF_Unc/(SelSF*SelSF) )
    

    if(abs(flav) == 13 ):
        SelSF = h_Mu_SF.GetBinContent(h_Mu_SF.GetXaxis().FindBin(eta),h_Mu_SF.GetYaxis().FindBin(min(pt,199)))
        SelSF_Unc = h_Mu_Unc.GetBinContent(h_Mu_Unc.GetXaxis().FindBin(eta),h_Mu_Unc.GetYaxis().FindBin(min(pt,199)))
        
        SF = SelSF
        SFError = SelSF_Unc/SelSF
    

    return SF, SFError


def main(raw_args=None): 
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--ifile', type=str, required=True, help="Input text file")
    parser.add_argument("-o", "--odir", type=str, required=True, help=" Directory to dump the output root file")
    parser.add_argument("-t", "--tree", type=str, default="eventTree", help= " Name of tree in root file")
    args = parser.parse_args(raw_args)

    inputText = args.ifile
    outputDir = args.odir
    treeName = args.tree


    if not outputDir.endswith("/"): 
        outputDir = outputDir+"/"

    input_trees = []

    outDict = {}
    
    



    

    # Parse text file for sample categorization: 
    with open(inputText) as file: 
        for line in file: 
            era = line.split("/")[-3]
            sample = line.split("/")[-2]
            rootFile_name = line.split("/")[-1].rstrip()
            SubFolder_Name = era + "/" + sample
            path_to_sample = line.rstrip()
            print(path_to_sample)
            print(SubFolder_Name)
            print(rootFile_name)
            t = uproot.open(path_to_sample)[treeName]

            if not os.path.exists(outputDir+SubFolder_Name):
                os.makedirs(outputDir+SubFolder_Name, mode=0o777, exist_ok=True)
            isCrack = False

            outDict["ZZCand_dataMCWeight"] = np.zeros(len(t["LepLepId"]))

            LepIDX_list = ak.concatenate([t["ZZCand_Z1l1Idx"].array(), t["ZZCand_Z1l2Idx"].array(), t["ZZCand_Z2l1Idx"].array(), t["ZZCand_Z2l2Idx"].array()], axis = 1)
            
            LepIDs = t["Lepton_pdgId"].array()
            LepPts = t["Lepton_pt"].array()
            LepEtas = t["Lepton_eta"].array()
            LepPhis = t['Lepton_phi'].array()
            Electron_deltaEtaSC = t["Electron_deltaEtaSC"].array()
            for i in tqdm(np.arange(len(t["LepLepId"].array())), miniters=int(len(t["LepLepId"].array())/100)): 
            # for i in trange(100):
                SF = 1 
                SFerror = 0 
                for k in np.arange(4):
                    lepID = abs(LepIDs[i][LepIDX_list[i][k]])
                    lepPt = LepPts[i][LepIDX_list[i][k]]
                    lepEta = LepEtas[i][LepIDX_list[i][k]]
                    lepPhi = LepPhis[i][LepIDX_list[i][k]]
                    mySCeta = lepEta
                    
                    if lepID == 11: 
                        # print("I HERE", i)
                        # print("K HERE", k)
                        mySCeta = lepEta + Electron_deltaEtaSC[i][LepIDX_list[i][k]]
                    
                    mySCeta = min(mySCeta,2.49)
                    mySCeta = max(mySCeta,-2.49)

                    pair = getSF(lepID, lepPt, lepEta, mySCeta, isCrack)
                    SF *= pair[0]
                    SFerror += pair[1]
                

                    if SF == 0: 
                        SF, SFerror = 1, 0.5
                
                outDict["ZZCand_dataMCWeight"] = np.append(outDict["ZZCand_dataMCWeight"], [SF])
                
                

            


            

            ### Write everything out: 
            N_total = t.num_entries
            # N_total = 100
            entries_per_chunk = 10000

            with uproot.recreate(outputDir+SubFolder_Name+"/"+rootFile_name) as fout: 
                first = True
                pbar = tqdm(total=N_total, desc=f"Writing combined tree for {sample}", unit="events")
                for start in range(0, N_total, entries_per_chunk):
                    stop = min(N_total, start + entries_per_chunk)
                    adding_chunk = t.arrays(entry_start=start, entry_stop=stop, library="np")
                    new_chunk = {}

                    for branch, array in outDict.items():
                        # Convert awkward array to numpy if needed, then slice
                        if isinstance(array, ak.Array):
                            new_chunk[branch] = array[start:stop]
                        else:
                            new_chunk[branch] = array[start:stop]

                    # Combine dictionaries for this chunk
                    combined_chunk = {**adding_chunk, **new_chunk}

                    if first:
                        fout[treeName] = combined_chunk
                        first = False
                    else:
                        fout[treeName].extend(combined_chunk)

                    pbar.update(stop - start)

                pbar.close()


                    # fout[treeName] = combined

                
            





if __name__ == "__main__":
    main(sys.argv[1:])


