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

from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
import vector 


import os
basePath = f"{os.environ['CMSSW_BASE']}/src/ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/"

year = "2022"

if year == "2022":
    #ID
    f_eleID          = basePath+"SF2022eleID_postEE.root"
    #RECO
    f_eleReco_highPt = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptAbove75.txt_EGM2D_2022preEE.root"
    f_eleReco_midPt  = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptBelow75.txt_EGM2D_2022preEE.root"
    f_eleReco_lowPt  = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptBelow20.txt_EGM2D_2022preEE.root"
    #MUON
    f_mu = basePath+"final_HZZ_SF_Run3_2022_mupogsysts_newLoose_abseta3_fix_BCD_RMS.root"

elif year == "2022EE":
    #ID
    f_eleID           = basePath+"SF2022eleID_postEE.root"
    #RECO
    f_eleReco_highPt = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptAbove75.txt_EGM2D_2022postEE.root"
    f_eleReco_midPt  = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptBelow75.txt_EGM2D_2022postEE.root"
    f_eleReco_lowPt  = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptBelow20.txt_EGM2D_2022postEE.root"
    
    f_mu = basePath+"final_HZZ_SF_Run3_2022_mupogsysts_newLoose_abseta3_fix_EFG_RMS.root"

elif year == "2023preBPix":
    #ID
    f_eleID          = basePath+"SF2023eleID_preBPix.root"
    #RECO
    f_eleReco_highPt = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptAbove75.txt_EGM2D_2023preBPix.root"
    f_eleReco_midPt  = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptBelow75.txt_EGM2D_2023preBPix.root"
    f_eleReco_lowPt  = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptBelow20.txt_EGM2D_2023preBPix.root"
    #MUON
    f_mu = basePath+"final_HZZ_SF_2023C_RMS_mupogsysts.root"

elif year == "2023postBPix":
    #ID
    f_eleID          = basePath+"SF2023eleID_postBPix.root"
    f_eleID_HoleBPix = basePath+"SF2023eleID_postBPix_Hole.root"
    #RECO
    f_eleReco_highPt = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptAbove75.txt_EGM2D_2023postBPix.root"
    f_eleReco_midPt  = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptBelow75.txt_EGM2D_2023postBPix.root"
    f_eleReco_lowPt  = "/eos/user/m/mmanoni/RECO_ID_original/egammaEffi_ptBelow20.txt_EGM2D_2023postBPix.root"
    #MUON
    f_mu = basePath+"final_HZZ_SF_2023D_RMS_mupogsysts.root"


root_file = ROOT.TFile.Open(f_eleID, "READ")
h_Ele_ID = root_file.Get("EGamma_SF2D").Clone("h_Ele_ID")
h_Ele_ID.SetDirectory(0)
root_file.Close()

'''h_Ele_ID_HoleBPix = None
if f_eleID_HoleBPix != "":
    root_file = ROOT.TFile.Open(f_eleID_HoleBPix, "READ")
    h_Ele_ID_HoleBPix = root_file.Get("EGamma_SF2D").Clone("h_Ele_ID_HoleBPix")
    h_Ele_ID_HoleBPix.SetDirectory(0)
    root_file.Close()'''

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
root_file = ROOT.TFile.Open(f_mu,"READ")
h_Mu_SF  = root_file.Get("FINAL").Clone("h_Mu_SF")
h_Mu_Unc = root_file.Get("ERROR").Clone("h_Mu_Unc")
h_Mu_SF.SetDirectory(0)
h_Mu_Unc.SetDirectory(0)
root_file.Close()


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
        

        SelSF = h_Ele_ID.GetBinContent(h_Ele_ID.FindFixBin(SCeta, min(pt, 499)))
        SelSF_Unc = h_Ele_ID.GetBinError(h_Ele_ID.FindFixBin(SCeta, min(pt, 499)))

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
    #--odir /eos/user/m/mmanoni/ --tree Events
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

            outDict["ZZCand_dataMCWeight"] = []

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
                
                outDict["ZZCand_dataMCWeight"].append(SF)
            
            outDict["ZZCand_dataMCWeight"] = np.array(outDict["ZZCand_dataMCWeight"])
            
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