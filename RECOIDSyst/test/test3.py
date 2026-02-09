import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# =====================================================================
# File
# =====================================================================
file_path = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/ggH125/ZZ4lAnalysis.root"
file_path2 = "/eos/user/m/mmanoni/LeptSyst_RECOID_skimmed_samples/2024_ggH125/ZZ4lAnalysis_SKIMMED_LepSyst.root"

#"/eos/user/m/mmanoni/prod_ggH/2022/2024_skimmed.root"
#"/eos/user/m/mmanoni/LeptSyst_RECOID_skimmed_samples/2024_ggH125/ZZ4lAnalysis_SKIMMED_LepSyst.root"


# =====================================================================
# Read both trees with uproot (one method only)
# =====================================================================
f = uproot.open(file_path)
f2 = uproot.open(file_path2)

# Events tree
t_ev = f["Events"]
mu_sf_ev   = ak.flatten(t_ev["Muon_dataMC"].array())
mu_err_ev  = ak.flatten(t_ev["Muon_dataMCUnc"].array())
el_sf_ev   = ak.flatten(t_ev["Electron_dataMC"].array())
el_err_ev  = ak.flatten(t_ev["Electron_dataMCUnc"].array())

# ZZTree/candTree
t_zzt = f2["ZZTree/candTree"]
mu_sf_zzt   = ak.flatten(t_zzt["Muon_SF"].array())
mu_err_zzt  = ak.flatten(t_zzt["Muon_SFUnc"].array())
el_sf_zzt   = ak.flatten(t_zzt["Electron_SF"].array())
el_err_zzt  = ak.flatten(t_zzt["Electron_SFUnc"].array())

# =====================================================================
# Histogram overlay function
# =====================================================================
def overlay_hist(d1, d2, l1, l2, title, filename, bins=60):
    plt.figure(figsize=(8,6))
    plt.hist(d1, bins=bins, density=True, histtype="step", label=l1, linewidth=2)
    plt.hist(d2, bins=bins, density=True, histtype="step", label=l2, linewidth=2)
    plt.grid(True)
    plt.xlabel(title)
    plt.ylabel("Normalized counts")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# =====================================================================
# Make the 4 comparison plots
# =====================================================================

overlay_hist(mu_sf_ev, mu_sf_zzt,
             "Muon SF (Events tree)", "Muon SF (ZZTree/candTree)",
             "Muon SF Comparison", "MuonSF_compare.png")

overlay_hist(mu_err_ev, mu_err_zzt,
             "Muon SF err (Events tree)", "Muon SF err (ZZTree/candTree)",
             "Muon SF Uncertainty Comparison", "MuonSFerr_compare.png")

overlay_hist(el_sf_ev, el_sf_zzt,
             "Electron SF (Events tree)", "Electron SF (ZZTree/candTree)",
             "Electron SF Comparison", "ElectronSF_compare.png")

overlay_hist(el_err_ev, el_err_zzt,
             "Electron SF err (Events tree)", "Electron SF err (ZZTree/candTree)",
             "Electron SF Uncertainty Comparison", "ElectronSFerr_compare.png")

print("Created:")
print("  MuonSF_compare.png")
print("  MuonSFerr_compare.png")
print("  ElectronSF_compare.png")
print("  ElectronSFerr_compare.png")
