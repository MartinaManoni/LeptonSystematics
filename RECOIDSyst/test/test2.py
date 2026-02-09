import uproot
import awkward as ak
import matplotlib.pyplot as plt

# ROOT file path
#file_path = '/eos/user/m/mmanoni/prod_ggH/2022/PROD_samplesNano_2022_MC_cc84ce40/ggH125/ZZ4lAnalysis.root'

#file_path = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/ggH125/ZZ4lAnalysis.root"
file_path ="/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2023preBPix_MC/ggH125/ZZ4lAnalysis.root"

# Open the ROOT file and read the Events tree
with uproot.open(file_path) as f:
    tree = f["Events"]

    # Load arrays as awkward arrays
    mu_pt = tree["Muon_pt"].arrays(library="ak")["Muon_pt"]
    mu_eta = tree["Muon_eta"].arrays(library="ak")["Muon_eta"]
    mu_sf = tree["Muon_dataMC"].arrays(library="ak")["Muon_dataMC"]

    el_pt = tree["Electron_pt"].arrays(library="ak")["Electron_pt"]
    el_eta = tree["Electron_eta"].arrays(library="ak")["Electron_eta"]
    el_sf = tree["Electron_dataMC"].arrays(library="ak")["Electron_dataMC"]

# ----- Mask SF=0 -----
mask_mu = mu_sf == 1
mask_el = el_sf == 1

# Flatten masked arrays using awkward
pt_zero_mu = ak.flatten(mu_pt[mask_mu])
eta_zero_mu = ak.flatten(mu_eta[mask_mu])

pt_zero_el = ak.flatten(el_pt[mask_el])
eta_zero_el = ak.flatten(el_eta[mask_el])

# ----- Scatter plots -----
if len(pt_zero_mu) > 0:
    plt.figure(figsize=(8,6))
    plt.scatter(eta_zero_mu, pt_zero_mu, s=15, alpha=0.7, color='blue', label='Muon SF=0')
    plt.xlabel("Muon η")
    plt.ylabel("Muon pT [GeV]")
    plt.title("Muons with SF=0")
    plt.ylim(0, 700)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("scatter_uproot_mu.png")
else:
    print("No muons with SF=0 found.")

if len(pt_zero_el) > 0:
    plt.figure(figsize=(8,6))
    plt.scatter(eta_zero_el, pt_zero_el, s=15, alpha=0.7, color='red', label='Electron SF=0')
    plt.xlabel("Electron η")
    plt.ylabel("Electron pT [GeV]")
    plt.title("Electrons with SF=0")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("scatter_uproot_el.png")
else:
    print("No electrons with SF=0 found.")


