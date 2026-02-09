import uproot
import awkward as ak

# ------------------------
# ROOT file path
# ------------------------
file_path = "/eos/user/m/mmanoni/HZZ_prod_ggH2022RECO/2022_MC/PROD_samplesNano_2022_MC_cc84ce40/ggH125/ZZ4lAnalysis.root"

# ------------------------
# Open ROOT file
# ------------------------
with uproot.open(file_path) as f:
    tree = f["Events"]

    # Load arrays as awkward arrays
    mu_pt  = tree["Muon_pt"].arrays(library="ak")["Muon_pt"]
    mu_eta = tree["Muon_eta"].arrays(library="ak")["Muon_eta"]
    mu_phi = tree["Muon_phi"].arrays(library="ak")["Muon_phi"]
    mu_sf  = tree["Muon_dataMC"].arrays(library="ak")["Muon_dataMC"]

    el_pt  = tree["Electron_pt"].arrays(library="ak")["Electron_pt"]
    el_eta = tree["Electron_eta"].arrays(library="ak")["Electron_eta"]
    el_phi = tree["Electron_phi"].arrays(library="ak")["Electron_phi"]
    el_sf  = tree["Electron_dataMC"].arrays(library="ak")["Electron_dataMC"]

# ------------------------
# Number of events
# ------------------------
n_events = len(mu_pt)
print(f"Total events in file: {n_events}")

# ------------------------
# Loop over events
# ------------------------
for evt_idx in range(n_events):
    # ------------------------
    # MUONS
    # ------------------------
    for i, (pt, eta, phi, sf) in enumerate(zip(mu_pt[evt_idx], mu_eta[evt_idx], mu_phi[evt_idx], mu_sf[evt_idx])):
        eta_ok = abs(eta) <= 2.5
        pt_ok = 5 <= pt <= 120

        if sf == 1 and eta_ok and pt_ok:
            print(f"\nEvent #{evt_idx}, lepton #{i}")
            print(f"  Lep ID   = 13 (muon)")
            print(f"  Pt       = {pt:.3f} GeV")
            print(f"  Eta      = {eta:.3f}")
            print(f"  Eta      = {eta:.3f}")
            print(f"  Phi      = {phi:.3f}")
            print(f"  SF_lep   = {sf:.6f}")
            print(f"  Full Muon_SF: {mu_sf[evt_idx]}")
            print(f"  Full Electron_SF: {el_sf[evt_idx]}")
            print(f"  eta_ok = {eta_ok}, pt_ok = {pt_ok}")

    # ------------------------
    # ELECTRONS
    # ------------------------
    for i, (pt, eta, phi, sf) in enumerate(zip(el_pt[evt_idx], el_eta[evt_idx], el_phi[evt_idx], el_sf[evt_idx])):
        eta_ok = abs(eta) <= 2.5
        pt_ok = 7 <= pt <= 500

        if sf == 1 and eta_ok and pt_ok:
            print(f"\nEvent #{evt_idx}, lepton #{i}")
            print(f"  Lep ID   = 11 (electron)")
            print(f"  Pt       = {pt:.3f} GeV")
            print(f"  Eta      = {eta:.3f}")
            print(f"  Phi      = {phi:.3f}")
            print(f"  SF_lep   = {sf:.6f}")
            print(f"  Full Muon_SF: {mu_sf[evt_idx]}")
            print(f"  Full Electron_SF: {el_sf[evt_idx]}")
            print(f"  eta_ok = {eta_ok}, pt_ok = {pt_ok}")

