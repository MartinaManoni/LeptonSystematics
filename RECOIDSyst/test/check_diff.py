import ROOT
import uproot
import awkward as ak

ROOT.gROOT.SetBatch(True)

year = "2022"

# ------------------------
# ROOT file configuration
# ------------------------
if year == "2022":
    folder1 = '/eos/user/m/mmanoni/HZZ_prod_ggH2022RECO/2022_MC/PROD_samplesNano_2022_MC_cc84ce40/'
    file_name1 = '/ZZ4lAnalysis_SKIMMED_BestCand_OK.root'
    list1 = ['ggH125']

    file2 = '/eos/user/m/mmanoni/HZZ_prod_ggH2022RECO/2022_MC/PROD_samplesNano_2022_MC_cc84ce40/ggH125/ZZ4lAnalysis.root'
else:
    raise ValueError(f"Unknown year: {year}")

# ------------------------
# Lepton selection function
# ------------------------
def lepton_passes(lep_id, pt, eta, sf):
    eta_ok = abs(eta) <= 2.5
    pt_ok = (5 <= pt <= 120 if lep_id == 13 else 7 <= pt <= 500)
    return eta_ok and pt_ok

# ------------------------
# Load ROOT candTree leptons
# ------------------------
mc = ROOT.TFile.Open(folder1 + list1[0] + file_name1)
tree = mc.Get("ZZTree/candTree")
n_events = tree.GetEntries()
print(f"Total events in ROOT file: {n_events}")

root_leptons = []
for evt_idx in range(n_events):
    tree.GetEntry(evt_idx)
    muon_idx = 0
    electron_idx = 0
    evt_leptons = []

    for i in range(4):
        lep_id = abs(tree.LepLepId[i])
        if lep_id not in [11, 13]:
            continue
        pt = tree.LepPt[i]
        eta = tree.LepEta[i]
        phi = tree.LepPhi[i]
        sf = tree.Muon_SF[muon_idx] if lep_id == 13 else tree.Electron_SF[electron_idx]
        #print("sf", sf)

        evt_leptons.append({"lep_id": lep_id, "pt": pt, "eta": eta, "phi": phi, "sf": sf})

        if lep_id == 13:
            muon_idx += 1
        else:
            electron_idx += 1

    root_leptons.append(evt_leptons)

# ------------------------
# Load uproot leptons + best candidate indices
# ------------------------
with uproot.open(file2) as f:
    tree2 = f["Events"]

    mu_pt  = tree2["Muon_pt"].arrays(library="ak")["Muon_pt"]
    mu_eta = tree2["Muon_eta"].arrays(library="ak")["Muon_eta"]
    mu_phi = tree2["Muon_phi"].arrays(library="ak")["Muon_phi"]
    mu_sf  = tree2["Muon_dataMC"].arrays(library="ak")["Muon_dataMC"]

    el_pt  = tree2["Electron_pt"].arrays(library="ak")["Electron_pt"]
    el_eta = tree2["Electron_eta"].arrays(library="ak")["Electron_eta"]
    el_phi = tree2["Electron_phi"].arrays(library="ak")["Electron_phi"]
    el_sf  = tree2["Electron_dataMC"].arrays(library="ak")["Electron_dataMC"]

    best_idx = tree2["bestCandIdx"].arrays(library="ak")["bestCandIdx"]
    Z1l1_idx = tree2["ZZCand_Z1l1Idx"].arrays(library="ak")["ZZCand_Z1l1Idx"]
    Z1l2_idx = tree2["ZZCand_Z1l2Idx"].arrays(library="ak")["ZZCand_Z1l2Idx"]
    Z2l1_idx = tree2["ZZCand_Z2l1Idx"].arrays(library="ak")["ZZCand_Z2l1Idx"]
    Z2l2_idx = tree2["ZZCand_Z2l2Idx"].arrays(library="ak")["ZZCand_Z2l2Idx"]



# NOTE: NanoAOD merges ELECTRONS FIRST, THEN MUONS (Electron concat Muon)
#       This is why your old code was mismatching.

lep_pt  = ak.concatenate([el_pt,  mu_pt], axis=1)
lep_eta = ak.concatenate([el_eta, mu_eta], axis=1)
lep_phi = ak.concatenate([el_phi, mu_phi], axis=1)
lep_sf  = ak.concatenate([el_sf,  mu_sf], axis=1)

lep_id = ak.concatenate([
    ak.ones_like(el_pt) * 11,
    ak.ones_like(mu_pt) * 13
], axis=1)


# ------------------------
# Compare events
# ------------------------
for evt_idx in range(n_events):

    if best_idx[evt_idx] < 0:
        continue

    tree.GetEntry(evt_idx)  # load ROOT event to access full SF arrays

    root_evt = root_leptons[evt_idx]

    # full ROOT SF arrays for this event
    root_mu_sfs = list(tree.Muon_SF)
    root_el_sfs = list(tree.Electron_SF)
    root_id = list(tree.LepLepId)

    # full NanoAOD SF lists for this event
    nano_mu_sfs = list(mu_sf[evt_idx])
    nano_el_sfs = list(el_sf[evt_idx])

    # ------------------------
    # Extract best-candidate indices
    # ------------------------
    cand = best_idx[evt_idx]

    idxs = [
        Z1l1_idx[evt_idx][cand],
        Z1l2_idx[evt_idx][cand],
        Z2l1_idx[evt_idx][cand],
        Z2l2_idx[evt_idx][cand]
    ]

    # ------------------------
    # Build best-candidate leptons
    # ------------------------
    best_leps = []
    for idx in idxs:
        best_leps.append({
            "lep_id": lep_id[evt_idx][idx],
            "pt":     float(lep_pt[evt_idx][idx]),
            "eta":    float(lep_eta[evt_idx][idx]),
            "phi":    float(lep_phi[evt_idx][idx]),
            "sf":     float(lep_sf[evt_idx][idx]),
        })
    
    '''best_leps = []

    for i, lep_index in enumerate(ZlIdx):
        # skip if index out of bounds
        if lep_index >= len(lep_pt[evt_idx]):
            continue

        lid = abs(tree.LepLepId[i])
        if lid == 11:
            if lep_index >= len(el_sf[evt_idx]):
                continue
            sf = float(el_sf[evt_idx][lep_index])
        elif lid == 13:
            if lep_index >= len(mu_sf[evt_idx]):
                continue
            sf = float(mu_sf[evt_idx][lep_index])
        else:
            continue

        best_leps.append({
            "lep_id": lid,
            "pt":     float(lep_pt[evt_idx][lep_index]),
            "eta":    float(lep_eta[evt_idx][lep_index]),
            "phi":    float(lep_phi[evt_idx][lep_index]),
            "sf":     sf,
        })
        '''
    # ------------------------
    # Compare ROOT vs Uproot leptons
    # and print full SF lists when sf==1
    # ------------------------
    for r_lep, u_lep in zip(root_evt, best_leps):

        if r_lep["sf"] == 1 or u_lep["sf"] == 1:

            print(f"\n=== Event #{evt_idx} ===")

            print("\n--- ROOT candTree full SF lists ---")
            print("ROOT Electron_SF:", root_el_sfs)
            print("ROOT Muon_SF    :", root_mu_sfs)
            print("ROOT Electron_SF:", root_id)

            print("\n--- Uproot NanoAOD full SF lists ---")
            print("NanoAOD Electron_dataMC:", nano_el_sfs)
            print("NanoAOD Muon_dataMC    :", nano_mu_sfs)

            print("\n--- Selected lepton comparison ---")
            print(f"Lepton ID: {r_lep['lep_id']}")
            print(f"  ROOT   -> pt: {r_lep['pt']:.3f}, eta: {r_lep['eta']:.3f}, "
                  f"phi: {r_lep['phi']:.3f}, sf: {r_lep['sf']:.6f}")
            print(f"  Uproot -> pt: {u_lep['pt']:.3f}, eta: {u_lep['eta']:.3f}, "
                  f"phi: {u_lep['phi']:.3f}, sf: {u_lep['sf']:.6f}")
            print(f"  pt_diff  = {abs(r_lep['pt'] - u_lep['pt']):.6f}")
            print(f"  eta_diff = {abs(r_lep['eta'] - u_lep['eta']):.6f}")
            print(f"  phi_diff = {abs(r_lep['phi'] - u_lep['phi']):.6f}")
            print(f"  sf_diff  = {abs(r_lep['sf'] - u_lep['sf']):.6f}")

