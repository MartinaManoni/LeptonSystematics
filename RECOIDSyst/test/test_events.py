import ROOT

ROOT.gROOT.SetBatch(True)

year = "2022"

if year == "2022":
    folder ='/eos/user/m/mmanoni/HZZ_prod_ggH2022RECO/2022_MC/PROD_samplesNano_2022_MC_cc84ce40/'
    file_name = '/ZZ4lAnalysis_SKIMMED.root'
    List = ['ggH125']

elif year == "2024":
    folder = '/eos/user/m/mmanoni/prod_ggH/2022/'
    file_name = '/2024_skimmed.root'
    List = ['ggH125']

else:
    raise ValueError(f"Unknown year: {year}")

for Type in List:
    mc = ROOT.TFile.Open(folder + Type + file_name)
    tree = mc.Get("ZZTree/candTree")

    n_events = tree.GetEntries()
    print(f"Total events in file: {n_events}")

    # ----------------------------
    # Loop over all events
    # ----------------------------
    for evt_idx in range(n_events):
        tree.GetEntry(evt_idx)

        muon_idx = 0
        electron_idx = 0

        for i in range(4):
            lep_id = abs(tree.LepLepId[i])

            # Compute eta/pt selection
            eta_ok = abs(tree.LepEta[i]) <= 2.5
            pt_ok = (
                (lep_id == 13 and 5 <= tree.LepPt[i] <= 120) or
                (lep_id == 11 and 7 <= tree.LepPt[i] <= 500)
            )

            # Retrieve SF for this lepton
            if lep_id == 13:
                SF_lep = tree.Muon_SF[muon_idx]
            elif lep_id == 11:
                SF_lep = tree.Electron_SF[electron_idx]
            else:
                continue  # skip unexpected lepton IDs

            # Print info only if SF == 1 AND lepton passes selection
            if SF_lep == 1 and (eta_ok or pt_ok):
                print(f"\nEvent #{evt_idx}, lepton #{i}")
                print(f"  Lep ID   = {lep_id} ({'muon' if lep_id==13 else 'electron'})")
                print(f"  Pt       = {tree.LepPt[i]:.3f} GeV")
                print(f"  Eta      = {tree.LepEta[i]:.3f}")
                print(f"  Phi      = {tree.LepPhi[i]:.3f}")
                print(f"  SF_lep   = {SF_lep:.6f}")
                print(f"  Full Muon_SF: {getattr(tree, 'Muon_SF', 'N/A')}")
                print(f"  Full Electron_SF: {getattr(tree, 'Electron_SF', 'N/A')}")
                print(f"  eta_ok = {eta_ok}, pt_ok = {pt_ok}")

            # Increment lepton-specific index
            if lep_id == 13:
                muon_idx += 1
            elif lep_id == 11:
                electron_idx += 1
