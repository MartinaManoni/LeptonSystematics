import ROOT

ROOT.gROOT.SetBatch(True)

year = "2022"

if year == "2022":
    folder ='/eos/user/m/mmanoni/HZZ_prod_ggH2022RECO/2022_MC/PROD_samplesNano_2022_MC_cc84ce40/'
    file_name = '/ZZ4lAnalysis_SKIMMED_BestCand.root'
    List = ['ggH125']

elif year == "2024":
    folder = '/eos/user/m/mmanoni/prod_ggH/2022/'
    file_name = '/2024_skimmed.root'
    List = ['ggH125']

else:
    raise ValueError(f"Unknown year: {year}")


event_to_print =  1369  # <<< CHANGE EVENT NUMBER HERE

for Type in List:
    mc = ROOT.TFile.Open(folder + Type + file_name)
    tree = mc.Get("ZZTree/candTree")

    # Check tree size
    n_events = tree.GetEntries()
    print(f"Total events in file: {n_events}")

    if event_to_print >= n_events:
        print(f"Requested event {event_to_print} exceeds total events!")
        continue

    # Get the specific event
    tree.GetEntry(event_to_print)

    print(f"\n========== Event {event_to_print} ==========")

    muon_idx = 0
    electron_idx = 0

    for i in range(4):
        lep_id = abs(tree.LepLepId[i])

        # -------- MUONS --------
        if lep_id == 13:
            sf = tree.Muon_SF[muon_idx]
            print(f"Muon {muon_idx}: SF = {sf}, eta = {tree.LepEta[i]}, phi = {tree.LepPhi[i]}, pt = {tree.LepPt[i]}")
            muon_idx += 1

        # -------- ELECTRONS --------
        elif lep_id == 11:
            sf = tree.Electron_SF[electron_idx]
            print(f"Electron {electron_idx}: SF = {sf}, eta = {tree.LepEta[i]}, phi = {tree.LepPhi[i]}, pt = {tree.LepPt[i]}")
            electron_idx += 1




