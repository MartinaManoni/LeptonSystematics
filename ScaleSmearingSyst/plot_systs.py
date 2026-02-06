import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import math
from scipy.interpolate import interp1d



plt.style.use(hep.style.CMS)

base = "FitResults/"
proc = "DY_2024_mu_"
name = "2024"

if name == "2022":
    lumi_dict = "7.98"
elif name == "2022EE":
    lumi_dict = "26.67"
elif name == "2023preBPix":
    lumi_dict = "18.06"
elif name == "2023postBPix":
    lumi_dict = "9.69"
elif name == "2024":
    lumi_dict = "108.8"

print("lumi_dict", lumi_dict)
# Eta bin edges and centers
eta_edges = np.array([-2.6, -1.95, -1.3, -0.65, 0, 0.65, 1.3, 1.95, 2.6])
eta_centers = 0.5*(eta_edges[:-1] + eta_edges[1:])
Nbins = len(eta_centers)

# Arrays for mu and sigma values and uncertainties
mu_nom, mu_sUp, mu_sDn, mu_mUp, mu_mDn = [np.zeros(Nbins) for _ in range(5)]
err_nom, err_sUp, err_sDn, err_mUp, err_mDn = [np.zeros(Nbins) for _ in range(5)]
sigma_nom, sigma_sUp, sigma_sDn, sigma_mUp, sigma_mDn = [np.zeros(Nbins) for _ in range(5)]
err_sigma_nom, err_sigma_sUp, err_sigma_sDn, err_sigma_mUp, err_sigma_mDn = [np.zeros(Nbins) for _ in range(5)]

# ------------------------------
# Read txt files
# ------------------------------
def read_mu_sigma_error(filename):
    """Reads txt file and returns mu, sigma, mu_err, sigma_err"""
    with open(filename) as f:
        vals = [float(x) for x in f.read().split()]
    mu, sigma = vals[0], vals[1]
    mu_err, sigma_err = vals[2], vals[3]
    return mu, sigma, mu_err, sigma_err

for i in range(1, Nbins+1):
    binName = f"etaBin{i}"

    mu_nom[i-1], sigma_nom[i-1], err_nom[i-1], err_sigma_nom[i-1] = read_mu_sigma_error(f"{base}{proc}{binName}.txt")
    mu_sUp[i-1], sigma_sUp[i-1], err_sUp[i-1], err_sigma_sUp[i-1] = read_mu_sigma_error(f"{base}{proc}scaleUp_{binName}.txt")
    mu_sDn[i-1], sigma_sDn[i-1], err_sDn[i-1], err_sigma_sDn[i-1] = read_mu_sigma_error(f"{base}{proc}scaleDn_{binName}.txt")
    mu_mUp[i-1], sigma_mUp[i-1], err_mUp[i-1], err_sigma_mUp[i-1] = read_mu_sigma_error(f"{base}{proc}smearUp_{binName}.txt")
    mu_mDn[i-1], sigma_mDn[i-1], err_mDn[i-1], err_sigma_mDn[i-1] = read_mu_sigma_error(f"{base}{proc}smearDn_{binName}.txt")

# ------------------------------
# Safe ratio calculation
# ------------------------------
def safe_ratio(nom, nom_err, var, var_err):
    ratio = var / nom
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio_err = ratio * np.sqrt((var_err/var)**2 + (nom_err/nom)**2)
        ratio_err = np.nan_to_num(ratio_err, nan=0.0)
    return ratio, ratio_err

# ------------------------------
# Generic plot function
# ------------------------------
def plot_quantity_with_ratio(eta_centers, nom, err_nom, up, err_up, down, err_down, lumi_dict,
                             ylabel, title, fname, var_type="Scale"):

    ratio_up, ratio_err_up = safe_ratio(nom, err_nom, up, err_up)
    ratio_dn, ratio_err_dn = safe_ratio(nom, err_nom, down, err_down)

    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(10,8), sharex=True,
                                   gridspec_kw={'height_ratios':[3,1]})
    
    # CMS label
    #hep.cms.label(ax=ax1, data=False)
    hep.cms.label(ax=ax1, data=True, label="Preliminary", lumi=f"{lumi_dict}")
    #ax1.text(0.65, 0.85, r"26.7 fb$^{-1}$ (13.6 TeV)", transform=ax1.transAxes,
             #fontsize=16, horizontalalignment='left')

    # Top pad
    ax1.errorbar(eta_centers, nom, yerr=err_nom, fmt='o', color='k', label='Nominal', capsize=3)
    ax1.errorbar(eta_centers, up, yerr=err_up, fmt='o', color='b', label=f'{var_type} Up', capsize=3)
    ax1.errorbar(eta_centers, down, yerr=err_down, fmt='o', color='r', label=f'{var_type} Down', capsize=3)
    ax1.set_ylabel(ylabel, fontsize=25)
    ax1.set_title(title, fontsize=16)
    ax1.grid(True)
    ax1.legend()

    # Bottom pad: ratio
    ax2.axhline(1.0, color='k', linestyle='--', linewidth=1)
    ax2.errorbar(eta_centers, ratio_up, yerr=ratio_err_up, fmt='o', color='b', capsize=3)
    ax2.errorbar(eta_centers, ratio_dn, yerr=ratio_err_dn, fmt='o', color='r', capsize=3)

    ax2.set_xlabel(r"$\eta$", fontsize=25)
    ax2.set_ylabel("Var./Nom.", fontsize=25)
    # dynamic y-limits
    ratio_min = min(ratio_up.min()-ratio_err_up.max(), ratio_dn.min()-ratio_err_dn.max())
    ratio_max = max(ratio_up.max()+ratio_err_up.max(), ratio_dn.max()+ratio_err_dn.max())
    margin = 0.01
    ax2.set_ylim(ratio_min-margin, ratio_max+margin)
    ax2.grid(True)
    ax2.set_xticks(eta_centers)
    ax2.set_xticklabels([f"{x:.2f}" for x in eta_centers])

    plt.tight_layout()
    plt.savefig(fname+".png")
    plt.savefig(fname+".pdf")
    plt.close()


def final_average_ratio(nom, up, down):
    ratios_up = up / nom
    ratios_dn = down / nom

    A_values = []

    # Loop over ALL eta bins
    for i in range(len(nom)):
        Ai = 0.5 * (np.abs(ratios_up[i]-1) + np.abs(ratios_dn[i]-1))
        A_values.append(Ai)
    # Final average across ALL eta bins
    final_value = np.mean(A_values)
    percentage_value= abs((final_value)*100)
    return np.array(A_values), final_value, percentage_value


def final_average_ratio_envelope(nom, up, down):
    ratio_up = up / nom
    ratio_dn = down / nom

    A_values = []
    for i in range(len(nom)):
        Ai = max(abs(ratio_up[i] - 1), abs(ratio_dn[i] - 1))
        A_values.append(Ai)

    final_value = np.mean(A_values)
    percentage_value = abs(final_value * 100)

    return np.array(A_values), final_value, percentage_value


def cms_envelope_systematic(nom, up, down):
    ratio_up = up / nom
    ratio_dn = down / nom

    # per-bin envelope deviations
    delta = np.maximum(np.abs(ratio_up - 1), np.abs(ratio_dn - 1))

    max_dev = np.max(delta)

    percentage = max_dev * 100.0
    return delta, max_dev, percentage

# ------------------------------
# Plot all 8 plots
# ------------------------------
plot_quantity_with_ratio(eta_centers, mu_nom, err_nom, mu_sUp, err_sUp, mu_sDn, err_sDn, lumi_dict,
                         r"$Z_{mm}\ \mu$ [GeV]", "", f"MuVsEta_Scale_{name}_mu", var_type="Scale")

plot_quantity_with_ratio(eta_centers, sigma_nom, err_sigma_nom, sigma_mUp, err_sigma_mUp, sigma_mDn, err_sigma_mDn, lumi_dict,
                         r"$Z_{mm}\ \sigma$ [GeV]", "", f"SigmaVsEta_Smear_{name}_mu", var_type="Smear")

print("DONE! All 8 plots saved:")

# Results Mean over eta bins
A_mu_scale, F_mu_scale, P_mu_scale= final_average_ratio(mu_nom, mu_sUp, mu_sDn)
A_sigma_smear, F_sigma_smear, P_sigma_smear = final_average_ratio(sigma_nom, sigma_mUp, sigma_mDn)
print("Run-III scale sys (Mean)=", P_mu_scale, "%")
print("Run-III smear sys (Mean)=", P_sigma_smear, "%")


# Results Envelope (take the max variation)
delta_scale, max_scale, pct_scale = cms_envelope_systematic(mu_nom, mu_sUp, mu_sDn)
delta_smear, max_smear, pct_smear = cms_envelope_systematic(sigma_nom, sigma_mUp, sigma_mDn)

print("Run-III scale sys (Envelope)=", pct_scale, "%")
print("Run-III smear sys (Envelope)=", pct_smear, "%")

