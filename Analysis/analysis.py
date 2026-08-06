import MDAnalysis as mda
from MDAnalysis.analysis import rms, rdf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os

os.makedirs("analysis", exist_ok=True)

# Load trajectory
u = mda.Universe("md.tpr", "md_noPBC.xtc")
print(f"Frames: {len(u.trajectory)}")
print(f"Atoms: {len(u.atoms)}")

# Atom selections
ethanol = u.select_atoms("resname ETHH")
ethanol_O = u.select_atoms("resname ETHH and name EO")
water_O = u.select_atoms("resname SOL and name OW")

# RMSD
print("Computing RMSD...")
R = rms.RMSD(ethanol, ethanol, select="resname ETHH", ref_frame=0)
R.run()
time_ps = R.results.rmsd[:, 1]
rmsd_vals = R.results.rmsd[:, 2]

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(time_ps, rmsd_vals, color="#2563EB", linewidth=1.5)
ax.axhline(np.mean(rmsd_vals), color="#DC2626", linestyle="--",
           linewidth=1, label=f"Mean = {np.mean(rmsd_vals):.2f} Å")
ax.set_xlabel("Time (ps)", fontsize=12)
ax.set_ylabel("RMSD (Å)", fontsize=12)
ax.set_title("Ethanol RMSD — 500 ps Production MD\nGROMOS54A7 / SPC Water / 300 K / 1 bar", fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("analysis/rmsd.png", dpi=150)
plt.close()
print(f"RMSD plot saved. Mean RMSD: {np.mean(rmsd_vals):.3f} Å")

# RDF
print("Computing RDF...")
rdf_analysis = rdf.InterRDF(
    ethanol_O, water_O,
    nbins=100,
    range=(0.0, 8.0)
)
rdf_analysis.run()
bins = rdf_analysis.results.bins
g_r = rdf_analysis.results.rdf

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(bins, g_r, color="#16A34A", linewidth=2.0)
ax.fill_between(bins, g_r, alpha=0.15, color="#16A34A")
ax.axhline(1.0, color="gray", linestyle="--", linewidth=1, label="g(r) = 1 (bulk)")
ax.set_xlabel("Distance r (Å)", fontsize=12)
ax.set_ylabel("g(r)", fontsize=12)
ax.set_title("Radial Distribution Function\nEthanol O — Water O | 500 ps | GROMOS54A7", fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 8)
plt.tight_layout()
plt.savefig("analysis/rdf.png", dpi=150)
plt.close()
print("RDF plot saved.")
print("Done. Check analysis/rmsd.png and analysis/rdf.png")
