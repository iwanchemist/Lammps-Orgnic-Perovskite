#!/usr/bin/env python3

import numpy as np
from pathlib import Path

# ============================================================
# SETTINGS
# ============================================================

RDF_DIR = Path("rdf_results")

# Number of atoms of the second species
N_PB = 32
N_I  = 96

# Simulation box
LX = 25.2
LY = 25.2
LZ = 32.0

VOLUME = LX * LY * LZ

# RDF pairs
pairs = {
    "C-Pb": N_PB,
    "C-I":  N_I,
    "N-Pb": N_PB,
    "N-I":  N_I,
    "S-Pb": N_PB,
    "S-I":  N_I,
}


# ============================================================
# FIND FIRST PEAK AND FIRST MINIMUM
# ============================================================

def find_first_shell(r, g):

    # Ignore unrealistically short distances
    mask = r >= 0.5

    r2 = r[mask]
    g2 = g[mask]

    # First major peak
    peak_index = np.argmax(g2)

    peak_r = r2[peak_index]
    peak_g = g2[peak_index]

    # Search for the first minimum after the peak
    indices = np.arange(
        peak_index + 1,
        len(g2)
    )

    # Search within 3 Å after the peak
    indices = indices[
        r2[indices] <= peak_r + 3.0
    ]

    if len(indices) == 0:
        return peak_r, peak_g, np.nan

    minimum_index = indices[
        np.argmin(g2[indices])
    ]

    cutoff = r2[minimum_index]

    return peak_r, peak_g, cutoff


# ============================================================
# COORDINATION NUMBER
# ============================================================

def calculate_cn(r, g, rho, cutoff):

    if np.isnan(cutoff):
        return np.nan

    dr = r[1] - r[0]

    mask = r <= cutoff

    cn = (
        4.0
        * np.pi
        * rho
        * np.sum(
            g[mask]
            * r[mask]**2
            * dr
        )
    )

    return cn


# ============================================================
# MAIN
# ============================================================

print("=" * 80)
print("COORDINATION NUMBER ANALYSIS")
print("Phenothiazine - CsPbI3")
print("=" * 80)

print(f"Simulation volume : {VOLUME:.3f} Å^3")
print()

results = []

for pair, n_b in pairs.items():

    filename = RDF_DIR / f"RDF_{pair.replace('-', '_')}.dat"

    if not filename.exists():

        print(f"WARNING: {filename} not found")
        continue

    data = np.loadtxt(
        filename,
        comments="#"
    )

    r = data[:, 0]
    g = data[:, 1]

    # Number density of species B
    rho_b = n_b / VOLUME

    # First peak and first minimum
    peak_r, peak_g, cutoff = find_first_shell(
        r,
        g
    )

    # Coordination number
    cn = calculate_cn(
        r,
        g,
        rho_b,
        cutoff
    )

    results.append([
        pair,
        peak_r,
        peak_g,
        cutoff,
        cn
    ])


# ============================================================
# SORT RESULTS
# ============================================================

# Desired order
order = [
    "C-Pb",
    "C-I",
    "N-Pb",
    "N-I",
    "S-Pb",
    "S-I"
]

results.sort(
    key=lambda x: order.index(x[0])
)


# ============================================================
# SAVE ALL RESULTS IN ONE DAT FILE
# ============================================================

dat_file = RDF_DIR / "coordination_numbers.dat"

with open(dat_file, "w") as f:

    f.write(
        "# Coordination number analysis: "
        "Phenothiazine - CsPbI3\n"
    )

    f.write(
        f"# Simulation box volume = "
        f"{VOLUME:.5f} Angstrom^3\n"
    )

    f.write(
        "#\n"
    )

    f.write(
        "# Pair  "
        "FirstPeak_r(Angstrom)  "
        "Peak_g(r)  "
        "FirstMinimum_r(Angstrom)  "
        "CoordinationNumber\n"
    )

    for row in results:

        pair, peak_r, peak_g, cutoff, cn = row

        f.write(
            f"{pair:5s} "
            f"{peak_r:20.6f} "
            f"{peak_g:12.6f} "
            f"{cutoff:25.6f} "
            f"{cn:18.6f}\n"
        )


# ============================================================
# SAVE ALL RESULTS IN ONE CSV FILE
# ============================================================

csv_file = RDF_DIR / "coordination_numbers.csv"

with open(csv_file, "w") as f:

    f.write(
        "Pair,"
        "FirstPeak_r_A,"
        "Peak_g_r,"
        "FirstMinimum_r_A,"
        "CoordinationNumber\n"
    )

    for row in results:

        pair, peak_r, peak_g, cutoff, cn = row

        f.write(
            f"{pair},"
            f"{peak_r:.6f},"
            f"{peak_g:.6f},"
            f"{cutoff:.6f},"
            f"{cn:.6f}\n"
        )


# ============================================================
# PRINT SUMMARY
# ============================================================

print(
    f"{'Pair':<8}"
    f"{'Peak r (Å)':>15}"
    f"{'Peak g(r)':>15}"
    f"{'Minimum (Å)':>18}"
    f"{'CN':>15}"
)

print("-" * 80)

for row in results:

    pair, peak_r, peak_g, cutoff, cn = row

    print(
        f"{pair:<8}"
        f"{peak_r:>15.3f}"
        f"{peak_g:>15.3f}"
        f"{cutoff:>18.3f}"
        f"{cn:>15.3f}"
    )

print()
print("=" * 80)
print("RESULT FILES")
print("=" * 80)
print(f"DAT : {dat_file}")
print(f"CSV : {csv_file}")
print("=" * 80)
