#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# USER SETTINGS
# ============================================================

TRAJ = "trajectory.lammpstrj"

# Atom types from system.data
TYPE_PB = 2
TYPE_I  = 3

TYPE_C  = 4
TYPE_N  = 5
TYPE_S  = 6

# RDF parameters
R_MAX = 12.0       # Angstrom
N_BINS = 240

# Analyze every frame
FRAME_STRIDE = 1

# First-shell cutoff used for coordination number.
# This is only a starting value; inspect RDF minima and adjust.
CN_CUTOFF = 4.0

# ============================================================
# LAMMPS TRAJECTORY READER
# ============================================================

def read_lammps_dump(filename):
    """
    Read a LAMMPS custom dump with:

    ITEM: TIMESTEP
    ITEM: NUMBER OF ATOMS
    ITEM: BOX BOUNDS
    ITEM: ATOMS id mol type q x y z vx vy vz
    """

    frames = []

    with open(filename, "r") as f:

        while True:

            line = f.readline()

            if not line:
                break

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(f.readline().strip())

            # NUMBER OF ATOMS
            line = f.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise RuntimeError("Unexpected trajectory format.")

            natoms = int(f.readline().strip())

            # BOX BOUNDS
            line = f.readline()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError("BOX BOUNDS not found.")

            bounds = []

            for _ in range(3):
                values = list(map(float, f.readline().split()))
                bounds.append(values[:2])

            bounds = np.array(bounds)

            # ATOMS
            line = f.readline()

            if not line.startswith("ITEM: ATOMS"):
                raise RuntimeError("ATOMS section not found.")

            columns = line.strip().split()[2:]

            data = []

            for _ in range(natoms):
                data.append(list(map(float, f.readline().split())))

            data = np.array(data)

            frames.append({
                "timestep": timestep,
                "natoms": natoms,
                "bounds": bounds,
                "columns": columns,
                "data": data
            })

    return frames


# ============================================================
# EXTRACT ATOM POSITIONS
# ============================================================

def get_positions(frame):

    data = frame["data"]
    columns = frame["columns"]

    type_idx = columns.index("type")
    x_idx = columns.index("x")
    y_idx = columns.index("y")
    z_idx = columns.index("z")

    atom_types = data[:, type_idx].astype(int)

    positions = data[:, [x_idx, y_idx, z_idx]]

    return atom_types, positions


# ============================================================
# MINIMUM IMAGE DISTANCE
# ============================================================

def minimum_image(delta, box):

    box_lengths = box[:, 1] - box[:, 0]

    delta -= box_lengths * np.round(delta / box_lengths)

    return delta


# ============================================================
# RDF CALCULATION
# ============================================================

def calculate_rdf(frames, type_a, type_b):

    edges = np.linspace(0.0, R_MAX, N_BINS + 1)

    r = 0.5 * (edges[:-1] + edges[1:])

    hist = np.zeros(N_BINS)

    total_pairs = 0

    nframes_used = 0

    for iframe, frame in enumerate(frames):

        if iframe % FRAME_STRIDE != 0:
            continue

        atom_types, positions = get_positions(frame)

        pos_a = positions[atom_types == type_a]
        pos_b = positions[atom_types == type_b]

        if len(pos_a) == 0 or len(pos_b) == 0:
            continue

        box = frame["bounds"]

        box_lengths = box[:, 1] - box[:, 0]

        # Calculate all A-B distances
        for pa in pos_a:

            delta = pos_b - pa

            delta = minimum_image(delta, box)

            distances = np.sqrt(np.sum(delta**2, axis=1))

            hist += np.histogram(
                distances,
                bins=edges
            )[0]

            total_pairs += len(distances)

        nframes_used += 1

    # --------------------------------------------------------
    # Normalize RDF
    #
    # g_AB(r) = counts /
    #           [N_A * rho_B * shell_volume * N_frames]
    # --------------------------------------------------------

    if nframes_used == 0:
        raise RuntimeError("No frames available.")

    # Number of A atoms per frame
    atom_types, _ = get_positions(frames[0])
    n_a = np.sum(atom_types == type_a)

    # Number of B atoms per frame
    n_b = np.sum(atom_types == type_b)

    box = frames[0]["bounds"]
    box_lengths = box[:, 1] - box[:, 0]
    volume = np.prod(box_lengths)

    rho_b = n_b / volume

    shell_volume = (
        4.0 / 3.0
        * np.pi
        * (edges[1:]**3 - edges[:-1]**3)
    )

    normalization = (
        nframes_used
        * n_a
        * rho_b
        * shell_volume
    )

    g_r = hist / normalization

    return r, g_r


# ============================================================
# COORDINATION NUMBER
# ============================================================

def coordination_number(r, g_r, n_b, volume, cutoff):

    rho_b = n_b / volume

    mask = r <= cutoff

    dr = r[1] - r[0]

    cn = 4.0 * np.pi * rho_b * np.sum(
        g_r[mask] * r[mask]**2 * dr
    )

    return cn


# ============================================================
# MAIN
# ============================================================

print("=" * 70)
print("Phenothiazine / CsPbI3 RDF Analysis")
print("=" * 70)

print(f"Trajectory : {TRAJ}")
print(f"Rmax       : {R_MAX:.2f} Å")
print(f"Bins       : {N_BINS}")
print()

frames = read_lammps_dump(TRAJ)

print(f"Frames read: {len(frames)}")

if len(frames) == 0:
    raise RuntimeError("No trajectory frames found.")

print(f"Atoms/frame: {frames[0]['natoms']}")

# ------------------------------------------------------------
# Atom mapping
# ------------------------------------------------------------

pairs = {
    "C-Pb": (TYPE_C, TYPE_PB),
    "C-I":  (TYPE_C, TYPE_I),

    "N-Pb": (TYPE_N, TYPE_PB),
    "N-I":  (TYPE_N, TYPE_I),

    "S-Pb": (TYPE_S, TYPE_PB),
    "S-I":  (TYPE_S, TYPE_I),
}

# ------------------------------------------------------------
# Create output directory
# ------------------------------------------------------------

outdir = Path("rdf_results")
outdir.mkdir(exist_ok=True)

results = {}

# ------------------------------------------------------------
# Calculate RDF
# ------------------------------------------------------------

for name, (type_a, type_b) in pairs.items():

    print(f"Calculating {name} ...")

    r, g = calculate_rdf(
        frames,
        type_a,
        type_b
    )

    results[name] = (r, g)

    # Save numerical RDF
    outfile = outdir / f"RDF_{name.replace('-', '_')}.dat"

    np.savetxt(
        outfile,
        np.column_stack([r, g]),
        header=f"r(Angstrom)   g(r)    {name}"
    )

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    plt.figure(figsize=(7, 5))

    plt.plot(
        r,
        g,
        linewidth=1.5
    )

    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")

    plt.title(f"{name} RDF")

    plt.xlim(0, R_MAX)

    plt.grid(alpha=0.25)

    plt.tight_layout()

    pngfile = outdir / f"RDF_{name.replace('-', '_')}.png"

    plt.savefig(
        pngfile,
        dpi=300
    )

    plt.close()


# ============================================================
# COMPARISON PLOTS
# ============================================================

# C vs Pb/I
plt.figure(figsize=(7, 5))

for name in ["C-Pb", "C-I"]:
    r, g = results[name]
    plt.plot(r, g, label=name, linewidth=1.5)

plt.xlabel("r (Å)")
plt.ylabel("g(r)")
plt.title("Phenothiazine C vs CsPbI3")
plt.xlim(0, R_MAX)
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(outdir / "RDF_C_comparison.png", dpi=300)
plt.close()


# N vs Pb/I
plt.figure(figsize=(7, 5))

for name in ["N-Pb", "N-I"]:
    r, g = results[name]
    plt.plot(r, g, label=name, linewidth=1.5)

plt.xlabel("r (Å)")
plt.ylabel("g(r)")
plt.title("Phenothiazine N vs CsPbI3")
plt.xlim(0, R_MAX)
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(outdir / "RDF_N_comparison.png", dpi=300)
plt.close()


# S vs Pb/I
plt.figure(figsize=(7, 5))

for name in ["S-Pb", "S-I"]:
    r, g = results[name]
    plt.plot(r, g, label=name, linewidth=1.5)

plt.xlabel("r (Å)")
plt.ylabel("g(r)")
plt.title("Phenothiazine S vs CsPbI3")
plt.xlim(0, R_MAX)
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(outdir / "RDF_S_comparison.png", dpi=300)
plt.close()


# ============================================================
# COMBINED PLOT
# ============================================================

plt.figure(figsize=(8, 6))

for name, (r, g) in results.items():

    plt.plot(
        r,
        g,
        label=name,
        linewidth=1.3
    )

plt.xlabel("r (Å)")
plt.ylabel("g(r)")
plt.title("Phenothiazine–CsPbI3 Partial RDFs")

plt.xlim(0, R_MAX)

plt.legend(
    ncol=2,
    fontsize=9
)

plt.grid(alpha=0.25)

plt.tight_layout()

plt.savefig(
    outdir / "RDF_all.png",
    dpi=300
)

plt.close()


# ============================================================
# PRINT FIRST PEAK
# ============================================================

print()
print("=" * 70)
print("FIRST RDF PEAKS")
print("=" * 70)

for name, (r, g) in results.items():

    # Ignore first 0.5 Å
    mask = r >= 0.5

    peak_index = np.argmax(
        g[mask]
    )

    r_valid = r[mask]
    g_valid = g[mask]

    peak_r = r_valid[peak_index]
    peak_g = g_valid[peak_index]

    print(
        f"{name:5s} : "
        f"peak r = {peak_r:6.3f} Å, "
        f"g(r) = {peak_g:10.4f}"
    )


print()
print("=" * 70)
print("Output written to:", outdir)
print("=" * 70)
