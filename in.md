# CsPbI3 slab + phenothiazine (PTZ) rigid-body MD
# real units, PBC x/y/z; z contains vacuum
units real
atom_style full
boundary p p p
read_data system.data

# --- Force fields ---
# One global long-range Coulomb term + Buckingham for CsPbI3 + LJ for PTZ.
# This avoids assigning multiple long-range substyles in a hybrid pair style.
pair_style hybrid/overlay coul/long 12.0 buck 12.0 lj/cut 12.0
kspace_style pppm 1.0e-5

# Long-range Coulomb applies to every atom pair.
pair_coeff * * coul/long

# CsPbI3 Buckingham parameters from the published CsPbI3 BC model.
# A and C converted from eV to kcal/mol; rho remains Angstrom.
pair_coeff 1 1 buck 8.2649004e+17 0.0843 5542.3721
pair_coeff 1 2 buck 7.625662e+12 0.10519 11069.063
pair_coeff 1 3 buck 113296.47 0.3814 11070.447
pair_coeff 2 2 buck 70360511 0.131258 0
pair_coeff 2 3 buck 103496.89 0.321737 0
pair_coeff 3 3 buck 22793.507 0.482217 696.95548

# PTZ self and PTZ-inorganic LJ parameters.
# PTZ self parameters are OPLS-like starting values.
# Inorganic LJ values/cross terms are screening placeholders.
pair_coeff 4 4 lj/cut 0.07000000 3.55000000
pair_coeff 5 5 lj/cut 0.17000000 3.25000000
pair_coeff 6 6 lj/cut 0.25000000 3.56000000
pair_coeff 7 7 lj/cut 0.03000000 2.42000000
pair_coeff 1 4 lj/cut 0.05916080 3.42500000
pair_coeff 1 5 lj/cut 0.09219544 3.27500000
pair_coeff 1 6 lj/cut 0.11180340 3.43000000
pair_coeff 1 7 lj/cut 0.03872983 2.86000000
pair_coeff 2 4 lj/cut 0.08366600 3.52500000
pair_coeff 2 5 lj/cut 0.13038405 3.37500000
pair_coeff 2 6 lj/cut 0.15811388 3.53000000
pair_coeff 2 7 lj/cut 0.05477226 2.96000000
pair_coeff 3 4 lj/cut 0.11832160 3.77500000
pair_coeff 3 5 lj/cut 0.18439089 3.62500000
pair_coeff 3 6 lj/cut 0.22360680 3.78000000
pair_coeff 3 7 lj/cut 0.07745967 3.21000000
# No explicit bonded terms: PTZ is treated as a rigid molecule.
group perovskite type 1 2 3
group ptz molecule 1

# Exclude intramolecular PTZ nonbonded interactions because rigid-body dynamics
# preserves the molecular geometry.
neigh_modify exclude group ptz ptz

# Freeze nothing: surface and PTZ are allowed to relax.
compute eint ptz group/group perovskite pair yes kspace yes
thermo 1000
thermo_style custom step temp pe ke etotal press vol c_eint

# Initial minimization: keep the rigid PTZ geometry fixed while the perovskite surface relaxes.
fix hold_ptz ptz setforce 0.0 0.0 0.0
min_style cg
minimize 1.0e-6 1.0e-8 5000 20000
unfix hold_ptz

# Save minimized configuration
write_data minimized.data

# Gentle thermalization at 300 K
velocity ptz set 0.0 0.0 0.0
velocity perovskite create 50.0 4928459 mom yes rot yes dist gaussian

timestep 0.25
fix int1 perovskite nvt temp 50.0 300.0 100.0
fix int2 ptz rigid/nvt molecule temp 50.0 300.0 100.0
run 40000
unfix int1
unfix int2

# Production NVT at 300 K
fix prod1 perovskite nvt temp 300.0 300.0 100.0
fix prod2 ptz rigid/nvt molecule temp 300.0 300.0 100.0

variable Eint equal c_eint
fix energyout all ave/time 100 10 1000 v_Eint file interaction_energy.dat
compute mdTemp all temp
dump mdtraj all custom 1000 trajectory.lammpstrj id mol type q x y z vx vy vz
fix tempout all ave/time 100 10 1000 c_mdTemp file temperature.dat


run 400000

write_data final.data
