# mean-field-dipolar
Simulate mean field dipolar array properties

Use to calculate phase diagram lobes for an array of dipolar bosons in an optical lattice.
The simulation is based on C Trefzger et al J. Phys. B: At. Mol. Opt. Phys. 44 193001 (2011).

# Installation

1. Install Rust
2. Clone repository
2. Build using `cargo build --release` in the repository directory
3. Run by typing `./target/release/mean-field-dipolar -r <RESULTS PATH>`

# Configuration

Configuration file which defines simulation sweep parameters is sim.toml in the main repository directory.
The parameters that are swept are interaction range (from range_start to range_end in steps of 1) and system size (from size_start to size_end in steps of 2).
Other initialized parameters are theta and phi, defining the dipole angle and u_onsite, defining the onsite interaction energy U.
Currently, the script supports patterns "filled" (fully filled lattice), "cb" (checkerboard lattice), and "hstripe" (horizontal stripe). 
