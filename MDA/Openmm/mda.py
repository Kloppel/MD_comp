import MDAnalysis as mda

u = mda.Universe("1eru.psf", "1eru_eq.dcd")

u.trajectory[-1]

u.atoms.write("1eru_restart_eq")