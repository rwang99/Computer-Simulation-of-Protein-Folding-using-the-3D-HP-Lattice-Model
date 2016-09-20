# Computer-Simulation-of-Protein-Folding-using-the-3D-HP-Lattice-Model
- All code written by Competition Entrant.
- Outputs a data.txt general data file providing information on energies in both forms, and trial number
- Outputs coordinate output_%i.xyz file for each trial that can be loaded with programs such as Visual Molecular Dynamics to view 3D fold conformation and path.

Main parameters:
- protein: HP sequence, capitalized
- numsteps: Number of evaluating steps per trial
- temperature: Starting temperature of simulation
- rate: simulated annealing rate, inversely proportional to numsteps
- numtrials: number of repetitive trials
- energytype: energy function type. -1 is distance-based energy function, 1 is contact-based energy function
