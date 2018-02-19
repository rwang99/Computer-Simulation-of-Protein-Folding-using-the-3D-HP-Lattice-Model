# Computer-Simulation-of-Protein-Folding-using-the-3D-HP-Lattice-Model

<h4> Author's note </h4>
<p> The code for this project was written in 2016, while I was still in high school and programmed without much knowledge of best coding practices and efficiency. As such, this code does not fully represent my current coding capabilities, but is still a big part of my journey to learn computer science so feel free to look and play around with it! Enjoy!</p>
<p> To learn more about the sciency aspects of my project, feel free to look at the poster or paper included in this repo! </p>


<h4> Instructions </h4>
<p> This program attempts to find the most stable fold for a protein chain sequence using the 3D Hydrophobic-Polar lattice model. To test different chains, change the protein variable at the top. Only modify the variables at the top to configure your run environment. Can run trials based on a set period of time (ex overnight for 10 hours), or run a certain number of trials regardless of time (ex 15 trials). Change the energy type variable to choose between either a contact-based energy function or distance-based energy function. </p>

<p> The program outputs a data.txt general data file providing information on energies in both forms, and trial number. Individual data/3D folds are outputed to output_(TRIALNUMBER).xyz for each trial that can be visualized with programs such as Visual Molecular Dynamics. Use this .xyz file to view progression of folds and for cool 3D graphics. </p>

![Image of sample fold](https://i.imgur.com/XoR82GG.png)

![Image of fold progression](https://i.imgur.com/imEgAJt.png)

Images taken from my poster! ^


Main parameters:
- protein: HP sequence, capitalized
- numsteps: Number of evaluating steps per trial
- temperature: Starting temperature of simulation
- rate: simulated annealing rate, inversely proportional to numsteps
- numtrials: number of repetitive trials
- energytype: energy function type. -1 is distance-based energy function, 1 is contact-based energy function

