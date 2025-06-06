# Computational Biology Homework 8

## Q1

Docking is a method to predict orientation of one molecule when it binds to a second molecule to form a stable complex. It is often used to model the the interaction between a small molecule and a larger macromolecule, typically a protein.

## Q2

- Docking predicts the most likely three-dimensional arrangement of the ligand within the receptor's binding site. This includes the position and conformation of the ligand and how it fits into the binding pocket.

- Docking uses scoring functions to estimate the strength of the interaction between the ligand and the receptor in the predicted binding mode (binding energy).

- Docking can help identify the specific atomic interactions that stabilize the ligand-receptor complex, such as hydrogen bonds.

- By visualizing how a molecule binds, docking can help formulate hypotheses about how it might activate or inhibit a protein's function, or how mutations in the protein might affect binding.

## Q3

1.  Receptor Preparation

    - Obtain 3D structure from database like PDB
    - Clean and preprocess: remove water molecules and other unwanted ions, co-factors and existing ligands
    - Add hydrogens
    - Assign charges andd atom types: determine protonation states of ionizable residues and atom partial charges at a given pH

2.  Ligand Preparation

    - Obtain or build structure: It can be sourced from databases like PubChem or drawn by chemical sketchers.
    - Generate 3D conformation: If it starts from 2D structures, 3D conformations need to be generated.
    - Optimize Geometry: using quantum mechanics or molecular mechanics methods
    - Assign protonation states, tautomers and partial charges: Similar to the receptor, appropriate protonation states and tautomeric forms for the ligand at the relevant pH are determined.

3.  Grid Generation

    - Define the docking region: This box specifies the volume within which the docking algorithm will search for ligand binding poses.
    - Grid calculation: decide which one to use, targeted docking or blind docking

4.  Running the Docking

    - Searching: Various search algorithms are used, such as genetic algorithms, Monte Carlo methods, or systematic searches.
    - Scoring: A scoring function estimates the binding affinity or the "goodness of fit" between the ligand and the receptor.
    - Ranking Poses: The docking program ranks the generated poses based on their scores, with the top-ranked poses representing the most probable binding modes.
    - Visual Inspection: The top-ranked poses are visually inspected to assess their plausibility.

## Q4

Classic docking tools:

- AutoDock
- DOCK
- Genetic Optimisation for Ligand Docking (GOLD)
- Glide

Co-folding tools:

- AlphaFold-Multimer
- RoseTTAFoldAll
- HADDOCK
- ColabFold

Best tool: AlphaFold3

## Q5

Typical patterns:

- Parallel-displaced: $3.4 \AA$ to $3.8 \AA$ vertically, $1.5 \AA$ to $1.8 \AA$ horizontally
- T-shaped
  - The distance from the center of the top ring (face) to the plane of the stem ring (edge) is $4.5 \AA$ to $5.0 \AA$, for benzene dimer is around $4.96 \AA$.
  - The distance from the interacting hydrogen atoms of the stem ring to the center of the top ring's face is $2.0 \AA$ to $2.5 \AA$
- - Sandwich configuration: The interplanar distance is

## Q6

Providing $\pi$ system: Trp, Tyr, Phe  
Providing cation: Lys, Arg, protonated His

## Q7

Halogen bond is a highly directional, noncovalent interaction where a halogen atom in a molecule R-X acts as an electrophilic species and interacts with a Lewis base.

From the weakest to the strongest: R-Cl < R-Br < R-I

## Q8

Amide-pi interaction is a type of noncovalent attractive force between an amide group and an aromatic $\pi$ system. It is stacked or T-shaped.

This interaction are generally considered to be relatively week, involving electrostatic interactions, dispersion forces and polarization effects.
