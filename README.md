# MCPBy-param-for-multiple-chains
Setting Up MCPBy Parameters for Pyruvate Kinase with Multiple Ligands and Metals
<br> This tutorial guides you through setting up MCPBy parameters for a Pyruvate Kinase system containing multiple ligands (PEP, FBP) and metal ions (Mg²⁺, K⁺).
<br>
Important Note: While some steps might seem unnecessary, they are crucial for successful parameter generation.

# Software Used:
Maestro (or similar protein preparation tool)
pdb4amber
Gaussian
MCPBy (MCPB.py)
Steps:
# Protein and Ligand Preparation:
Clean and minimize: Use Maestro or a similar tool to prepare the protein structure (FBP_PEP_8IAX.pdb). Assign protonation states at pH 7. Save the final complex without hydrogens (FBP_PEP_8IAX_prep_noh.pdb) using pdb4amber.

# Extract Ligands and Protein:
Extract FBP, PEP, and protein from the complex PDB file, resulting in four separate files: <br>
FBP_all.pdb
PEP_all.pdb
MG_all.pdb # (Magnesium ion file)
protein_noH.pdb # (Protein without hydrogens)

# Add Hydrogens:
Independently add hydrogens to each extracted file (except MG_all.pdb as metals don't have hydrogens) using pdb4amber: <br>
FBP_all_H.pdb
PEP_all_H.pdb
protein_H.pdb

# Open Ligands in Gaussian (Optional but highly recommended):
Open each ligand file (FBP_all.pdb, PEP_all.pdb) in Gaussian to check and potentially modify protonation states if needed.


# Merge and Renumber:

Combine all prepared files into a single file named complex.pdb:
<br>
bash : cat protein_H.pdb FBP_all_H.pdb PEP_all_H.pdb MG_all.pdb > complex.pdb
Use code with caution.
content_copy
Renumber the atoms in complex.pdb using pdb4amber:

Bash
pdb4amber -i complex.pdb -o complex_prep.pdb
Use code with caution.
content_copy
<br>
# Identify Metal Ion IDs and generate fake parameters:

Use grep MG complex_prep.pdb to find the residue IDs for the Magnesium ions.
Generate Initial Input File:
original_pdb complex_prep.pdb

<br>
Create an input file (input.in) for MCPBy specifying the system details and metal ion IDs.
( group_name PYK84
cut_off 2.8
ion_ids 31079 31080 31081 31082
software_version g09
ion_mol2files MG.mol2
naa_mol2files PEP.mol2 HOH.mol2 FBP.mol2
frcmod_files PEP.frcmod )
<br>
Run MCPBy (MCPB.py -i input.in -s 1) to generate a sample input file for further editing.
Energy Minimization (simple SCF calculation) and Parameter Generation: these are not real parameters we are going to use.
Stage 1: Run MCPBy with the sample input file (MCPB.py -i input.in -s 2) to perform energy minimization and generate initial parameter files.

# Stage 2 (Repeat for MG1-MG4): these will be real parameters. 
Create four directories named MG1, MG2, MG3, and MG4.
<br> Copy all files generated in Stage 1 into each directory (MG1-MG4).
Edit the input.in file in each directory to specify a different metal ion (MG1, MG2, MG3, MG4) based on their IDs.
Run MCPBy in each directory for parameter generation specific to each metal ion: MCPB.py -i input.in -s 1
This will generate separate parameter files for each metal ion (PYK84_mcpby.pdb) in their respective folders.
**Geometry Optimization:**

Follow a stepwise geometry optimization process for each metal ion directory (MG1-MG4) using Gaussian:
Step 1: PM7 with harmonic restraint (1000 kcal/mol*Å²)
Step 2: PM7 with frequency calculation (freq) and a lower harmonic restraint (100 kcal/mol*Å²)
Analyze frequencies in Gaussian to check for clashes.
If necessary, modify coordinates slightly to avoid inappropriate hydrogen bonding.
Step 3: Re-optimize at the same level (PM7)
Step 4: DFT optimisation with tight convergence criteria optimization with frequency calculation (opt freq)


# MCPBY file generatio for the individual metal unit:
After successfully able to complete optimisation and frequency calculation without any imaginary frequency, you are ready to go to next step.
MCPBY MCPB.py -i input.in -s 2
MCPBY MCPB.py -i input.in -s 3
MCPBY MCPB.py -i input.in -s 4

# Merge files: Copy all optimized molecule files (*.mol2) from MG1-MG4 folders into a single "merge" folder.
Copy one of the PYK84_mcpby.frcmod files ( force constant and geometrical parameters).

# Final parameterisation of all metal units in all chains1: now we have everything in merge folder, we will work on all the files in this folder.
In PYK84_mcpby.pdb (copied and generated in step 7), we have ligands as Ap1,AP2,AP3,AP4 and similarly for MG1, PP and GU their correspodning atomtypes are availabel in mol2 files. check the correspondence of *.mol2 and .frcmod file and the coordinate.

MCPB.py -i input.in -s 2
MCPB.py -i input.in -s 3
MCPB.py -i input.in -s 4

# load the correct forcefield in tleap and set the box, insert parameters for FBP and ready to generate the final PYK84_solv.inpcrd and PYK84_solv.prmtop files
calculate the buffer ions according to the electrolyte concenetration of the simulation box. for me, system has -80 charge and i added these lines in tleap script.
solvatebox protein OPCBOX 14
addions protein Na+ 227
addions protein Cl- 147
 
