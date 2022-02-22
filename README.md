# Protein-Hydrogen-Bond-Analysis

A python based routine for the analysis of the hydrogen bond network in a protein structure.

The intended use is for analysis of Molecular Dynamics trajectories, though in principle, a single protein structure file will work as well (the only caveat being that the lifetime value would be meaningless). 

The routine focuses on protein backbone and side-chain hydrogen bonds (presumably it could be easily amended to include protein bound water, should one have the interest).

As input, the routine takes protein structure coordinate file(s) in the PDB.

The results include the Hydrogen Bond energy and Bond lifetime.
