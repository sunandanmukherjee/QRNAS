QRNA 0.3 - Quick Refinement of Nucleic Acids 0.3
Julek Stasiewicz 2012

The command to build multiprocessor version of QRNA is:
  $ ./qrnamake parallel
It should work only on POSIX-compatible systems.
Single-core machines and non-POSIX users should build single-threaded version by typing:
  $ ./qrnamake sequential
It should go quietly without warnings.

Example use cases:

  $ ./QRNA -i pdbfile.pdb
  It minimizes pdbfile.pdb and writes pdbfile_out.pdb every 100 steps. 
  All default parameters are used, which should be fine in most cases.

  $ ./QRNA -i pdbfile.pdb -o outfile.pdb
  It minimizes pdbfile.pdb and writes outfile.pdb every 100 steps. 
  All default parameters are used, which should be fine in most cases.

  $ ./QRNA -m restrfile.txt -i pdbfile.pdb
  It minimizes pdbfile.pdb and writes pdbfile_out.pdb every 100 steps.
  Restraints on secondary structure can be specified in restrfile.txt.
  Pairwise restraints on atoms' distances can be also specified there.
  All default parameters are used.
  
  $ ./QRNA -c configfile.txt
  In this case `INPUTPDB' and `OUTPUTPDB' must be specified in configfile.txt 
  (perhaps among other parameters). `RESTRFILE' can be also specified there.

  $ ./QRNA -i pdbfile.pdb -o outfile.pdb -c configfile.txt
  It minimizes pdbfile.pdb and writes pdbfile_out.pdb.
  It uses parameters from configfile.txt, but overrides `INPUTPDB' and `OUTPUTPDB' 
  (if they collide with pdbfile.pdb or outfile.pdb respectively). 
  Same thing applies to commandline option `-m restrfile.txt' vs. `RESTRFILE otherrestrfile.txt' in configfile.txt.
  
  $ ./QRNA -P -i pdbfile.pdb
  It minimizes pdbfile.pdb and writes pdbfile_out.pdb every 100 steps. 
  All default parameters are used.
  Occupancy and beta-factors in pdbfile.pdb are treated as restraint information exactly as by SimRNA:
  Occupancy = 1.00        => no restraint
  0.00 < Occupancy < 1.00 => restraint enabled; then restrain_spring_constant = B-factor ^-1 [kcal/mol/A^2]
  Occupancy = 0.00        => atom completely fixed.

For details on config file format please see configfile.txt
For details on restraints file format please see restraints.txt

Please note that QRNA executable cannot be moved around without the accompanying ./forcefield directory!
