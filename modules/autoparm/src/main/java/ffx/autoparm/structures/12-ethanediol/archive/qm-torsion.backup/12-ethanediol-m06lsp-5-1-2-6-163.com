%RWF=/scratch/Gau-12-ethanediol/,64GB
%Chk=12-ethanediol-m06lsp-5-1-2-6-163.chk
%Mem=778484KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=64GB

12-ethanediol Rotatable Bond SP Calculation on node7.bme.utexas.edu

0 1
 C    0.427811    0.625562    0.095959
 C   -0.427811   -0.625562    0.095959
 H    0.065507    2.524830   -0.051274
 H   -0.065507   -2.524830   -0.051274
 O   -0.427811    1.719234   -0.106113
 O    0.427811   -1.719234   -0.106113
 H    0.972238    0.710963    1.032739
 H    1.153262    0.542060   -0.708313
 H   -0.972238   -0.710963    1.032739
 H   -1.153262   -0.542060   -0.708313

