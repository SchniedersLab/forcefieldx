%RWF=/scratch/Gau-12-ethanediol/,32GB
%Nosave
%Chk=12-ethanediol-opt-5-1-2-6-193.chk
%Mem=389242KB
%Nproc=1
#HF/6-31G* Opt=ModRed MaxDisk=32GB

12-ethanediol Rotatable Bond Optimization on node9.bme.utexas.edu

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

5 1 2 6     192.64 F

