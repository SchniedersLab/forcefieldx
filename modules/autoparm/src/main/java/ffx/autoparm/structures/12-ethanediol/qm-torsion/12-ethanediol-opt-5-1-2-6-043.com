%RWF=/scratch/Gau-12-ethanediol/,32GB
%Nosave
%Chk=12-ethanediol-opt-5-1-2-6-043.chk
%Mem=389242KB
%Nproc=1
#HF/6-31G* Opt=ModRed MaxDisk=32GB

12-ethanediol Rotatable Bond Optimization on node9.bme.utexas.edu

0 1
 C    0.280458    0.701252    0.523262
 C   -0.280458   -0.701252    0.523262
 H    0.209636    2.209163   -0.690751
 H   -0.209636   -2.209163   -0.690751
 O   -0.280458    1.413052   -0.544781
 O    0.280458   -1.413052   -0.544781
 H    0.039906    1.177670    1.472709
 H    1.362249    0.644112    0.436718
 H   -0.039906   -1.177670    1.472709
 H   -1.362249   -0.644112    0.436718

5 1 2 6      42.64 F

