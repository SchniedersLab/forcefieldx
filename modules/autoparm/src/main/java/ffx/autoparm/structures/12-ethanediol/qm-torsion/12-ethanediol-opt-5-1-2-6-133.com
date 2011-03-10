%RWF=/scratch/Gau-12-ethanediol/,32GB
%Nosave
%Chk=12-ethanediol-opt-5-1-2-6-133.chk
%Mem=389242KB
%Nproc=1
#HF/6-31G* Opt=ModRed MaxDisk=32GB

12-ethanediol Rotatable Bond Optimization on node8.bme.utexas.edu

0 1
 C    0.354013    0.673872    0.407494
 C   -0.354013   -0.673872    0.407494
 H    0.143069    2.355805   -0.536582
 H   -0.143069   -2.355805   -0.536582
 O   -0.354013    1.558163   -0.423555
 O    0.354013   -1.558163   -0.423555
 H    0.393872    1.071071    1.419614
 H    1.371302    0.526550    0.060448
 H   -0.393872   -1.071071    1.419614
 H   -1.371302   -0.526550    0.060448

5 1 2 6     132.64 F

