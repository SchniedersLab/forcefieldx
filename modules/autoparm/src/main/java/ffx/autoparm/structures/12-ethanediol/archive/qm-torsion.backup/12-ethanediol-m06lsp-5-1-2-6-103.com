%RWF=/scratch/Gau-12-ethanediol/,64GB
%Chk=12-ethanediol-m06lsp-5-1-2-6-103.chk
%Mem=778484KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=64GB

12-ethanediol Rotatable Bond SP Calculation on node7.bme.utexas.edu

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

