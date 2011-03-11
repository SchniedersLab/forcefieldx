%RWF=/scratch/Gau-12-ethanediol/,64GB
%Chk=12-ethanediol-m06lsp-5-1-2-6-193.chk
%Mem=778484KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=64GB

12-ethanediol Rotatable Bond SP Calculation on node7.bme.utexas.edu

0 1
 C   -0.429680    0.623579    0.069804
 C    0.429680   -0.623579    0.069804
 H   -0.068690    2.525587   -0.033356
 H    0.068690   -2.525587   -0.033356
 O    0.429680    1.722518   -0.077625
 O   -0.429680   -1.722518   -0.077625
 H   -1.131666    0.558761   -0.756915
 H   -1.001184    0.682361    0.992446
 H    1.131666   -0.558761   -0.756915
 H    1.001184   -0.682361    0.992446

