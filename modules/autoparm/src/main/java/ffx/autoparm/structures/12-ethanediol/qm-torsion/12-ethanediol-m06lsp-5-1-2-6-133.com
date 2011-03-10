%RWF=/scratch/Gau-12-ethanediol/,32GB
%Nosave
%Chk=12-ethanediol-m06lsp-5-1-2-6-133.chk
%Mem=389242KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=32GB

12-ethanediol Rotatable Bond SP Calculation on node9.bme.utexas.edu

0 1
 C    0.403141    0.647792    0.260116
 C   -0.403141   -0.647792    0.260116
 H    0.072160    2.484222   -0.264829
 H   -0.072160   -2.484222   -0.264829
 O   -0.403141    1.665751   -0.275974
 O    0.403141   -1.665751   -0.275974
 H    0.730658    0.904314    1.264420
 H    1.284821    0.489869   -0.352497
 H   -0.730658   -0.904314    1.264420
 H   -1.284821   -0.489869   -0.352497

