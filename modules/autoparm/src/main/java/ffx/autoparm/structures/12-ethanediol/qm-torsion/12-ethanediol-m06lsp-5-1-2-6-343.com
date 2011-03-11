%RWF=/scratch/Gau-12-ethanediol/,32GB
%Nosave
%Chk=12-ethanediol-m06lsp-5-1-2-6-343.chk
%Mem=389242KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=32GB

12-ethanediol Rotatable Bond SP Calculation on node9.bme.utexas.edu

0 1
 C   -0.074320    0.765840    0.618886
 C    0.074320   -0.765840    0.618886
 H   -0.010458    2.217283   -0.638982
 H    0.010458   -2.217283   -0.638982
 O    0.074320    1.275260   -0.672496
 O   -0.074320   -1.275260   -0.672496
 H   -1.057653    1.028348    1.004854
 H    0.660847    1.188323    1.300778
 H    1.057653   -1.028348    1.004854
 H   -0.660847   -1.188323    1.300778

