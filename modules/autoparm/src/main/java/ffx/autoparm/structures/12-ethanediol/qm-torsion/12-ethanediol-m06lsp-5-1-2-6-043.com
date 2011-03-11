%RWF=/scratch/Gau-12-ethanediol/,32GB
%Nosave
%Chk=12-ethanediol-m06lsp-5-1-2-6-043.chk
%Mem=389242KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=32GB

12-ethanediol Rotatable Bond SP Calculation on node9.bme.utexas.edu

0 1
 C    0.177033    0.738911    0.589283
 C   -0.177033   -0.738911    0.589283
 H    0.080904    2.221395   -0.638077
 H   -0.080904   -2.221395   -0.638077
 O   -0.177033    1.310956   -0.634430
 O    0.177033   -1.310956   -0.634430
 H   -0.333797    1.214723    1.425139
 H    1.248080    0.845857    0.752682
 H    0.333797   -1.214723    1.425139
 H   -1.248080   -0.845857    0.752682

