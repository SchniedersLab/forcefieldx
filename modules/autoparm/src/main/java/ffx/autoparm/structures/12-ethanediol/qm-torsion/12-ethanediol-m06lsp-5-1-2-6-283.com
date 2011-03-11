%RWF=/scratch/Gau-12-ethanediol/,32GB
%Nosave
%Chk=12-ethanediol-m06lsp-5-1-2-6-283.chk
%Mem=389242KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=32GB

12-ethanediol Rotatable Bond SP Calculation on node9.bme.utexas.edu

0 1
 C   -0.294045    0.696082    0.508405
 C    0.294045   -0.696082    0.508405
 H   -0.207627    2.221625   -0.684579
 H    0.207627   -2.221625   -0.684579
 O    0.294045    1.434936   -0.527415
 O   -0.294045   -1.434936   -0.527415
 H   -1.370308    0.619516    0.382471
 H   -0.096175    1.165306    1.471004
 H    1.370308   -0.619516    0.382471
 H    0.096175   -1.165306    1.471004

