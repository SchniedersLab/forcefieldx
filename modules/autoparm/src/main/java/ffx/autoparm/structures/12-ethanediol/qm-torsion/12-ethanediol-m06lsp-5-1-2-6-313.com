%RWF=/scratch/Gau-12-ethanediol/,32GB
%Nosave
%Chk=12-ethanediol-m06lsp-5-1-2-6-313.chk
%Mem=389242KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=32GB

12-ethanediol Rotatable Bond SP Calculation on node9.bme.utexas.edu

0 1
 C   -0.194940    0.732875    0.581177
 C    0.194940   -0.732875    0.581177
 H   -0.103191    2.220285   -0.642746
 H    0.103191   -2.220285   -0.642746
 O    0.194940    1.322443   -0.623533
 O   -0.194940   -1.322443   -0.623533
 H   -1.273624    0.810890    0.704748
 H    0.273566    1.213155    1.439200
 H    1.273624   -0.810890    0.704748
 H   -0.273566   -1.213155    1.439200

