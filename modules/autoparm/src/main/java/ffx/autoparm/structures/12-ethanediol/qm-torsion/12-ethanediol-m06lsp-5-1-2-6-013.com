%RWF=/scratch/Gau-12-ethanediol/,32GB
%Nosave
%Chk=12-ethanediol-m06lsp-5-1-2-6-013.chk
%Mem=389242KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=32GB

12-ethanediol Rotatable Bond SP Calculation on node9.bme.utexas.edu

0 1
 C    0.054268    0.769013    0.621655
 C   -0.054268   -0.769013    0.621655
 H    0.006924    2.216437   -0.640567
 H   -0.006924   -2.216437   -0.640567
 O   -0.054268    1.272644   -0.675931
 O    0.054268   -1.272644   -0.675931
 H   -0.720851    1.175581    1.267840
 H    1.011870    1.058592    1.050247
 H    0.720851   -1.175581    1.267840
 H   -1.011870   -1.058592    1.050247

