%RWF=/scratch/Gau-12-ethanediol/,64GB
%Chk=12-ethanediol-m06lsp-5-1-2-6-223.chk
%Mem=778484KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=64GB

12-ethanediol Rotatable Bond SP Calculation on node7.bme.utexas.edu

0 1
 C   -0.408555    0.643721    0.234947
 C    0.408555   -0.643721    0.234947
 H   -0.066036    2.496548   -0.221966
 H    0.066036   -2.496548   -0.221966
 O    0.408555    1.678113   -0.250635
 O   -0.408555   -1.678113   -0.250635
 H   -1.265099    0.492147   -0.414080
 H   -0.776244    0.874531    1.231446
 H    1.265099   -0.492147   -0.414080
 H    0.776244   -0.874531    1.231446

