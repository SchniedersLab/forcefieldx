%RWF=/scratch/Gau-12-ethanediol/,64GB
%Chk=12-ethanediol-m06lsp-5-1-2-6-073.chk
%Mem=778484KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=64GB

12-ethanediol Rotatable Bond SP Calculation on node7.bme.utexas.edu

0 1
 C    0.280458    0.701252    0.523262
 C   -0.280458   -0.701252    0.523262
 H    0.209636    2.209163   -0.690751
 H   -0.209636   -2.209163   -0.690751
 O   -0.280458    1.413052   -0.544781
 O    0.280458   -1.413052   -0.544781
 H    0.039906    1.177670    1.472709
 H    1.362249    0.644112    0.436718
 H   -0.039906   -1.177670    1.472709
 H   -1.362249   -0.644112    0.436718

