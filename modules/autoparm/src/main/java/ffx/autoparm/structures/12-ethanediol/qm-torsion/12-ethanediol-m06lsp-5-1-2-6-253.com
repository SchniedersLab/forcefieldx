%RWF=/scratch/Gau-12-ethanediol/,32GB
%Nosave
%Chk=12-ethanediol-m06lsp-5-1-2-6-253.chk
%Mem=389242KB
%Nproc=1
#M06L/6-31G** SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=32GB

12-ethanediol Rotatable Bond SP Calculation on node9.bme.utexas.edu

0 1
 C    0.645614    0.385937   -0.405105
 C   -0.645614    0.385937    0.405105
 H    2.368627   -0.496813   -0.280285
 H   -2.368627   -0.496813    0.280285
 O    1.598682   -0.401999    0.262258
 O   -1.598682   -0.401999   -0.262258
 H    0.427689   -0.004788   -1.393415
 H    1.018020    1.401970   -0.515755
 H   -0.427689   -0.004788    1.393415
 H   -1.018020    1.401970    0.515755

