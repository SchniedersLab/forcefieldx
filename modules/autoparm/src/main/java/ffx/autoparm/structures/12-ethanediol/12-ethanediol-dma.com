%RWF=/scratch/,-1,/home/scratch/,500GB
%Nosave
%Chk=12-ethanediol-dma.chk
%Mem=678484KB
%Nproc=1
#MP2/6-311G** Sp Density=MP2 MaxDisk=500GB

 12-ethanediol Gaussian SP Calculation on node9.bme.utexas.edu

0  1
C          	       0.280458    0.701252    0.523262
C         	      -0.280458   -0.701252    0.523262
H          	       0.209636    2.209163   -0.690751
H         	      -0.209636   -2.209163   -0.690751
O         	      -0.280458    1.413052   -0.544781
O          	       0.280458   -1.413052   -0.544781
H          	       0.039906    1.177670    1.472709
H          	       1.362249    0.644112    0.436718
H         	      -0.039906   -1.177670    1.472709
H         	      -1.362249   -0.644112    0.436718
	
