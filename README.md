The package includes the R source code of drug-target interaction prediction algorithm (RLS-KF) [1], which combines the NII [2] 
and Kernel Fusion [3] technique.

How To Run:

(1) Download benchmark data sets from http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/ [4].
  
  Total 12 files:
 
 [e_admat_dgc.txt, e_simmat_dc.txt, e_simmat_dg.txt, gpcr_admat_dgc.txt, gpcr_simmat_dc.txt,
  gpcr_simmat_dg.txt, ic_admat_dgc.txt, ic_simmat_dc.txt, ic_simmat_dg.txt, nr_admat_dgc.txt,
  nr_simmat_dc.txt, nr_simmat_dg.txt]
 
(2) Compile the *.cpp files. For example, on Windows 7, one can use 'Rtools: https://cran.r-project.org/bin/windows/Rtools/'

(3) Install required packages

(4) setwd("dir including source codes")

(5) source("demo-RLS_KF.R") 
 
 One can modify the variable "partfn" value in the 'demo-RLS_KF.R' file.
 
 Default: partfn = "nr", which gives the results for NR data set.

 These codes were tested on Windows 7, but should also work on Linux.

Attention: 
- This package is free for academic usage ONLY. 

- This package was developed by Dr. Ming Hao (kevin.m.hao@gmail.com). For any problem about these codes, 
  please feel free to contact Dr. Hao.

References:

[1] M Hao, et al., Improved prediction of drug-target interactions using regularized least squares integrating with kernel fusion technique. Analytica Chimica Acta 909 (2016) 41-50.

[2] JP Mei, et al., Bioinformatics 29 (2013) 238-245.

[3] B Wang, et al., Nat. Methods 11 (2014) 333-337.

[4] Y Yamanishi, et al., Bioinformatcs 24 (2008) i232-i240.
