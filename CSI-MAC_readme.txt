Updated May 3rd, 2016
Zi-Ming Zhao
ziming.gt@gmail.com


Compile [ OS X 10.9 system]
make
or 
g++ -O3 -o CSI-MAC cMACprf.cpp cPRFCluster.cpp base.cpp -w


Compile [ before OS X 10.9 system]
c++ -O3 -o CSI-MAC cMACprf.cpp cPRFCluster.cpp base.cpp -w


For help:
./CIS-MAC -h



Example runs: 


./CSI-MAC -DC TP53_mut_scaled3_brca_pc_1182_NM_000546_59.txt -DN 771 -S 1 -SC 3 -SR 4.561607e-06 -RF TP53_nonsyn_scaled3_brca_pc_1182_NM_000546_recurrent_50.txt >TP53.brca_pc_louise.CSIMAC.output.txt &



./CIS-MAC -DC NRAS_mut_scaled1_Gilead_Melanoma_570_NM_002524_5.txt -DN 188 -S 1 -SC 1 -SR 8.868393e-06 -RF NRAS_nonsyn_scaled1_Gilead_Melanoma_570_NM_002524_recurrent_4.txt >NRAS.Gilead_Melanoma_CISMAC.output.txt &


./CSI-MAC -DC NF2_mut_scaled1_Meningiomas_1788_NM_000268_1.txt -DN 58 -S 1 -SC 1 -SR 1.249919e-006 -RF NF2_nonsyn_scaled1_Meningiomas_1788_NM_000268_recurrent_0.txt >NF2.Meningiomas.cMACPRF.output.txt &
