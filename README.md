# omp_for
there are two files in this branch


"RK45_w_tur_Cfin.py" is the one to execute and call the function in "cf_genmat_para_fin.so"


"cf_genmat_para_fin.c" should be compiled by the following command 
"icc -fPIC -shared -o cf_genmat_para_fin.so cf_genmat_para_fin.c -qopenmp" 
thus generate "cf_genmat_para_fin.so"


detailed about "RK45_w_tur_Cfin.py"
