#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
//#include <complex.h>
/*icc -fPIC -shared -o cf_genmat_para_fin.so cf_genmat_para_fin.c -qopenmp*/
// icc cf_genmat.c -o exe_cf_genmat
//#define ms 3 //this will make all the ms to 3 without the need of declaration



///////////////////////////////////////////////////////////
#define dk   1000  //
#define hk 600000 ///6*10**5
#define lk -600000
#define nk (2*hk/dk+1)
#define ms (nk * 4)
//#define ConsMu 0 // 0 false 1 true       0 means false for simplicity
//#define ConsMu
//--------------------------------------
#define CC 0.001
#define kcut 0.01
#define alpha 5.0/3.0
///////////////////////////////////////////////////////////
int rms(int i){
    return ms;
}
int rdk(int i){
    return dk;
}
int rhk(int i){
    return hk;
}
int rnk(int i){
    return nk;
}
int rConsMu(int i){
    int rcs;
    #ifdef ConsMu
        rcs = 1;
    #else
        rcs = 0;
    #endif
    return rcs;
}
/*double rmuz(double z){
    double mu = pow(10.0,4.0) * exp(-0.3 * z);
    return mu;
}*/
/////////////////////////////////////////////////////////////////
double fmuz(double z){
    double out;
    /*if (ConsMu == 0){
        out = pow(10.0,4.0) * exp(-0.3 * z);
    }
    else{
        out = pow(10.0,3.0);
    }*/
    #ifdef ConsMu
        out = pow(10.0,3.0);
    #else
        out = pow(10.0,4.0) * exp(-0.3 * z);
    #endif
    return out;
}

double lda0(double z){
    return 300 * fmuz(z);
}

double ldkx(double kx, double z){
    //double muz_ldkx = fmuz(z);
    double out;
    /*if(kx == 0.0){
        out = lda0(z) / M_PI / 2.0;
    }
    else{
        //out = CC * lda0(z) * (pow((abs(kx) / kcut) , -alpha/2.0)) /2.0 /M_PI;

    }*/
    //out = CC * 300 * (pow((abs(kx) / kcut) , -alpha/2.0)) /2.0 /M_PI;
    return out;
}

/////////////////////////////////////
typedef struct pointda {
    double A1d[ms];
    //double y;
} POINTDA;
POINTDA *get_pointda(
    double *y, double z, double ar, double omg, double omga, double Omg, double vxp, double vxn, double vz
    )  // pointer to double array const false real
{   //clock_t st = clock();
    //time_t st = time(NULL);
    //printf("para_fin\n");
    double muz = fmuz(z);
    double hvv  = muz * (1-(vz * vz + vxp * vxn));
    double hvva = -hvv * ar;
    double Lda = 2.0 * (1-ar) * hvv;
    const double one_2pi = 1.0 /2.0 /M_PI;
    const double one_kcut = 1.0/ kcut;
    const double two_alpha = -alpha/2.0;
    double factors = CC * 300 * muz * one_2pi;
    //------------------------------------------------
    double* pLHSM = calloc(ms * ms, sizeof(double));   // init to 0 automatically
    int i,j,k;
    POINTDA *dqdt = malloc(sizeof(POINTDA));

    //printf("pre para\n");
    #pragma omp parallel private( i, j, k)
    {
        //printf("%i\n",k);

        //////////////////////////////////////////
        #pragma omp for
        for(i=0;i<ms;i+=4){

            j = i;
            int modk = lk + (i/4) * dk;
            //back up 1/////////
            pLHSM[i * ms + j] =  (-Omg + omg + vxp * modk) + Lda;
            pLHSM[i * ms + j + 1] = -hvv;
            pLHSM[i * ms + j + 3] = -hvva;

            //////////////
            pLHSM[(i+1) * ms + j +1] = (-Omg + omg + vxn * modk) + Lda;
            pLHSM[(i+1) * ms + j+1 - 1] = -hvv;
            pLHSM[(i+1) * ms + j+1 +1] = -hvva;

            ////////////
            pLHSM[(i+2) * ms + j+2] = (-Omg + omga + vxp * modk) + Lda;
            pLHSM[(i+2) * ms + j+2 - 1] = -hvv;
            pLHSM[(i+2) * ms + j+2 + 1] = -hvva;

            //////////
            pLHSM[(i+3) * ms + j+3] =  (-Omg + omga + vxn * modk) + Lda;
            pLHSM[(i+3) * ms + j+3 - 1] = -hvva;
            pLHSM[(i+3) * ms + j+3 - 3] = -hvv;

            for (j=0;j<ms;j+=4){

                //int modki = lk + (i/4) * dk;
                //int modkj = lk + (j/4) * dk;
                double MID = abs( (i-j)/4*dk );
                //printf("%.17g abs MID\n ",MID);
                for(k=0;k<4;k++){
                    #ifdef ConsMu
                        pLHSM[(i+k) * ms + j+k] += factors * (pow(((MID+0.00000251188643150981) * one_kcut) , two_alpha));
                    #else                                                 //        lda0  fmuz             +0.00000251189
                                                                    //the ugly number is need for the modification of the shift of the x
                        pLHSM[(i+k) * ms + j+k] += factors * (pow(((MID+0.00000251188643150981) * one_kcut) , two_alpha));
                    #endif
                }

            }
        }
        //printf("1st for\n");
        // LHSM is built
        //////////////////////////////////////////////////////////////////////////////
        // start to mutilplied 2 matrix
        #pragma omp for
        for(i=0;i<ms;i++){
            dqdt->A1d[i] = 0.0;
        }
        //printf("2nd for\n");

        #pragma omp for
        for(i=0;i<ms;i++){
            double A1d_i = 0;
            //#pragma omp for reduction(+:A1d_i)
            for(j = 0;j<ms;j++){
                //printf("%f, ",pLHSM[i * ms + j]);
                A1d_i += pLHSM[i * ms + j] * y[j] / vz;
                //printf("%.17g,  ",pLHSM[i * ms + j]);
            }
            //printf("\n");
            //printf("\n");
            dqdt->A1d[i] = A1d_i;
        }
    }
    free(pLHSM);
    /*for(i=0;i<ms;i++){
        printf("%f\n",dqdt.A1d[i]);
    }*/

    /*clock_t et = clock();
    double tuse = (double)(et - st) / CLOCKS_PER_SEC;
    printf("%f time used\n",tuse);*/
    //printf("here hrer\n");
    return dqdt;
}

void free_pointda(POINTDA *pda)
{
    //printf("%s\n","deleteed");
    free(pda);
}

////////////////////////////////////////////////////////////////////////////////////////
int main(){
    int i,j;
    /*for(i=0;i<ms;i++){
        for(j=0;j<ms;j++){
            printf("%f\n",LHSM[i][j]);
        }
    }*/

    printf("%s\n","started");
    printf("%s\n","started");
    return 0;
}
