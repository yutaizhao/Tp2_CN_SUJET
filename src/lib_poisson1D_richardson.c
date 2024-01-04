/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>

void eig_poisson1D(double* eigval, int *la){
    for(int j = 1; j<=*la ; ++j){
        eigval[j-1] = 2*(1-cos((M_PI*j)/ *la));
    }
}

double eigmax_poisson1D(int *la){
    return 2*(1-cos((M_PI*(*la))/ *la));
}

double eigmin_poisson1D(int *la){
    return 2*(1-cos((M_PI* 1)/ *la));
}

double richardson_alpha_opt(int *la){
  return 2/(eigmin_poisson1D(la)+eigmax_poisson1D(la));
}


void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {


    
    double *b = (double *)malloc((*la) * sizeof(double));
    memset(b,0, sizeof(double) * (*la));
    
    //calcul residu initial r0
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, *alpha_rich, AB, *lab, X, 1, 0, b, 1);
    

    //calcul residu relatif
    double r_relatif = relative_forward_error(b,RHS,la);
    
    

    //Richardson
    while (r_relatif > *tol && *nbite < *maxit) {
        
        //calcul sol
        for (int i = 0; i < *la; ++i) {
            X[i] = X[i] + (*alpha_rich) * (RHS[i]-b[i]);
        }

        //residu update
        r_relatif = relative_forward_error(b,RHS,la);
        resvec[*nbite] = r_relatif;
        
        //b_exp update
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, *alpha_rich, AB, *lab, X, 1, 0, b, 1);
        
        
        *nbite = *nbite + 1;
    }
    
    free(b);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

