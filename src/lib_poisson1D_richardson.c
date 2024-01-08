/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>

static int LAB_init = 3; //lab de AB
static int KU_init = 1; //ku de AB
static int KL_init = 1; //kl de AB

void eig_poisson1D(double* eigval, int *la){
    for(int j = 1; j<=*la ; ++j){
        eigval[j-1] = 2*(1-cos((M_PI*j)/ (*la+1)));
    }
}

double eigmax_poisson1D(int *la){
    return 2*(1-cos((M_PI*(*la))/ (*la+1)));
}

double eigmin_poisson1D(int *la){
    return 2*(1-cos((M_PI* 1)/ (*la+1)));
}

double richardson_alpha_opt(int *la){
    return 2/(eigmin_poisson1D(la)+eigmax_poisson1D(la));
}


void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    
    double *b = (double *)malloc((*la) * sizeof(double));
    memset(b,0, sizeof(double) * (*la));
    
    //calcul b_exp initial b0
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1, AB, *lab, X, 1, 0, b, 1);
    
    
    //calcul residu relatif initial
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
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1, AB, *lab, X, 1, 0, b, 1);
        
        
        *nbite = *nbite + 1;
    }
    
    free(b);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    int ii, jj, kk;
    for (jj=0;jj<(*la);jj++){
        kk = jj*(*lab);
        MB[kk+ *kv]=1/AB[kk+ *kv];
    }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    
    /*calcul (D-E)*/
    
    int ii, jj, kk;
    for (jj=0;jj<(*la);jj++){
        kk = jj*(*lab);
        if (*kv>=0){
            for (ii=0;ii< *kv;ii++){
                MB[kk+ii]=0.0;
            }
        }
        //diag 
        MB[kk+ *kv]=AB[kk+ *kv];
        for (int i=1;i<= *kl;i++){
            //sub diagonales
            MB[kk+ *kv +i]= AB[kk+ *kv +i]; 
        }
    }
    
    /*calcul (D-E)^-1*/
    
    
    /*Initialization*/
    
    //facteurs de permutation pour la décomposition LU
    int *ipiv = (int *)malloc((*la)* sizeof(int));
    memset(ipiv,0, sizeof(int) * (*la));
    
    //matrice idendite & resultat
    double *I_inv = (double *)malloc((*la+1) * (*la) * sizeof(double));
    memset(I_inv,0, sizeof(double) * (*la+1)*(*la));
    for (int jj=0;jj<(*la);jj++){
        I_inv[jj*(*la+1)]=1;
    }
    
    
    //autres
    int info = 0;
    int ku_new = 0;
    
    //Partant de la DEF de Matrice Inverse
    //Résoudre le système d'équations AX = I_inv, ici A = MB
    //la solution X = inverse de A, stocke dans I_inv
    
    dgbtrf_(la, la, kl, &ku_new, MB, lab, ipiv, &info);
    
    if (info != 0) {
        perror("MB non LU ou error !\n");
        printf("info = %d \n", info);
        return;
    }
    
    dgbtrs_("N", la, kl, &ku_new, la, MB, lab, ipiv, I_inv, la, &info, 1);
    if (info != 0) {
        perror("MB non inversible ou error !\n");
        printf("info = %d \n", info);
        return;
    }  
    
    
    /*DONE I_inv = (D-E)^-1*/
    
    //update MB
    for (ii=0;ii<(*la);ii++){
        for (jj=0;jj<(*la);jj++){
            MB[ii*(*la)+jj] = I_inv[ii*(*la+1)+jj];
        }
    }
    
    free(ipiv);
    free(I_inv);  
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    
    double *b = (double *)malloc((*la) * sizeof(double));
    memset(b,0, sizeof(double) * (*la));
    
    //calcul b_exp initial b0
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, KL_init, KU_init, 1, AB, LAB_init, X, 1, 0, b, 1);
    
    
    //calcul residu relatif initial
    double r_relatif = relative_forward_error(b,RHS,la);
    
    //Richardson
    double r[*la];
    
    while (r_relatif > *tol && *nbite < *maxit) {
        
        //calcul b - AB*X = RHS - b
        for (int i = 0; i < *la; ++i) {
            r[i] = RHS[i]-b[i];
        }
        
        //calcul sol X
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1, MB, *lab, r, 1, 1, X, 1);
        
        
        //residu update
        r_relatif = relative_forward_error(b,RHS,la);
        resvec[*nbite] = r_relatif;
        
        //b_exp update
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, KL_init, KU_init, 1, AB, LAB_init, X, 1, 0, b, 1);
        
        
        *nbite = *nbite + 1;
    }
    
    
    free(b);
}

