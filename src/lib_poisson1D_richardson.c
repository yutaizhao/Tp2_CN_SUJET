/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>

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
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	MB[kk+ii]=0.0;
      }
    }
    MB[kk+ *kv]=AB[kk+ *kv];
  }
  MB[0]=0.0;
  MB[(*lab)*(*la)-1]=0.0;
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){

 int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	MB[kk+ii]=0.0;
      }
    }
    MB[kk+ *kv]=AB[kk+ *kv];
    for (int i=1;i<= *kl;i++){
	MB[kk+ *kv -i]= - AB[kk+ *kv -i];
      }
  }
  MB[0]=0.0;
  MB[(*lab)*(*la)-1]=0.0;

}

/*
void richardson_MB_jac(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
 	
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
            printf("%f \n",1/MB[i*(*lab)+(*kl)]);
            X[i] = X[i] + (1/MB[i*(*lab)+(*kl)]) * (RHS[i]-b[i]);
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
*/

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
 	
    

    //calcul matrice inverse



    // Facteurs de permutation pour la décomposition LU

    int *ipiv = (int *)malloc((*la)* sizeof(int));

    memset(ipiv,0, sizeof(int) * (*la));

   

    //matrice idendite & resultat



    double *I_inv = (double *)malloc((*lab) * (*la) * sizeof(double));

    int ii, jj, kk;

  for (jj=0;jj<(*la);jj++){

    kk = jj*(*lab);

    I_inv[kk+ 1]=1;

  }

    // Résoudre le système d'équations AX = I_inv, ici A = MB

    //la solution X = inverse de A, stocke dans I_inv

    int info = 0;
    int ku_new = 0;
    int lab_new = 0;
    dgbsv_(la, kl, &ku_new, la, MB, lab, ipiv, I_inv , la, &info);


    if (info != 0) {

    	perror("MB non inversible ou error !\n");

    	printf("info = %d \n", info);

    	return;

    }    
    
    
  


    double *b = (double *)malloc((*la) * sizeof(double));

    memset(b,0, sizeof(double) * (*la));

        

    //calcul b_exp initial b0

    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1, AB, *lab, X, 1, 0, b, 1);





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

        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1, I_inv, *lab, r, 1, 1, X, 1);
    


        //residu update

        r_relatif = relative_forward_error(b,RHS,la);

        resvec[*nbite] = r_relatif;

        

        //b_exp update

        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1, AB, *lab, X, 1, 0, b, 1);

        

        

        *nbite = *nbite + 1;

    }



    free(b);

    free(ipiv);

    free(I_inv);

}

