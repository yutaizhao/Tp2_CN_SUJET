/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include <string.h>
#include <time.h>

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */

{

remove("PERF_TRFnTRS"); //files names for perf
remove("PERF_SV");


  for(int n = 10; n < 1000 ; ++n){
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB, *AB_init;

  double relres;
  
  struct timespec t1_tr, t2_tr, t1_sv, t2_sv;
  double elapsed = 0.0;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]); //0 if no entry
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS=1;
  nbpoints=n; //10 initialement
  la=nbpoints-2; //colonne
  T0=-5.0;
  T1=5.0;

  //printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  /*
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");
  */

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1; //ligne

  AB = (double *) malloc(sizeof(double)*lab*la);
  
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  //write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  
  

  //printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));

  /* LU Factorization */
  if (IMPLEM == TRF) {
    clock_gettime(CLOCK_MONOTONIC_RAW, &t1_tr);
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    if (info==0){
        /* Solution (Triangular) */
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info,1);
      clock_gettime(CLOCK_MONOTONIC_RAW, &t2_tr);
      if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
      printf("\n INFO = %d\n",info);
    }
  }

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  if (IMPLEM == TRI) {
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    if (info==0){
        /* Solution (Triangular) */
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info,1);
      
      if (info!=0){
      
      		printf("\n INFO DGBTRS = %d\n",info);
      		
      }else{
          if(nbpoints == 10){
          	
          	/*Validation dgbtrftridiag*/
      		double r = relative_forward_error(RHS, EX_SOL, &la); //RHS is now the sol x
      		printf("\nEX6 : dgbtrftridiag validation pour n =10: r = error entre Sol calculee a partir de dgbtrftridiag et sol exacte \n");
  		printf("r = %e\n",r);
  	}
      }
    }else{
      printf("\n INFO = %d\n",info);
    }
  }


  /* It can also be solved with dgbsv */
  if (IMPLEM == SV) {
    clock_gettime(CLOCK_MONOTONIC_RAW, &t1_sv);
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    clock_gettime(CLOCK_MONOTONIC_RAW, &t2_sv);
  }

  //write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  //write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  if(nbpoints == 10){
  relres = relative_forward_error(RHS, EX_SOL, &la); //RHS is now the sol x
  printf("\nThe relative forward error is relres = %e\n",relres);
  }
  
  /* dgbmv validation */
  if(nbpoints == 10){
  AB_init = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D(AB_init, &lab, &la, &kv);
  
  double *EXP_RHS=(double *) malloc(sizeof(double)*la);
  memset(EXP_RHS, 0, sizeof(double) * la);
  double *EX_RHS=(double *) malloc(sizeof(double)*la);
  memset(EX_RHS, 0, sizeof(double) * la);
  
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1, AB_init, lab-1, RHS, 1, 1, EXP_RHS,1); //A*x = y => AB*RHS = EXP_RHS
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1, AB_init, lab-1, EX_SOL, 1, 1, EX_RHS,1); //A*x = y => AB*EX_SOL = EX_RHS
  
  double r = relative_forward_error(EXP_RHS, EX_RHS, &la); //RHS is now the sol x
  
  printf("\nEX4 : dgbmv validation pour n =10: r = erreur entre AB*X(calcule) et AB*X(exacte) \n");
  printf("r = %e\n",r);
  
  free(EXP_RHS);
  free(EX_RHS);
  free(AB_init);
  }
  
  
  if (IMPLEM == TRF) {
    
    FILE *file;
    file = fopen("PERF_TRFnTRS", "a");
    
    if (file == NULL) {
        fprintf(stderr, "Error opening the file.\n");
        return 1; 
    }
     elapsed = (double)(t2_tr.tv_nsec - t1_tr.tv_nsec);
     fprintf(file, "%d %f\n", n,elapsed);
     fclose(file);
  }
  
  
  
  if (IMPLEM == SV) {
    
    FILE *file;
    file = fopen("PERF_SV", "a");
    
    if (file == NULL) {
        fprintf(stderr, "Error opening the file.\n");
        return 1; 
    }
      elapsed = (double)(t2_sv.tv_nsec - t1_sv.tv_nsec);
      fprintf(file, "%d %f\n", n,elapsed);
     fclose(file);
  }
  
  
  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  //printf("\n\n--------- End -----------\n");
  }
}
