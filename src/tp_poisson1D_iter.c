/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include <time.h>

#define ALPHA 0
#define JAC 1
#define GS 2

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
    
    remove("PERF_GS_time"); //files names for perf
    remove("PERF_JAC_time");
    remove("PERF_ALPHA_time");
    
    int maxit = 500;
    int IMPLEM = 0;
    
    //variation de taille de matrice
    for(int n = 10; n < 500 ; n=n+50){
        
        struct timespec t1_a, t2_a, t1_jac, t2_jac , t1_gs, t2_gs;
        double elapsed = 0.0;
        
        //nombre d'iteration
        for(int it = 0; it < maxit ; ++it){
            
            int ierr;
            int jj;
            int nbpoints, la;
            int ku, kl, lab, kv;
            int *ipiv;
            int info;
            int NRHS;
            double T0, T1;
            double *RHS, *SOL, *EX_SOL, *X;
            double *AB;
            double *MB;
            
            double temp, relres;
            
            double opt_alpha;
            
            if (argc == 2) {
                IMPLEM = atoi(argv[1]);
            } else if (argc > 2) {
                perror("Application takes at most one argument");
                exit(1);
            }
            
            /* Size of the problem */
            NRHS=1;
            nbpoints=n;
            la=nbpoints-2;
            
            /* Dirichlet Boundary conditions */
            T0=5.0;
            T1=20.0;
            
            //printf("--------- Poisson 1D ---------\n\n");
            RHS=(double *) malloc(sizeof(double)*la);
            SOL=(double *) calloc(la, sizeof(double));
            EX_SOL=(double *) malloc(sizeof(double)*la);
            X=(double *) malloc(sizeof(double)*la);
            
            /* Setup the Poisson 1D problem */
            /* General Band Storage */
            set_grid_points_1D(X, &la);
            set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
            set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
            
            //write_vec(RHS, &la, "RHS.dat");
            //write_vec(EX_SOL, &la, "EX_SOL.dat");
            //write_vec(X, &la, "X_grid.dat");
            
            kv=0;
            ku=1;
            kl=1;
            lab=kv+kl+ku+1;
            
            AB = (double *) malloc(sizeof(double)*lab*la);
            set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
            
            /* uncomment the following to check matrix A */
            //write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
            
            /********************************************/
            /* Solution (Richardson with optimal alpha) */
            
            /* Computation of optimum alpha */
            opt_alpha = richardson_alpha_opt(&la);
            //printf("Optimal alpha for simple Richardson iteration is : %lf",opt_alpha);
            
            /* Solve */
            double tol=1e-3;
            int maxit=1000;
            double *resvec;
            int nbite=0;
            
            resvec=(double *) calloc(maxit, sizeof(double));
            
            /* Solve with Richardson alpha */
            if (IMPLEM == ALPHA) {
                clock_gettime(CLOCK_MONOTONIC_RAW, &t1_a);
                richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
                clock_gettime(CLOCK_MONOTONIC_RAW, &t2_a);
            }
            
            /* Richardson General Tridiag */
            
            /* get MB (:=M^-1, D^-1 for Jacobi, (D-E)^-1 for Gauss-seidel) */
            kv = 1;
            ku = 1;
            kl = 1;
            
            /* Solve with General Richardson */
            if (IMPLEM == JAC) {
                MB = (double *) malloc(sizeof(double)*(lab)*la);
                clock_gettime(CLOCK_MONOTONIC_RAW, &t1_jac);
                extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
                //write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "MB.dat");
                richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
                clock_gettime(CLOCK_MONOTONIC_RAW, &t2_jac);
                
            } else if (IMPLEM == GS) {
                MB = (double *) malloc(sizeof(double)*(la)*la);
                clock_gettime(CLOCK_MONOTONIC_RAW, &t1_gs);
                extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
                //write_GB_operator_colMajor_poisson1D(MB, &la, &la, "MB.dat");
                kl = la-1;
                ku = 0;
                richardson_MB(AB, RHS, SOL, MB, &la, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
                clock_gettime(CLOCK_MONOTONIC_RAW, &t2_gs);
            }
            
            /* Write solution */
            //write_vec(SOL, &la, "SOL.dat");
            
            /* Write convergence history */
            //write_vec(resvec, &nbite, "RESVEC.dat");
            
            /* Erreur par rapport a la sol analytique
             double e = relative_forward_error(SOL,EX_SOL,&la);
             printf("\nl’erreur par rapport a la solution analytique = %e\n",e);
             */
            
            free(resvec);
            free(RHS);
            free(SOL);
            free(EX_SOL);
            free(X);
            free(AB);
            free(MB);
            //printf("\n\n--------- End -----------\n");
        }
        
        if (IMPLEM == ALPHA) {
            
            FILE *file;
            file = fopen("PERF_ALPHA_time", "a");
            
            if (file == NULL) {
                fprintf(stderr, "Error opening the file.\n");
                return 1;
            }
            elapsed = (double)(t2_a.tv_nsec - t1_a.tv_nsec)/maxit;
            fprintf(file, "%d %f\n", n,elapsed);
            fclose(file);
        }
        
        
        
        if (IMPLEM == JAC) {
            
            FILE *file;
            file = fopen("PERF_JAC_time", "a");
            
            if (file == NULL) {
                fprintf(stderr, "Error opening the file.\n");
                return 1;
            }
            elapsed = (double)(t2_jac.tv_nsec - t1_jac.tv_nsec)/maxit;
            fprintf(file, "%d %f\n", n,elapsed);
            fclose(file);
        }
        
        
        if (IMPLEM == GS) {
            
            FILE *file;
            file = fopen("PERF_GS_time", "a");
            
            if (file == NULL) {
                fprintf(stderr, "Error opening the file.\n");
                return 1;
            }
            elapsed = (double)(t2_gs.tv_nsec - t1_gs.tv_nsec)/maxit;
            fprintf(file, "%d %f\n", n,elapsed);
            fclose(file);
        }
        
    }
}
