//
//  main.c
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 3/12/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "precompute.h"
#include "processmesh.h"
#include "calcRes.h"
#include "fluxfunctions.h"
#include "outputcalc.h"


int main(int argc, const char * argv[]) {
    // boundary conditions, variables
    struct Parameters bc;
    bc.gamma = 1.4;
    bc.Minf = 0.5;
    bc.alpha = 0.0;
    bc.R = 1.0;
    bc.Pinf = 1.0;
    bc.Ttot = 1.0 + bc.Minf*bc.Minf*(bc.gamma-1)/2;
    bc.Ptot = pow(bc.Ttot,bc.gamma/(bc.gamma-1));
    bc.h = 0.0625;
    bc.tolerance = pow(10,-7);
    bc.freestate[0] = 1.0;
    bc.freestate[1] = bc.Minf*cos(bc.alpha);
    bc.freestate[2] = bc.Minf*sin(bc.alpha);
    bc.freestate[3] = 1/(bc.gamma*(bc.gamma-1))+bc.Minf*bc.Minf/2;
    
    
    // read in mesh
    int bumplevel = 0;
    char grifile[100]; sprintf(grifile,"meshes/bump%d.gri",bumplevel);
    printf("Processing mesh...\n");
    struct MeshData mesh = processmesh(grifile);
    
    // get quadrature parameters
    int p = 2;
    struct ResData resdata = getResData(p);
    double CFL = 0.5/(p+1);
    
    
    // get jacobian matrices
    printf("Precomputing jacobians...\n");
    struct JData jcb = getJData(mesh, resdata);
    
    
    
    // get mass matrix & inverse
    double ***M = malloc(mesh.Nelem * sizeof(double **));
    double ***Minv = malloc(mesh.Nelem * sizeof(double **));
    for (int k=0; k<mesh.Nelem; k++){
        M[k] = malloc(resdata.np * sizeof(double *));
        Minv[k] = malloc(resdata.np * sizeof(double *));
        for (int i=0; i<resdata.np; i++){
            M[k][i] = calloc(resdata.np, sizeof(double *));
            Minv[k][i] = calloc(resdata.np, sizeof(double *));
        }
    }
    getMassMatrix(M, mesh, resdata, jcb);
    // test: print out mass matrix
//    for (int k=0; k<mesh.Nelem; k++){
//        for (int i=0; i<resdata.np; i++){
//            for (int j=0; j<resdata.np; j++){
//                printf("%e ",M[k][i][j]);
//            }
//            printf("\n");
//        }
//        printf("\n");
//    }
    
    
    getMMinverse(M, Minv, mesh, resdata);
    // test: print inverse mass matrix
//    for (int k=0; k<mesh.Nelem; k++){
//        for (int i=0; i<resdata.np; i++){
//            for (int j=0; j<resdata.np; j++){
//                printf("%e ",Minv[k][i][j]);
//            }
//            printf("\n");
//        }
//        printf("\n");
//    }
    
    
    // create state vector, and temps
    double **U = malloc(mesh.Nelem*resdata.np * sizeof(double *));
    double **Uinterm = malloc(mesh.Nelem*resdata.np * sizeof(double *));
    for (int i=0; i<(mesh.Nelem*resdata.np); i++){
        U[i] = (double *) calloc(4, sizeof(double));
        Uinterm[i] = (double *) calloc(4, sizeof(double));
    }
    
//    // initialize state to get coefficients without doing timestepping
//    // read from file
//    char statefname_convg[100]; sprintf(statefname_convg,"outputs/U_p%d_bump%d.txt",p,bumplevel);
//    FILE *statefile_convg = fopen(statefname_convg,"r");
//    for (int k=0; k<mesh.Nelem; k++){
//        for (int i=0; i<resdata.np; i++){
//            fscanf(statefile_convg, "%lf %lf %lf %lf\n", &U[k*resdata.np+i][0], &U[k*resdata.np+i][1], &U[k*resdata.np+i][2], &U[k*resdata.np+i][3]);
//        }
//    }
//    // get outputs
//    getWallCoeffs(U, mesh, resdata, jcb, bc, bumplevel);
//    getEntropyError(U, mesh, resdata, jcb, bc, bumplevel);
    
    
    
    // initialize the state
    if (p == 0){
        // set state to freestream
        for (int i=0; i<(mesh.Nelem*resdata.np); i++){
            U[i][0] = bc.freestate[0];
            U[i][1] = bc.freestate[1];
            U[i][2] = bc.freestate[2];
            U[i][3] = bc.freestate[3];
        }
    } else {
        // read from file
        char statefname[100]; sprintf(statefname,"outputs/U_p0_bump%d.txt",bumplevel);
        FILE *statefile = fopen(statefname,"r");
        for (int k=0; k<mesh.Nelem; k++){
            fscanf(statefile, "%lf %lf %lf %lf\n", &U[k*resdata.np][0], &U[k*resdata.np][1], &U[k*resdata.np][2], &U[k*resdata.np][3]);
            for (int i=1; i<resdata.np; i++){
                U[k*resdata.np+i][0] = U[k*resdata.np][0];
                U[k*resdata.np+i][1] = U[k*resdata.np][1];
                U[k*resdata.np+i][2] = U[k*resdata.np][2];
                U[k*resdata.np+i][3] = U[k*resdata.np][3];
            }
        }
        fclose(statefile);
    }
    
    // create residual vectors, zero out entries
    double **Res1 = malloc(mesh.Nelem*resdata.np * sizeof(double *));
    double **Res2 = malloc(mesh.Nelem*resdata.np * sizeof(double *));
    double **Res3 = malloc(mesh.Nelem*resdata.np * sizeof(double *));
    double **Res4 = malloc(mesh.Nelem*resdata.np * sizeof(double *));
    for (int i=0; i<(mesh.Nelem*resdata.np); i++){
        Res1[i] = (double *) calloc(4, sizeof(double));
        Res2[i] = (double *) calloc(4, sizeof(double));
        Res3[i] = (double *) calloc(4, sizeof(double));
        Res4[i] = (double *) calloc(4, sizeof(double));
    }
    // create wavespeed vector for each cell
    double *wavespeed = calloc(mesh.Nelem, sizeof(double));
    double *wavespeedFE = calloc(mesh.Nelem, sizeof(double));
    double *timestep = malloc(mesh.Nelem * sizeof(double));
    
    // track which stage of RK4
    int stage = 0;
    
    
    // OPEN OUTPUT FILES FOR PRINTING

    // write residual norm history to file
    char resnorm_str[100]; sprintf(resnorm_str,"outputs/resnorm_p%d_bump%d.txt",p,bumplevel);
//    char resnorm_str[100]; sprintf(resnorm_str,"outputs/resnorm_fstream_p%d_bump%d.txt",p,bumplevel);
    FILE *resnormfile = fopen(resnorm_str,"w+");

    
    //--------------------------------------------------------------------------------------------------------------
    // start timestepping loop here
    //--------------------------------------------------------------------------------------------------------------
    int t = 0;
    double resnorm = 1.0;
    while (resnorm > bc.tolerance){
        
        // FIRST STAGE UPDATE
        // calculate residual
        stage = 1;
        //printf("stage %d\n",stage);
        calcRes(Res1,wavespeed,U,mesh,resdata,jcb,bc);
        
        // get maximum residual
        resnorm = getLinfNorm(Res1, mesh.Nelem*resdata.np, 4);
        if (isnan(resnorm)){
            printf("---NaN---\n");
            break;
        }
        // write residual norm to file
        fprintf(resnormfile, "%.20lf\n", resnorm);
        // output progress to console
        printf("Iteration %d, |Res| = %e\n",t,resnorm);
        
        // set starting point U
        for (int i=0; i<mesh.Nelem*resdata.np; i++){
            for (int j=0; j<4; j++){
                Uinterm[i][j] = U[i][j];
            }
        }
        
        // update 2nd stage state
        for (int k=0; k<mesh.Nelem; k++){
            // calculate local timestep
            timestep[k] = 2*CFL * mesh.Area[k] / (wavespeed[k]*mesh.Per[k]);
            
            // get forward Euler state update
            // loop through basis functions i,j
            for (int i=0; i<resdata.np; i++){
                for (int j=0; j<resdata.np; j++){
                    // repeat for each state variable
                    for (int var=0; var<4; var++){
                        Uinterm[k*resdata.np+i][var] +=  - Minv[k][i][j]*timestep[k]*Res1[k*resdata.np+j][var]/2;
                    }
                }
            }
        }
        
        // SECOND STAGE UPDATE
        // calculate residual
        stage = 2;
//        printf("stage %d\n",stage);
        calcRes(Res2, wavespeedFE, Uinterm, mesh, resdata, jcb, bc);

        // set to starting point U
        for (int i=0; i<mesh.Nelem*resdata.np; i++){
            for (int j=0; j<4; j++){
                Uinterm[i][j] = U[i][j];
            }
        }

        // update 3rd stage state
        for (int k=0; k<mesh.Nelem; k++){
            // loop through basis functions i,j
            for (int i=0; i<resdata.np; i++){
                for (int j=0; j<resdata.np; j++){
                    // set U for next timestep using RK2 update
                    // repeat for each state variable
                    for (int var=0; var<4; var++){
                        Uinterm[k*resdata.np+i][var] +=  - Minv[k][i][j]*timestep[k]*Res2[k*resdata.np+j][var]/2;
                    }
                }
            }
        }
        
        // THIRD STAGE UPDATE
        // calculate residual
        stage = 3;
//        printf("stage %d\n",stage);
        calcRes(Res3, wavespeedFE, Uinterm, mesh, resdata, jcb, bc);
        
        // set to starting point U
        for (int i=0; i<mesh.Nelem*resdata.np; i++){
            for (int j=0; j<4; j++){
                Uinterm[i][j] = U[i][j];
            }
        }
        
        // update 3rd stage state
        for (int k=0; k<mesh.Nelem; k++){
            // loop through basis functions i,j
            for (int i=0; i<resdata.np; i++){
                for (int j=0; j<resdata.np; j++){
                    // set U for next timestep using RK2 update
                    // repeat for each state variable
                    for (int var=0; var<4; var++){
                        Uinterm[k*resdata.np+i][var] +=  - Minv[k][i][j]*timestep[k]*Res3[k*resdata.np+j][var];
                    }
                }
            }
        }
        
        // FOURTH STAGE UPDATE
        // calculate residual
        stage = 4;
//        printf("stage %d\n",stage);
        calcRes(Res4, wavespeedFE, Uinterm, mesh, resdata, jcb, bc);

        
        // FULL STATE UPDATE
        // set to starting point U
        for (int i=0; i<mesh.Nelem*resdata.np; i++){
            for (int j=0; j<4; j++){
                Uinterm[i][j] = U[i][j];
            }
        }
        
        // apply RK4 state update
        for (int k=0; k<mesh.Nelem; k++){
            // loop through basis functions i,j
            for (int i=0; i<resdata.np; i++){
                for (int j=0; j<resdata.np; j++){
                    // set U for next timestep using RK2 update
                    // repeat for each state variable
                    for (int var=0; var<4; var++){
                        Uinterm[k*resdata.np+i][var] +=  - Minv[k][i][j]*timestep[k]*(Res1[k*resdata.np+j][var] + 2*Res2[k*resdata.np+j][var] + 2*Res3[k*resdata.np+j][var] + Res4[k*resdata.np+j][var])/6;
                    }
                }
            }
        }
        
        // set U to the RK4 update
        for (int i=0; i<mesh.Nelem*resdata.np; i++){
            for (int j=0; j<4; j++){
                U[i][j] = Uinterm[i][j];
            }
        }
        
        // go back to top of loop, compute Res with new U
        t++;
    }
    fclose(resnormfile);
    //--------------------------------------------------------------------------------------------------------------
    // end timestepping loop
    //--------------------------------------------------------------------------------------------------------------
    
    

    
    
    //--------------------------------------------------------------------------------------------------------------
    // write state vector, outputs to file
    //--------------------------------------------------------------------------------------------------------------
    // write final U state to file
//    char u_str[100]; sprintf(u_str,"outputs/U_p%d_bump%d.txt",p,bumplevel);
////    char u_str[100]; sprintf(u_str,"kfid_U.txt");
//    FILE *ufile = fopen(u_str,"w+");
//    for (int i=0; i<mesh.Nelem*resdata.np; i++){
//        fprintf(ufile,"%.20lf %.20lf %.20lf %.20lf\n",U[i][0],U[i][1],U[i][2],U[i][3]);
//    }
//    fclose(ufile);
    
    
    // get outputs
    getWallCoeffs(U, mesh, resdata, jcb, bc, bumplevel);

    getEntropyError(U, mesh, resdata, jcb, bc, bumplevel);
    
    //--------------------------------------------------------------------------------------------------------------
    // free any allocated arrays
    //--------------------------------------------------------------------------------------------------------------
    // free mass matrix
    for (int k=0; k<mesh.Nelem; k++){
        for (int i=0; i<resdata.np; i++){
            free(M[k][i]);
            free(Minv[k][i]);
        }
        free(M[k]);
        free(Minv[k]);
    }
    free(M);
    free(Minv);
    
    
    
    // free residual vector
    for (int i=0; i<(mesh.Nelem*resdata.np); i++){
        free(Res1[i]);
        free(Res2[i]);
        free(Res3[i]);
        free(Res4[i]);
    }
    free(Res1);
    free(Res2);
    free(Res3);
    free(Res4);
    // free wavespeed vector
    free(wavespeed);
    free(wavespeedFE);
    // free timesteps
    free(timestep);

    
    // free state vector
    for (int i=0; i<(mesh.Nelem*resdata.np); i++){
        free(U[i]);
        free(Uinterm[i]);
    }
    free(U);
    free(Uinterm);
    
    
    // free quadrature parameter arrays
    freeResData(resdata);
    
    // free mesh data
    freeMesh(mesh);
    
    // free jacobian arrays
    freeJData(jcb);
    
    return 0;
}
