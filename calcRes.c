//
//  calcRes.c
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 3/14/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#include "calcRes.h"


void calcRes(double **Res, double *wavespeed, double **U, const struct MeshData mesh, const struct ResData resdata, struct JData jcb, const struct Parameters bc){
    int np = resdata.np; // for ease of writing indices later
    
    // zero out the residual
    for (int k=0; k<mesh.Nelem; k++){
        for (int j=0; j<resdata.np; j++){
            for (int var=0; var<4; var++){
                Res[k*np+j][var] = 0.0;
            }
        }
    }
    
    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------
    // loop over elements to calculate interior contributions
    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    double aflux[4][2]; // analytical flux vector
    double vel[2];      // velocity vector
    double P;           // pressure
    double J;           // determinant of Jacobian
    double Jinv[2][2];  // Jacobian inverse matrix
    int iscurved = 0;   // raise flag if element is curved
    int kv = 0;         // curved element indexing
    
    // basis function gradient
    double **gradPhi = malloc(resdata.np * sizeof(double *));
    for (int i=0; i<resdata.np; i++){
        gradPhi[i] = malloc(2 * sizeof(double));
    }
    // temporary variable to store conversion to global space
    double globalGPhi[2];
    
    // basis function
    double *phi = malloc(resdata.np * sizeof(double));
    // state u in reference space at one element
    double *Uref = malloc(4 * sizeof(double));
    
    // loop over elements
    for (int k=0; k<mesh.Nelem; k++){
        // calculate jacobian determinant and  inverse jacobian matrix
        // check if element is curved or not
        iscurved = 0;
        for (int i=0; i<jcb.Ncurvelem; i++){
            if (jcb.elemNum[i] == (k+1)){
                iscurved = 1;
                kv = i;
                break;
            }
        }
        
        
        // loop over quadrature points
        for (int q=0; q<resdata.n2; q++){
            // get Jacobian determinant and inverse at the quad point
            if (iscurved == 1){ // if this element is curved, use jacobian array
                J = jcb.detJint[kv*resdata.n2+q];
                Jinv[0][0] = jcb.Jintinv[kv*resdata.n2+q][0][0]; Jinv[0][1] = jcb.Jintinv[kv*resdata.n2+q][0][1];
                Jinv[1][0] = jcb.Jintinv[kv*resdata.n2+q][1][0]; Jinv[1][1] = jcb.Jintinv[kv*resdata.n2+q][1][1];
            } else{ // use linear jacobian
                J = mesh.detJ[k];
                Jinv[0][0] = mesh.Jinv[k][0][0]; Jinv[0][1] = mesh.Jinv[k][0][1];
                Jinv[1][0] = mesh.Jinv[k][1][0]; Jinv[1][1] = mesh.Jinv[k][1][1];
            }
            
            // get basis function gradient values at this quadrature point
            calcGradPhi(gradPhi, resdata.p, resdata.xq2[q][0], resdata.xq2[q][1]);
            
            // get state Uref at this element & quadrature point
            calcPhi(phi, resdata.p, resdata.xq2[q][0], resdata.xq2[q][1]);
            // initialize Uref to zero
            for (int i=0; i<4; i++){
                Uref[i] = 0.0;
            }
            // sum up phi*u for each basis function
            for (int j=0; j<resdata.np; j++){
                // repeat for each of 4 state variables
                for (int var=0; var<4; var++){
                    Uref[var] += phi[j] * U[k*resdata.np+j][var];
                }
            }
            
            // inner loop over basis gradient functions per element
            for (int i=0; i<resdata.np; i++){
                
                // calculate analytical flux of element (in reference space)
                vel[0] = Uref[1]/Uref[0]; vel[1] = Uref[2]/Uref[0];
                P = (bc.gamma-1)*(Uref[3]-0.5*Uref[0]*(vel[0]*vel[0]+vel[1]*vel[1]));
                aflux[0][0] = Uref[1];                                      aflux[0][1] = Uref[2];
                aflux[1][0] = Uref[1]*Uref[1]/Uref[0]+P;                    aflux[1][1] = Uref[1]*Uref[2]/Uref[0];
                aflux[2][0] = Uref[1]*Uref[2]/Uref[0];                      aflux[2][1] = Uref[2]*Uref[2]/Uref[0]+P;
                aflux[3][0] = Uref[1]*Uref[3]/Uref[0]+Uref[1]*P/Uref[0];    aflux[3][1] = Uref[2]*Uref[3]/Uref[0]+Uref[2]*P/Uref[0];
                
                // convert gradient from ref space to global space
                globalGPhi[0] = gradPhi[i][0]*Jinv[0][0] + gradPhi[i][1]*Jinv[1][0];
                globalGPhi[1] = gradPhi[i][0]*Jinv[0][1] + gradPhi[i][1]*Jinv[1][1];
                gradPhi[i][0] = globalGPhi[0];
                gradPhi[i][1] = globalGPhi[1];
                
                // calculate dot product integral and add to residual (minus sign)
                Res[k*np+i][0] -= (gradPhi[i][0]*aflux[0][0] + gradPhi[i][1]*aflux[0][1]) * J * resdata.wq2[q];
                Res[k*np+i][1] -= (gradPhi[i][0]*aflux[1][0] + gradPhi[i][1]*aflux[1][1]) * J * resdata.wq2[q];
                Res[k*np+i][2] -= (gradPhi[i][0]*aflux[2][0] + gradPhi[i][1]*aflux[2][1]) * J * resdata.wq2[q];
                Res[k*np+i][3] -= (gradPhi[i][0]*aflux[3][0] + gradPhi[i][1]*aflux[3][1]) * J * resdata.wq2[q];
            }
        }
//        printf("%e %e %e %e\n",Res[k][0],Res[k][1],Res[k][2],Res[k][3]);
    }
    // free gradPhi array
    for (int i=0; i<resdata.np; i++){
        free(gradPhi[i]);
    }
    free(gradPhi);
    // free phi array
    free(phi);
    free(Uref);
    
//    double interiornorm = getLinfNorm(Res, mesh.Nelem*resdata.np, 4);
//    printf("interior resnorm = %lf\n",interiornorm);
    
    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------
    // loop over edges to get interface contributions
    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------
    
    // left/right state vectors
    double *uL = calloc(4, sizeof(double));
    double *uR = calloc(4, sizeof(double));
    // numerical flux
    double *fhat = calloc(4, sizeof(double));
    // adjacent element numbers
    int elemL, elemR;
    // wavespeed check
    double wavespeed_new = 0.0;
    // basis functions
    double *phiL = malloc(resdata.np * sizeof(double));
    double *phiR = malloc(resdata.np * sizeof(double));
    // quadrature point coordinates in 2D reference space
    double xi = 0.0; double eta = 0.0;
    
    
    // TESTING FLUXES
//    double *nvec = malloc(2 * sizeof(double)); nvec[0] = 0.0; nvec[1] = -1.0;
//    double *state1 = malloc(4 * sizeof(double)); double *state2 = malloc(4 * sizeof(double));
//    state1[0] = 1.10554037; state1[1] = 0.52349369; state1[2] = 0.03032289; state1[3] = 2.17947479;
//    state2[0] = 1.02594646; state2[1] = 0.5046758 ; state2[2] = 0.01323541; state2[3] = 1.97556506;
////    double temp = roeflux(fhat, state1, state2, nvec, bc.gamma);
//    double temp = outflow(fhat, state1, nvec, bc);
//    for (int i=0; i<4; i++){
//        printf("%lf ",fhat[i]);
//    }
//    printf("\n"); printf("%e\n",temp);
    
    //--------------------------------------------------------------------------------------------------------------
    // INTERIOR EDGE CONTRIBUTIONS
    //--------------------------------------------------------------------------------------------------------------
    
    // interior edges loop
    for (int edge=0; edge<mesh.niedge; edge++){
        // adjacent element numbers
        elemL = mesh.I2E[edge][0];
        elemR = mesh.I2E[edge][2];
        
        // loop over 1D edge quadrature points
        for (int q=0; q<resdata.n1; q++){
            // calculate ref space coordinates for this quad point
            // ... for the left element
            if (mesh.I2E[edge][1] == 1){
                xi = 1 - resdata.xq1[q];
                eta = resdata.xq1[q];
            } else if (mesh.I2E[edge][1] == 2){
                xi = 0.0;
                eta = 1 - resdata.xq1[q];
            } else if (mesh.I2E[edge][1] == 3){
                xi = resdata.xq1[q];
                eta = 0.0;
            }
            // calculate left basis functions at quadrature point
            calcPhi(phiL, resdata.p, xi, eta);
            
            // ... for the right element
                if (mesh.I2E[edge][3] == 1){
                    xi = resdata.xq1[q];
                    eta = 1 - resdata.xq1[q];
                } else if (mesh.I2E[edge][3] == 2){
                    xi = 0.0;
                    eta = resdata.xq1[q];
                } else if (mesh.I2E[edge][3] == 3){
                    xi = 1 - resdata.xq1[q];
                    eta = 0.0;
                }
            // calculate right basis functions at quadrature points
            calcPhi(phiR, resdata.p, xi, eta);
            
            // calculate cell states (L,R) to pass into numerical flux function (Roe)
            // intialize sum to zero
            uL[0] = 0.0; uL[1] = 0.0; uL[2] = 0.0; uL[3] = 0.0;
            uR[0] = 0.0; uR[1] = 0.0; uR[2] = 0.0; uR[3] = 0.0;
            // sum up phi*u for each basis function
            for (int j=0; j<resdata.np; j++){
                // repeat for each of 4 state variables
                for (int var=0; var<4; var++){
                    uL[var] += phiL[j] * U[(elemL-1)*resdata.np+j][var];
                    uR[var] += phiR[j] * U[(elemR-1)*resdata.np+j][var];
                }
            }
            
            // calculate numerical flux using edge cell states
            wavespeed_new = roeflux(fhat, uL, uR, mesh.IN[edge], bc.gamma);
            // check maximum wavespeed
            if (wavespeed_new > wavespeed[elemL-1]){
                wavespeed[elemL-1] = wavespeed_new;
            }
            if (wavespeed_new > wavespeed[elemR-1]){
                wavespeed[elemR-1] = wavespeed_new;
            }
            
            // update residual with basis functions times fhat
            for (int i=0; i<resdata.np; i++){
                // do all state variables
                for (int var=0; var<4; var++){
                    Res[(elemL-1)*np+i][var] += phiL[i]*( fhat[var])*mesh.IN[edge][2]*resdata.wq1[q];
                    Res[(elemR-1)*np+i][var] += phiR[i]*(-fhat[var])*mesh.IN[edge][2]*resdata.wq1[q];
                }
            }
            
        }// end of loop over quadrature points
    }// end of interior edges loop

    
    
    //--------------------------------------------------------------------------------------------------------------
    // BOUNDARY EDGE CONTRIBUTIONS
    //--------------------------------------------------------------------------------------------------------------

    int freestreamtest = 0;// set to zero to run actual boundary conditions

    // start boundary edges loop
    double tangent[2] = {0,0};
    double edgelen = 0;
    double *nvector = malloc(2 * sizeof(double));
    double dxi = 0; double deta = 0;
    kv = 0; // curved edge counter
    for (int edge=0; edge<mesh.nbedge; edge++){
        // interior element number
        elemL = mesh.B2E[edge][0];
        iscurved = 0;
        // check if the edge is curved
        if (mesh.B2E[edge][2] == 4){
            iscurved = 1;
        }
        
        // loop over 1D edge quadrature points
        for (int q=0; q<resdata.n1; q++){
            // calculate ref space coordinates for this quadrature point
            if (mesh.B2E[edge][1] == 1){
                xi = 1 - resdata.xq1[q];
                eta = resdata.xq1[q];
                dxi = -1; deta = 1;
            } else if (mesh.B2E[edge][1] == 2){
                xi = 0.0;
                eta = 1 - resdata.xq1[q];
                dxi = 0; deta = -1;
            } else if (mesh.B2E[edge][1] == 3){
                xi = resdata.xq1[q];
                eta = 0.0;
                dxi = 1; deta = 0.0;
            }
            // get normal
            if (iscurved){//} && (resdata.Q>1)){
                tangent[0] = jcb.Jedge[kv*resdata.n1+q][0][0]*dxi + jcb.Jedge[kv*resdata.n1+q][0][1]*deta;
                tangent[1] = jcb.Jedge[kv*resdata.n1+q][1][0]*dxi + jcb.Jedge[kv*resdata.n1+q][1][1]*deta;
                edgelen = sqrt(tangent[0]*tangent[0] + tangent[1]*tangent[1]);
                nvector[0] =  tangent[1]/edgelen;
                nvector[1] = -tangent[0]/edgelen;
//                printf("%lf %lf %lf \n",nvector[0],nvector[1],sqrt(nvector[0]*nvector[0] + nvector[1]*nvector[1]));
                
            } else{
                nvector[0] = mesh.BN[edge][0];
                nvector[1] = mesh.BN[edge][1];
                edgelen = mesh.BN[edge][2];
            }
            
            // calculate basis functions at quadrature point
            calcPhi(phiL, resdata.p, xi, eta);
            // calculate cell state on edge at quad point
            // intialize sum to zero
            uL[0] = 0.0; uL[1] = 0.0; uL[2] = 0.0; uL[3] = 0.0;
            // sum up phi*u for each basis function
            for (int j=0; j<resdata.np; j++){
                // repeat for each of 4 state variables
                for (int var=0; var<4; var++){
                    uL[var] += phiL[j] * U[(elemL-1)*resdata.np+j][var];
                }
            }
            
            // check which boundary conditions to apply for numerical flux
            if (freestreamtest != 0){ // freestream test, ignores other boundary conditions
                // set "right" state equal to freestream
                uR[0] = bc.freestate[0];
                uR[1] = bc.freestate[1];
                uR[2] = bc.freestate[2];
                uR[3] = bc.freestate[3];
                wavespeed_new = roeflux(fhat, uL, uR, nvector, bc.gamma);
                
            } else { // actual boundary conditions
                if (mesh.B2E[edge][2] == 1){ // left boundary
                    // apply subsonic inflow conditions
                    wavespeed_new = inflow(fhat, uL, nvector, bc);
                } else if (mesh.B2E[edge][2] == 2){ // right boundary
                    // apply subsonic outflow conditions
                    wavespeed_new = outflow(fhat, uL, nvector, bc);
                } else if ((mesh.B2E[edge][2] == 3) || (mesh.B2E[edge][2] == 4)){ // boundaries 3,4 (top, bottom)
                    // apply inviscid wall conditions
                    wavespeed_new = inviscidWall(fhat, uL, nvector, bc);
                }
            }
            // check wavespeed
            if (wavespeed_new > wavespeed[elemL-1]){
                wavespeed[elemL-1] = wavespeed_new;
            }
            
            // update residual with basis functions times fhat
            for (int i=0; i<resdata.np; i++){
                // do all state variables
                for (int var=0; var<4; var++){
                    Res[(elemL-1)*np+i][var] += phiL[i]*fhat[var]*edgelen*resdata.wq1[q];
                }
            }
        } // end of quadrature points loop
        // if edge was curved, increment curve edge counter
        if (mesh.B2E[edge][2] == 4){
            kv++;
        }
    } // end of boundary edges loop
    
//    printf("bEdge kv = %d\n",kv);
    
    free(phiL);free(phiR);
    // free left/right state vectors
    free(uL); free(uR);
    // free numerical flux
    free(fhat);
    // free edge normal
    free(nvector);
 
    // END OF CALCRES FUNCTION
}




void calcPhi(double *phi, int p, double xi, double eta){
    // pass in pointer to basis function values
    if (p == 0){
        phi[0] = 1.0;
    } else if (p == 1){
        phi[0] = 1 - xi - eta;
        phi[1] = xi;
        phi[2] = eta;
    } else if (p == 2){
        phi[0] = 1.0 - 3*xi + 2*xi*xi - 3*eta + 4*xi*eta + 2*eta*eta;
        phi[1] = 4*xi - 4*xi*xi - 4*xi*eta;
        phi[2] = -xi + 2*xi*xi;
        phi[3] = 4*eta - 4*xi*eta - 4*eta*eta;
        phi[4] = 4*xi*eta;
        phi[5] = -eta + 2*eta*eta;
    } else if (p == 3){
        // 1    x    x^2    x^3    y   xy    x^2y    y^2    xy^2    y^3
        phi[0] = 1.0 - 5.5*xi + 9*xi*xi     - 4.5*xi*xi*xi  - 5.5*eta    + 18*xi*eta    - 13.5*xi*xi*eta  + 9*eta*eta  - 13.5*xi*eta*eta - 4.5*eta*eta*eta;
        phi[1] =         9*xi - 22.5*xi*xi  - 13.5*xi*xi*xi              - 22.5*xi*eta  + 27*xi*xi*eta                 + 13.5*xi*eta*eta;
        phi[2] =      -4.5*xi + 18*xi*xi    - 13.5*xi*xi*xi              + 4.5*xi*eta   - 13.5*xi*xi*eta;
        phi[3] =       1.0*xi - 4.5*xi*xi   + 4.5*xi*xi*xi;
        phi[4] =                                                9*eta   - 22.5*xi*eta + 13.5*xi*xi*eta  - 22.5*eta*eta + 27*xi*eta*eta  + 13.5*eta*eta*eta;
        phi[5] =                                                            27*xi*eta -   27*xi*xi*eta                 - 27*xi*eta*eta;
        phi[6] =                                                          -4.5*xi*eta + 13.5*xi*xi*eta;
        phi[7] =                                             -4.5*eta    + 4.5*xi*eta                     + 18*eta*eta - 13.5*xi*eta*eta - 13.5*eta*eta*eta;
        phi[8] =                                                          -4.5*xi*eta + 13.5*xi*eta*eta;
        phi[9] =                                              1.0*eta                                    - 4.5*eta*eta                   + 4.5*eta*eta*eta;
    }
}



void calcGradPhi(double **gradPhi, int p, double xi, double eta){
    // pass in pointer to gradient of basis function values
    if (p == 0){
        gradPhi[0][0] = 0.0;                    gradPhi[0][1] = 0.0;
    } else if (p == 1){
        gradPhi[0][0] = -1.0;                   gradPhi[0][1] = -1.0;
        gradPhi[1][0] = 1.0;                    gradPhi[1][1] = 0.0;
        gradPhi[2][0] = 0.0;                    gradPhi[2][1] = 1.0;
    } else if (p == 2){
        gradPhi[0][0] = -3 + 4*xi + 4*eta;      gradPhi[0][1] = -3 + 4*xi + 4*eta;
        gradPhi[1][0] = 4 - 8*xi - 4*eta;       gradPhi[1][1] = -4*xi;
        gradPhi[2][0] = -1 + 4*xi;              gradPhi[2][1] = 0.0;
        gradPhi[3][0] = -4*eta;                 gradPhi[3][1] = 4 - 4*xi - 8*eta;
        gradPhi[4][0] = 4*eta;                  gradPhi[4][1] = 4*xi;
        gradPhi[5][0] = 0.0;                    gradPhi[5][1] = -1 + 4*eta;
    } else if (p == 3){
        gradPhi[0][0] = -5.5 + 2*9*xi - 3*4.5*xi*xi + 18*eta - 2*13.5*xi*eta - 13.5*eta*eta;    gradPhi[0][1] = -5.5 + 18*xi - 13.5*xi*xi + 2*9*eta - 2*13.5*xi*eta - 3*4.5*eta*eta;
        gradPhi[1][0] = 9 - 2*22.5*xi + 3*13.5*xi*xi - 22.5*eta + 2*27*xi*eta + 13.5*eta*eta;   gradPhi[1][1] = -22.5*xi + 27*xi*xi + 2*13.5*xi*eta;
        gradPhi[2][0] = -4.5 + 2*18*xi - 3*13.5*xi*xi + 4.5*eta - 2*13.5*xi*eta;                gradPhi[2][1] = 4.5*xi - 13.5*xi*xi;
        gradPhi[3][0] = 1 - 2*4.5*xi + 3*4.5*xi*xi;                                             gradPhi[3][1] = 0.0;
        gradPhi[4][0] = -22.5*eta + 2*13.5*xi*eta + 27*eta*eta;                                 gradPhi[4][1] = 9.0-22.5*xi+13.5*xi*xi - 2*22.5*eta + 2*27*xi*eta + 3*13.5*eta*eta;
        gradPhi[5][0] = 27*eta - 2*27*xi*eta - 27*eta*eta;                                      gradPhi[5][1] = 27*xi - 27*xi*xi - 2*27*xi*eta;
        gradPhi[6][0] = -4.5*eta + 2*13.5*xi*eta;                                               gradPhi[6][1] = -4.5*xi + 13.5*xi*xi;
        gradPhi[7][0] = 4.5*eta - 13.5*eta*eta;                                                 gradPhi[7][1] = -4.5 + 4.5*xi - 2*18*eta - 2*13.5*xi*eta - 3*13.5*eta*eta;
        gradPhi[8][0] = -4.5*eta + 13.5*eta*eta;                                                gradPhi[8][1] = -4.5*xi + 2*13.5*xi*eta;
        gradPhi[9][0] = 0.0;                                                                    gradPhi[9][1] = 1 - 2*4.5*eta + 3*4.5*eta*eta;
    }
}


