//
//  precompute.c
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 3/19/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#include "precompute.h"

struct ResData getResData(int p){
    
    // setup quadrature data
    struct ResData resdata;
    
    // approximation order
    resdata.p = p;
    resdata.np = (p+1)*(p+2)/2;
    // geometry order
    resdata.Q = 2;//p+1;
    resdata.nQ = (resdata.Q+1)*(resdata.Q+2)/2;
    
    if (p == 0){ // Order 1 Gauss-Legendre, Dunavant Points
        // 1D quadrature
        resdata.n1 = 1;
        resdata.xq1 = malloc(resdata.n1 * sizeof(double));
        resdata.xq1[0] = 0.500000000000000;
        resdata.wq1 = malloc(resdata.n1 * sizeof(double));
        resdata.wq1[0] = 1.000000000000000;
        // 2D quadrature
        resdata.n2 = 1;
        resdata.xq2 = malloc(resdata.n2 * sizeof(double *));
        resdata.xq2[0] = malloc(2 * sizeof(double));
        resdata.xq2[0][0] = 0.333333333333333; resdata.xq2[0][1] = 0.333333333333333;
        resdata.wq2 = malloc(resdata.n2 * sizeof(double));
        resdata.wq2[0] = 0.500000000000000;
        
//    } else if (p == 1){ // Order 3 Gauss-Legendre, Dunavant Points
//        // 1D quadrature
//        resdata.n1 = 2;
//        resdata.xq1 = malloc(resdata.n1 * sizeof(double));
//        resdata.xq1[0] = 0.211324865405187;
//        resdata.xq1[1] = 0.788675134594813;
//        resdata.wq1 = malloc(resdata.n1 * sizeof(double));
//        resdata.wq1[0] = 0.5;
//        resdata.wq1[1] = 0.5;
//        // 2D quadrature
//        resdata.n2 = 4;
//        resdata.xq2 = malloc(resdata.n2 * sizeof(double *));
//        for (int i=0; i<resdata.n2; i++){
//            resdata.xq2[i] = malloc(2 * sizeof(double));
//        }
//        resdata.xq2[0][0] = 0.333333333333333; resdata.xq2[0][1] = 0.333333333333333;
//        resdata.xq2[1][0] = 0.600000000000000; resdata.xq2[1][1] = 0.200000000000000;
//        resdata.xq2[2][0] = 0.200000000000000; resdata.xq2[2][1] = 0.200000000000000;
//        resdata.xq2[3][0] = 0.200000000000000; resdata.xq2[3][1] = 0.600000000000000;
//        resdata.wq2 = malloc(resdata.n2 * sizeof(double));
//        resdata.wq2[0] = -0.281250000000000;
//        resdata.wq2[1] = 0.260416666666667;
//        resdata.wq2[2] = 0.260416666666667;
//        resdata.wq2[3] = 0.260416666666667;
        
        
    } else if (p == 1){ // Order 5 Gauss-Legendre, Dunavant Points
        // 1D quadrature
        resdata.n1 = 3;
        resdata.xq1 = malloc(resdata.n1 * sizeof(double));
        resdata.xq1[0] = 0.112701665379258;
        resdata.xq1[1] = 0.500000000000000;
        resdata.xq1[2] = 0.887298334620742;
        resdata.wq1 = malloc(resdata.n1 * sizeof(double));
        resdata.wq1[0] = 0.277777777777778;
        resdata.wq1[1] = 0.444444444444444;
        resdata.wq1[2] = 0.277777777777778;
        // 2D quadrature
        resdata.n2 = 7;
        resdata.xq2 = malloc(resdata.n2 * sizeof(double *));
        for (int i=0; i<resdata.n2; i++){
            resdata.xq2[i] = malloc(2 * sizeof(double));
        }
        resdata.xq2[0][0] = 0.333333333333333; resdata.xq2[0][1] = 0.333333333333333;
        resdata.xq2[1][0] = 0.059715871789770; resdata.xq2[1][1] = 0.470142064105115;
        resdata.xq2[2][0] = 0.470142064105115; resdata.xq2[2][1] = 0.470142064105115;
        resdata.xq2[3][0] = 0.470142064105115; resdata.xq2[3][1] = 0.059715871789770;
        resdata.xq2[4][0] = 0.797426985353087; resdata.xq2[4][1] = 0.101286507323456;
        resdata.xq2[5][0] = 0.101286507323456; resdata.xq2[5][1] = 0.101286507323456;
        resdata.xq2[6][0] = 0.101286507323456; resdata.xq2[6][1] = 0.797426985353087;
        resdata.wq2 = malloc(resdata.n2 * sizeof(double));
        resdata.wq2[0] = 0.112500000000000;
        resdata.wq2[1] = 0.066197076394253;
        resdata.wq2[2] = 0.066197076394253;
        resdata.wq2[3] = 0.066197076394253;
        resdata.wq2[4] = 0.062969590272414;
        resdata.wq2[5] = 0.062969590272414;
        resdata.wq2[6] = 0.062969590272414;
        
    } else if (p == 2){
        // 1d quadrature
        resdata.n1 = 5;
        resdata.xq1 = malloc(resdata.n1 * sizeof(double));
        resdata.xq1[0] = 0.046910077030668;
        resdata.xq1[1] = 0.230765344947158;
        resdata.xq1[2] = 0.500000000000000;
        resdata.xq1[3] = 0.769234655052841;
        resdata.xq1[4] = 0.953089922969332;
        resdata.wq1 = malloc(resdata.n1 * sizeof(double));
        resdata.wq1[0] = 0.118463442528095;
        resdata.wq1[1] = 0.239314335249683;
        resdata.wq1[2] = 0.284444444444444;
        resdata.wq1[3] = 0.239314335249683;
        resdata.wq1[4] = 0.118463442528095;
        // 2D quadrature
        resdata.n2 = 19;
        resdata.xq2 = malloc(resdata.n2 * sizeof(double *));
        for (int i=0; i<resdata.n2; i++){
            resdata.xq2[i] = malloc(2 * sizeof(double));
        }
        resdata.xq2[0][0] = 0.333333333333333; resdata.xq2[0][1] = 0.333333333333333;
        resdata.xq2[1][0] = 0.020634961602525; resdata.xq2[1][1] = 0.489682519198738;
        resdata.xq2[2][0] = 0.489682519198738; resdata.xq2[2][1] = 0.489682519198738;
        resdata.xq2[3][0] = 0.489682519198738; resdata.xq2[3][1] = 0.020634961602525;
        resdata.xq2[4][0] = 0.125820817014127; resdata.xq2[4][1] = 0.437089591492937;
        resdata.xq2[5][0] = 0.437089591492937; resdata.xq2[5][1] = 0.437089591492937;
        resdata.xq2[6][0] = 0.437089591492937; resdata.xq2[6][1] = 0.125820817014127;
        resdata.xq2[7][0] = 0.623592928761935; resdata.xq2[7][1] = 0.188203535619033;
        resdata.xq2[8][0] = 0.188203535619033; resdata.xq2[8][1] = 0.188203535619033;
        resdata.xq2[9][0] = 0.188203535619033; resdata.xq2[9][1] = 0.623592928761935;
        resdata.xq2[10][0] = 0.910540973211095; resdata.xq2[10][1] = 0.044729513394453;
        resdata.xq2[11][0] = 0.044729513394453; resdata.xq2[11][1] = 0.044729513394453;
        resdata.xq2[12][0] = 0.044729513394453; resdata.xq2[12][1] = 0.910540973211095;
        resdata.xq2[13][0] = 0.036838412054736; resdata.xq2[13][1] = 0.221962989160766;
        resdata.xq2[14][0] = 0.221962989160766; resdata.xq2[14][1] = 0.741198598784498;
        resdata.xq2[15][0] = 0.741198598784498; resdata.xq2[15][1] = 0.036838412054736;
        resdata.xq2[16][0] = 0.221962989160766; resdata.xq2[16][1] = 0.036838412054736;
        resdata.xq2[17][0] = 0.741198598784498; resdata.xq2[17][1] = 0.221962989160766;
        resdata.xq2[18][0] = 0.036838412054736; resdata.xq2[18][1] = 0.741198598784498;
        resdata.wq2 = malloc(resdata.n2 * sizeof(double));
        resdata.wq2[0] = 0.048567898141400;
        resdata.wq2[1] = 0.015667350113570;
        resdata.wq2[2] = 0.015667350113570;
        resdata.wq2[3] = 0.015667350113570;
        resdata.wq2[4] = 0.038913770502387;
        resdata.wq2[5] = 0.038913770502387;
        resdata.wq2[6] = 0.038913770502387;
        resdata.wq2[7] = 0.039823869463605;
        resdata.wq2[8] = 0.039823869463605;
        resdata.wq2[9] = 0.039823869463605;
        resdata.wq2[10] = 0.012788837829349;
        resdata.wq2[11] = 0.012788837829349;
        resdata.wq2[12] = 0.012788837829349;
        resdata.wq2[13] = 0.021641769688645;
        resdata.wq2[14] = 0.021641769688645;
        resdata.wq2[15] = 0.021641769688645;
        resdata.wq2[16] = 0.021641769688645;
        resdata.wq2[17] = 0.021641769688645;
        resdata.wq2[18] = 0.021641769688645;
    } else{
        printf("%d is not a supported p value\n",p);
    }
    
    
    return resdata;
}



void freeResData(struct ResData resdata){
    free(resdata.xq1);
    for (int i=0; i<resdata.n2; i++){
        free(resdata.xq2[i]);
    }
    free(resdata.xq2);
    free(resdata.wq1);
    free(resdata.wq2);
}





struct JData getJData(struct MeshData mesh, struct ResData resdata){
    // precomputes the Jacobian at each quadrature point in the curved elements
    int n1 = resdata.n1; int n2 = resdata.n2;
    
    // initialize the struct
    struct JData jcb;
    
    jcb.n1 = resdata.n1;
    jcb.n2 = resdata.n2;
    
    // interior: identify how many elements are curved
    int Ncurvelem = 0;
    for (int edge=0; edge<mesh.nbedge; edge++){
        if (mesh.B2E[edge][2] == 4){
            Ncurvelem++;
        }
    }
    jcb.Ncurvelem = Ncurvelem;
    
    // allocate the Jacobians for curved elements
    jcb.Jint = malloc(Ncurvelem*jcb.n2 * sizeof(double **));
    jcb.Jintinv = malloc(Ncurvelem*jcb.n2 * sizeof(double **));
    jcb.detJint = malloc(Ncurvelem*jcb.n2 * sizeof(double));
    for (int i=0; i<Ncurvelem*jcb.n2; i++){
        jcb.Jint[i] = malloc(2 * sizeof(double *));
        jcb.Jintinv[i] = malloc(2 * sizeof(double *));
        for (int j=0; j<2; j++){
            jcb.Jint[i][j] = malloc(2 * sizeof(double));
            jcb.Jintinv[i][j] = malloc(2 * sizeof(double));
        }
    }
    // allocate Jacobians for curved edges
    jcb.Jedge = malloc(Ncurvelem*jcb.n1 * sizeof(double **));
    jcb.Jedgeinv = malloc(Ncurvelem*jcb.n1 * sizeof(double **));
    jcb.detJedge = calloc(Ncurvelem*jcb.n1, sizeof(double));
    for (int i=0; i<Ncurvelem*jcb.n1; i++){
        jcb.Jedge[i] = malloc(2 * sizeof(double *));
        jcb.Jedgeinv[i] = malloc(2 * sizeof(double *));
        for (int j=0; j<2; j++){
            jcb.Jedge[i][j] = malloc(2 * sizeof(double));
            jcb.Jedgeinv[i][j] = malloc(2 * sizeof(double));
        }
    }
    
    // variables for evaluating quadrature points, geometry node coordinates
    double xi = 0.0; double eta = 0.0;
    int k = 0; // index of the element/edge that is curved
    int nd1 = 0; int nd2 = 0; int nd3 = 0; // node indices
    int elem = 0; // element number in the overall mesh
    double **gradPhi = malloc(resdata.nQ * sizeof(double *));
    double **XQ = malloc(resdata.nQ * sizeof(double *)); // geometry order node coords
    for (int i=0; i<resdata.nQ; i++){
        gradPhi[i] = malloc(2 * sizeof(double));
        XQ[i] = malloc(2 * sizeof(double));
    }

    
    // loop through boundary edges
    for (int edge=0; edge<mesh.nbedge; edge++){
        // check edge boundary type
        if (mesh.B2E[edge][2] == 4){
            // if this edge is on the bottom boundary, do calculations for Jacobian
            elem = mesh.B2E[edge][0];
            nd1 = mesh.E[elem-1][0]-1; nd2 = mesh.E[elem-1][1]-1; nd3 = mesh.E[elem-1][2]-1;
            // get the coordinates of geometry nodes
            
            if (resdata.Q == 1){
                XQ[0][0] = mesh.V[nd1][0]; XQ[0][1] = mesh.V[nd1][1];
                XQ[1][0] = mesh.V[nd2][0]; XQ[1][1] = mesh.V[nd2][1];
                XQ[2][0] = mesh.V[nd3][0]; XQ[2][1] = mesh.V[nd3][1];
            }
            // LOCAL EDGE 1
            else if (mesh.B2E[edge][1] == 1){
                if (resdata.Q == 2){
                    // six nodes total
                    // original nodes
                    XQ[0][0] = mesh.V[nd1][0]; XQ[0][1] = mesh.V[nd1][1];
                    XQ[2][0] = mesh.V[nd2][0]; XQ[2][1] = mesh.V[nd2][1];
                    XQ[5][0] = mesh.V[nd3][0]; XQ[5][1] = mesh.V[nd3][1];
                    // new nodes - unchanged
                    XQ[3][0] = 0.5*(mesh.V[nd1][0]+mesh.V[nd3][0]); XQ[3][1] = 0.5*(mesh.V[nd1][1]+mesh.V[nd3][1]);
                    XQ[4][0] = 0.5*(mesh.V[nd2][0]+mesh.V[nd3][0]); XQ[4][1] = exactBumpY(XQ[4][0]);
                    // new nodes - adjusted
                    XQ[1][0] = 0.5*(mesh.V[nd1][0]+mesh.V[nd2][0]); XQ[1][1] = 0.5*(mesh.V[nd1][1]+mesh.V[nd2][1]);
                } else if (resdata.Q == 3){
                    // ten nodes total
                    // original nodes
                    XQ[0][0] = mesh.V[nd1][0]; XQ[0][1] = mesh.V[nd1][1];
                    XQ[3][0] = mesh.V[nd2][0]; XQ[3][1] = mesh.V[nd2][1];
                    XQ[9][0] = mesh.V[nd3][0]; XQ[9][1] = mesh.V[nd3][1];
                    // new nodes - uncurved
                    XQ[4][0] = mesh.V[nd1][0]+(1/3)*(mesh.V[nd3][0]-mesh.V[nd1][0]); XQ[4][1] = mesh.V[nd1][1]+(1/3)*(mesh.V[nd3][1]-mesh.V[nd1][1]);
                    XQ[7][0] = mesh.V[nd1][0]+(2/3)*(mesh.V[nd3][0]-mesh.V[nd1][0]); XQ[7][1] = mesh.V[nd1][1]+(2/3)*(mesh.V[nd3][1]-mesh.V[nd1][1]);
                    XQ[6][0] = mesh.V[nd2][0]+(1/3)*(mesh.V[nd3][0]-mesh.V[nd2][0]); XQ[6][1] = exactBumpY(XQ[6][0]);
                    XQ[8][0] = mesh.V[nd2][0]+(2/3)*(mesh.V[nd3][0]-mesh.V[nd2][0]); XQ[8][1] = exactBumpY(XQ[8][0]);
                    // new nodes - adjusted
                    XQ[1][0] = mesh.V[nd1][0]+(1/3)*(mesh.V[nd2][0]-mesh.V[nd1][0]); XQ[1][1] = mesh.V[nd1][1]+(1/3)*(mesh.V[nd2][1]-mesh.V[nd1][1]);
                    XQ[2][0] = mesh.V[nd1][0]+(2/3)*(mesh.V[nd2][0]-mesh.V[nd1][0]); XQ[2][1] = mesh.V[nd1][1]+(2/3)*(mesh.V[nd2][1]-mesh.V[nd1][1]);
                    // centroid
                    XQ[5][0] = (1/3)*(mesh.V[nd1][0]+mesh.V[nd2][0]+mesh.V[nd3][0]); XQ[5][1] = (1/3)*(mesh.V[nd1][1]+mesh.V[nd2][1]+mesh.V[nd3][1]);
                }
                
            // LOCAL EDGE 2
            } else if (mesh.B2E[edge][1] == 2){
                if (resdata.Q == 2){
                    // six nodes total
                    // original nodes
                    XQ[0][0] = mesh.V[nd1][0]; XQ[0][1] = mesh.V[nd1][1];
                    XQ[2][0] = mesh.V[nd2][0]; XQ[2][1] = mesh.V[nd2][1];
                    XQ[5][0] = mesh.V[nd3][0]; XQ[5][1] = mesh.V[nd3][1];
                    // new nodes - unchanged
                    XQ[3][0] = 0.5*(mesh.V[nd1][0]+mesh.V[nd3][0]); XQ[3][1] = exactBumpY(XQ[3][0]);
                    XQ[4][0] = 0.5*(mesh.V[nd2][0]+mesh.V[nd3][0]); XQ[4][1] = 0.5*(mesh.V[nd2][1]+mesh.V[nd3][1]);
                    // new nodes - adjusted
                    XQ[1][0] = 0.5*(mesh.V[nd1][0]+mesh.V[nd2][0]); XQ[1][1] = 0.5*(mesh.V[nd1][1]+mesh.V[nd2][1]);
                } else if (resdata.Q == 3){
                    // ten nodes total
                    // original nodes
                    XQ[0][0] = mesh.V[nd1][0]; XQ[0][1] = mesh.V[nd1][1];
                    XQ[3][0] = mesh.V[nd2][0]; XQ[3][1] = mesh.V[nd2][1];
                    XQ[9][0] = mesh.V[nd3][0]; XQ[9][1] = mesh.V[nd3][1];
                    // new nodes - uncurved
                    XQ[4][0] = mesh.V[nd1][0]+(1/3)*(mesh.V[nd3][0]-mesh.V[nd1][0]); XQ[4][1] = exactBumpY(XQ[4][0]);
                    XQ[7][0] = mesh.V[nd1][0]+(2/3)*(mesh.V[nd3][0]-mesh.V[nd1][0]); XQ[7][1] = exactBumpY(XQ[7][0]);
                    XQ[6][0] = mesh.V[nd2][0]+(1/3)*(mesh.V[nd3][0]-mesh.V[nd2][0]); XQ[6][1] = mesh.V[nd2][1]+(1/3)*(mesh.V[nd3][1]-mesh.V[nd2][1]);
                    XQ[8][0] = mesh.V[nd2][0]+(2/3)*(mesh.V[nd3][0]-mesh.V[nd2][0]); XQ[8][1] = mesh.V[nd2][1]+(2/3)*(mesh.V[nd3][1]-mesh.V[nd2][1]);
                    // new nodes - adjusted
                    XQ[1][0] = mesh.V[nd1][0]+(1/3)*(mesh.V[nd2][0]-mesh.V[nd1][0]); XQ[1][1] = mesh.V[nd1][1]+(1/3)*(mesh.V[nd2][1]-mesh.V[nd1][1]);
                    XQ[2][0] = mesh.V[nd1][0]+(2/3)*(mesh.V[nd2][0]-mesh.V[nd1][0]); XQ[2][1] = mesh.V[nd1][1]+(2/3)*(mesh.V[nd2][1]-mesh.V[nd1][1]);
                    // centroid
                    XQ[5][0] = (1/3)*(mesh.V[nd1][0]+mesh.V[nd2][0]+mesh.V[nd3][0]); XQ[5][1] = (1/3)*(mesh.V[nd1][1]+mesh.V[nd2][1]+mesh.V[nd3][1]);
                }
                
            // LOCAL EDGE 3
            } else if (mesh.B2E[edge][1] == 3){
                if (resdata.Q == 2){
                    // six nodes total
                    // original nodes
                    XQ[0][0] = mesh.V[nd1][0]; XQ[0][1] = mesh.V[nd1][1];
                    XQ[2][0] = mesh.V[nd2][0]; XQ[2][1] = mesh.V[nd2][1];
                    XQ[5][0] = mesh.V[nd3][0]; XQ[5][1] = mesh.V[nd3][1];
                    // new nodes - uncurved
                    XQ[3][0] = 0.5*(mesh.V[nd1][0]+mesh.V[nd3][0]); XQ[3][1] = 0.5*(mesh.V[nd1][1]+mesh.V[nd3][1]);
                    XQ[4][0] = 0.5*(mesh.V[nd2][0]+mesh.V[nd3][0]); XQ[4][1] = 0.5*(mesh.V[nd2][1]+mesh.V[nd3][1]);
                    // new nodes - adjusted
                    XQ[1][0] = 0.5*(mesh.V[nd1][0]+mesh.V[nd2][0]); XQ[1][1] = exactBumpY(XQ[1][0]);
                } else if (resdata.Q == 3){
                    // ten nodes total
                    // original nodes
                    XQ[0][0] = mesh.V[nd1][0]; XQ[0][1] = mesh.V[nd1][1];
                    XQ[3][0] = mesh.V[nd2][0]; XQ[3][1] = mesh.V[nd2][1];
                    XQ[9][0] = mesh.V[nd3][0]; XQ[9][1] = mesh.V[nd3][1];
                    // new nodes - uncurved
                    XQ[4][0] = mesh.V[nd1][0]+(1/3)*(mesh.V[nd3][0]-mesh.V[nd1][0]); XQ[4][1] = mesh.V[nd1][1]+(1/3)*(mesh.V[nd3][1]-mesh.V[nd1][1]);
                    XQ[7][0] = mesh.V[nd1][0]+(2/3)*(mesh.V[nd3][0]-mesh.V[nd1][0]); XQ[7][1] = mesh.V[nd1][1]+(2/3)*(mesh.V[nd3][1]-mesh.V[nd1][1]);
                    XQ[6][0] = mesh.V[nd2][0]+(1/3)*(mesh.V[nd3][0]-mesh.V[nd2][0]); XQ[6][1] = mesh.V[nd2][1]+(1/3)*(mesh.V[nd3][1]-mesh.V[nd2][1]);
                    XQ[8][0] = mesh.V[nd2][0]+(2/3)*(mesh.V[nd3][0]-mesh.V[nd2][0]); XQ[8][1] = mesh.V[nd2][1]+(2/3)*(mesh.V[nd3][1]-mesh.V[nd2][1]);
                    // new nodes - adjusted
                    XQ[1][0] = mesh.V[nd1][0]+(1/3)*(mesh.V[nd2][0]-mesh.V[nd1][0]); XQ[1][1] = exactBumpY(XQ[1][0]);
                    XQ[2][0] = mesh.V[nd1][0]+(2/3)*(mesh.V[nd2][0]-mesh.V[nd1][0]); XQ[2][1] = exactBumpY(XQ[2][0]);
                    // centroid
                    XQ[5][0] = (1/3)*(mesh.V[nd1][0]+mesh.V[nd2][0]+mesh.V[nd3][0]); XQ[5][1] = (1/3)*(mesh.V[nd1][1]+mesh.V[nd2][1]+mesh.V[nd3][1]);
                }
            }
            
            
            
            // loop over interior quadrature points
            for (int q=0; q<resdata.n2; q++){
                // calculate gradPhi at this quadrature point
                calcGradPhi(gradPhi, resdata.Q, resdata.xq2[q][0], resdata.xq2[q][1]);
                // zero out jacobian entries
                jcb.Jint[k*n2+q][0][0] = 0.0; jcb.Jint[k*n2+q][0][1] = 0.0;
                jcb.Jint[k*n2+q][1][0] = 0.0; jcb.Jint[k*n2+q][1][1] = 0.0;
                // loop through basis functions to calculate Jacobian
                for (int i=0; i<resdata.nQ; i++){
                    jcb.Jint[k*n2+q][0][0] += XQ[i][0]*gradPhi[i][0]; jcb.Jint[k*n2+q][0][1] += XQ[i][0]*gradPhi[i][1];
                    jcb.Jint[k*n2+q][1][0] += XQ[i][1]*gradPhi[i][0]; jcb.Jint[k*n2+q][1][1] += XQ[i][1]*gradPhi[i][1];
                }
                // calculate the determinant of the jacobian
                jcb.detJint[k*n2+q] = (jcb.Jint[k*n2+q][0][0]*jcb.Jint[k*n2+q][1][1]) - (jcb.Jint[k*n2+q][0][1]*jcb.Jint[k*n2+q][1][0]);
                // inverse jacobian matrix
                jcb.Jintinv[k*n2+q][0][0] =   jcb.Jint[k*n2+q][1][1]/jcb.detJint[k*n2+q]; jcb.Jintinv[k*n2+q][0][1] = - jcb.Jint[k*n2+q][0][1]/jcb.detJint[k*n2+q];
                jcb.Jintinv[k*n2+q][1][0] = - jcb.Jint[k*n2+q][1][0]/jcb.detJint[k*n2+q]; jcb.Jintinv[k*n2+q][1][1] =   jcb.Jint[k*n2+q][0][0]/jcb.detJint[k*n2+q];
            } // end of interior quadrature point loop
            
            // compute Jacobian of edge quadrature points
            for (int q=0; q<resdata.n1; q++){
                // calculate ref space coordinates for this quadrature point
                if (mesh.B2E[edge][1] == 1){
                    xi = 1 - resdata.xq1[q];
                    eta = resdata.xq1[q];
                } else if (mesh.B2E[edge][1] == 2){
                    xi = 0.0;
                    eta = 1 - resdata.xq1[q];
                } else if (mesh.B2E[edge][1] == 3){
                    xi = resdata.xq1[q];
                    eta = 0.0;
                }
                // calculate gradPhi at this quadrature point
                calcGradPhi(gradPhi, resdata.Q, xi, eta);
                // zero out jacobian entries
                jcb.Jedge[k*n1+q][0][0] = 0.0; jcb.Jedge[k*n1+q][0][1] = 0.0;
                jcb.Jedge[k*n1+q][1][0] = 0.0; jcb.Jedge[k*n1+q][1][1] = 0.0;
                // loop through basis functions to calculate jacobian
                for (int i=0; i<resdata.nQ; i++){
                    jcb.Jedge[k*n1+q][0][0] += XQ[i][0]*gradPhi[i][0]; jcb.Jedge[k*n1+q][0][1] += XQ[i][0]*gradPhi[i][1];
                    jcb.Jedge[k*n1+q][1][0] += XQ[i][1]*gradPhi[i][0]; jcb.Jedge[k*n1+q][1][1] += XQ[i][1]*gradPhi[i][1];
                }
                // also calculate determinant
                jcb.detJedge[k*n1+q] = (jcb.Jedge[k*n1+q][0][0]*jcb.Jedge[k*n1+q][1][1]) - (jcb.Jedge[k*n1+q][0][1]*jcb.Jedge[k*n1+q][1][0]);
                // also calculate Jinverse
                jcb.Jedgeinv[k*n1+q][0][0] =   jcb.Jedge[k*n1+q][1][1]/jcb.detJedge[k*n1+q]; jcb.Jedgeinv[k*n1+q][0][1] = - jcb.Jedge[k*n1+q][0][1]/jcb.detJedge[k*n1+q];
                jcb.Jedgeinv[k*n1+q][1][0] = - jcb.Jedge[k*n1+q][1][0]/jcb.detJedge[k*n1+q]; jcb.Jedgeinv[k*n1+q][1][1] =   jcb.Jedge[k*n1+q][0][0]/jcb.detJedge[k*n1+q];
            } // end of edge quadrature point loop
            
            // increment curved element/edge counter
//            printf("%d\n",k);
            k++;
        } // END OF IF CHECK BTYPE == 4
        
    } // END OF BOUNDARY EDGES LOOP
    
    // create references to the right elements/edges
    jcb.elemNum = malloc((k+1) * sizeof(int));
    jcb.bedgeNum = malloc((k+1) * sizeof(int));
    k = 0;
    for (int j=0; j<mesh.nbedge; j++){
        if (mesh.B2E[j][2] == 4){
            jcb.elemNum[k] = mesh.B2E[j][0];
            jcb.bedgeNum[k] = j;
            k++;
        }
    }
    
    // free gradPhi array
    for (int i=0; i<resdata.nQ; i++){
        free(gradPhi[i]);
        free(XQ[i]);
    }
    free(gradPhi);
    free(XQ);
    
    // return the completed Jacobian matrices
    return jcb;
}



void freeJData(struct JData jcb){
    // free the Jacobians for curved elements
    for (int i=0; i<jcb.Ncurvelem*jcb.n2; i++){
        for (int j=0; j<2; j++){
            free(jcb.Jint[i][j]);
            free(jcb.Jintinv[i][j]);
        }
        free(jcb.Jint[i]);
        free(jcb.Jintinv[i]);
    }
    free(jcb.Jint);
    free(jcb.Jintinv);
    free(jcb.detJint);
    
    // free Jacobians for curved edges

    for (int i=0; i<jcb.Ncurvelem*jcb.n1; i++){
        for (int j=0; j<2; j++){
            free(jcb.Jedge[i][j]);
            free(jcb.Jedgeinv[i][j]);
        }
        free(jcb.Jedge[i]);
        free(jcb.Jedgeinv[i]);
    }
    free(jcb.Jedge);
    free(jcb.Jedgeinv);
    free(jcb.detJedge);
    
    // free element and bedge numbers
    free(jcb.elemNum);
    free(jcb.bedgeNum);
}




double exactBumpY(double x){
    double y = 0.0625*exp(-25*x*x);
    return y;
}




double getLinfNorm(double **Res, int nrows, int ncols){
    // evaluates the Linfinity Norm of the residual array
    double maxval = 0.0;
    // loop through rows
    for (int i=0; i<nrows; i++){
        // loop through each entry in row
        for (int j=0; j<ncols; j++){
            if (fabs(Res[i][j]) > maxval){
                maxval = fabs(Res[i][j]);
            }
        }
    }
    return maxval;
}



void getMassMatrix(double ***M, struct MeshData mesh, struct ResData resdata, struct JData jcb){
    // M should be initialized to zeros with calloc prior to being passed into this function
    
    // basis function
    double *phi = malloc(resdata.np * sizeof(double));
    
    int iscurved = 0; // raise to 1 if the element is curved
    
    // loop through elements (linear only)
    for (int k=0; k<mesh.Nelem; k++){
        
        // if the solution requires curved elements:
        if (resdata.p > 0){
            // check if the element is curved or not
            for (int a=0; a<jcb.Ncurvelem; a++){
                if (jcb.elemNum[a] == (k+1)){
                    iscurved = 1;
                }
            }
            // skip the linear jacobian on this element
            if (iscurved == 1){
                iscurved = 0;
                continue;
            }
        }
        
        // loop over each quadrature point
        for (int q=0; q<resdata.n2; q++){
            // calculate basis functions at this quadrature point
            calcPhi(phi, resdata.p, resdata.xq2[q][0], resdata.xq2[q][1]);
            
            // loop through basis functions i,j
            for (int i=0; i<resdata.np; i++){
                for (int j=0; j<resdata.np; j++){
                    M[k][i][j] += phi[i] * mesh.detJ[k] * resdata.wq2[q] * phi[j];
                }
            } // end of loop over basis functions
        } // end of loop over quadrature points
    } // end of loop over elements
    
    // loop through curved elements
    int elem = 0;
    for (int a=0; a<jcb.Ncurvelem; a++){
        elem = jcb.elemNum[a];
        // loop over each quadrature point
        for (int q=0; q<resdata.n2; q++){
            // calculate basis functions at this quadrature point
            calcPhi(phi, resdata.p, resdata.xq2[q][0], resdata.xq2[q][1]);
            
            // loop through basis functions i,j
            for (int i=0; i<resdata.np; i++){
                for (int j=0; j<resdata.np; j++){
                    M[elem-1][i][j] += phi[i] * jcb.detJint[a*resdata.n2+q] * resdata.wq2[q] * phi[j];
                }
            } // end of loop over basis functions
        } // end of loop over quadrature points
    } // end of loop over curved elements
    
    free(phi);
}


void getMMinverse(double ***M, double ***Minv, struct MeshData mesh, struct ResData resdata){
    // entire matrix is NxN = (Nelem*np)x(Nelem*np)
    // element-wise blocks M[k] are (np)x(np)
    double detM;
    for (int k=0; k<mesh.Nelem; k++){
        detM = Determinant(M[k], resdata.np);
        CoFactor(M[k], resdata.np, Minv[k]);
        Transpose(Minv[k], resdata.np);
        for (int i=0; i<resdata.np; i++){
            for (int j=0; j<resdata.np; j++){
                Minv[k][i][j] = Minv[k][i][j] / detM;
            }
        }
    }
}

/*
 Recursive definition of determinate using expansion by minors.
 */
double Determinant(double **a,int n){
    int i,j,j1,j2;
    double det = 0;
    double **m = NULL;
    
    if (n < 1) { /* Error */
        det = 1.0;
    } else if (n == 1) { /* Shouldn't get used */
        det = a[0][0];
    } else if (n == 2) {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {
        det = 0;
        for (j1=0;j1<n;j1++) {
            m = malloc((n-1)*sizeof(double *));
            for (i=0;i<n-1;i++)
                m[i] = malloc((n-1)*sizeof(double));
            for (i=1;i<n;i++) {
                j2 = 0;
                for (j=0;j<n;j++) {
                    if (j == j1)
                        continue;
                    m[i-1][j2] = a[i][j];
                    j2++;
                }
            }
            det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
            for (i=0;i<n-1;i++)
                free(m[i]);
            free(m);
        }
    }
    return(det);
}

/*
 Find the cofactor matrix of a square matrix
 */
void CoFactor(double **a,int n,double **b){
    int i,j,ii,jj,i1,j1;
    double det;
    double **c;
    
    c = malloc((n-1)*sizeof(double *));
    for (i=0;i<n-1;i++)
        c[i] = malloc((n-1)*sizeof(double));
    
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) {
            
            /* Form the adjoint a_ij */
            i1 = 0;
            for (ii=0;ii<n;ii++) {
                if (ii == i)
                    continue;
                j1 = 0;
                for (jj=0;jj<n;jj++) {
                    if (jj == j)
                        continue;
                    c[i1][j1] = a[ii][jj];
                    j1++;
                }
                i1++;
            }
            
            /* Calculate the determinate */
            det = Determinant(c,n-1);
            
            /* Fill in the elements of the cofactor */
            b[i][j] = pow(-1.0,i+j+2.0) * det;
        }
    }
    for (i=0;i<n-1;i++)
        free(c[i]);
    free(c);
}

/*
 Transpose of a square matrix, do it in place
 */
void Transpose(double **a,int n){
    int i,j;
    double tmp;
    
    for (i=1;i<n;i++) {
        for (j=0;j<i;j++) {
            tmp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = tmp;
        }
    }
}


