//
//  outputcalc.c
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 4/3/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#include "outputcalc.h"

// function to calculate and write outputs along bottom wall
void getWallCoeffs(double **U, struct MeshData mesh, struct ResData resdata, struct JData jcb, struct Parameters bc, int bumplevel){
    int np = resdata.np;
    double state[4]; int k = 0;
    // allocate basis function array
    double *phi = malloc(resdata.np * sizeof(double));
    double xi = 0.0; double eta = 0.0;
    double dxi = 0.0; double deta = 0.0;
    
    // normal vector stuff
    double tangent[2]; double nvector[2]; double edgelen;
    
    // calculate free stream dynamic pressure
    double p,v;
    double Pfstream = 0.5*bc.gamma*bc.Pinf*bc.Minf*bc.Minf;
    
    // initialize variables
    double CL = 0.0;
    double CD = 0.0;
    double *CP = malloc(jcb.Ncurvelem*resdata.n1 * sizeof(double));
    double *xwall = malloc(jcb.Ncurvelem*resdata.n1 * sizeof(double));
    int nd1 = 0; int nd2 = 0; double xstart = 0.0;
    int kv = 0; int edge = 0;
    
    // loop through boundary edges
    for (kv = 0; kv<jcb.Ncurvelem; kv++){
        edge = jcb.bedgeNum[kv];
        k = jcb.elemNum[kv] - 1;
        // integrate along edge
        // loop over edge quadrature points
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
            // calculate basis functions at the quadrature point
            calcPhi(phi, resdata.p, xi, eta);
            
            // calculate normal vector
            tangent[0] = jcb.Jedge[kv*resdata.n1+q][0][0]*dxi + jcb.Jedge[kv*resdata.n1+q][0][1]*deta;
            tangent[1] = jcb.Jedge[kv*resdata.n1+q][1][0]*dxi + jcb.Jedge[kv*resdata.n1+q][1][1]*deta;
            edgelen = sqrt(tangent[0]*tangent[0] + tangent[1]*tangent[1]);
            nvector[0] =  tangent[1]/edgelen;
            nvector[1] = -tangent[0]/edgelen;
                
            // loop over basis functions at each quad point to get the state at quad point
            for (int var=0; var<4; var++){
                state[var] = 0.0;
            }
            for (int i=0; i<resdata.np; i++){
                // do for each state variable
                for (int var=0; var<4; var++){
                    state[var] += U[k*np+i][var] * phi[i];
                }
            }
            
            // calculate lift, drag coefficient contributions
            v = (state[1]/state[0])*(state[1]/state[0]) + (state[2]/state[0])*(state[2]/state[0]);
            p = (bc.gamma-1.0)*(state[3] - 0.5*state[0]*v);
            CL += (p-bc.Pinf)*nvector[1]*resdata.wq1[q]*edgelen / Pfstream / bc.h;
            CD += (p-bc.Pinf)*nvector[0]*resdata.wq1[q]*edgelen / Pfstream / bc.h;
            
            // calculate pressure coefficient (just for quad point)
            CP[kv*resdata.n1+q] = (p-bc.Pinf) / Pfstream;
            
            // calculate xcoordinate of this quadrature point
            // first get the nodes of the edge
            if (mesh.B2E[edge][1] == 1){
                nd1 = mesh.E[k][1]; nd2 = mesh.E[k][2];
            } else if (mesh.B2E[edge][1] == 2){
                nd1 = mesh.E[k][2]; nd2 = mesh.E[k][0];
            } else if (mesh.B2E[edge][1] == 3){
                nd1 = mesh.E[k][0]; nd2 = mesh.E[k][1];
            }
            // find left-most node xcoord
            xstart = mesh.V[nd1-1][0];
//            if (mesh.V[nd1-1][0] < mesh.V[nd2-1][0]){
//                xstart = mesh.V[nd1-1][0];
//            } else {
//                xstart = mesh.V[nd2-1][0];
//            }
            
            // now set the x coord of the quadrature point
            xwall[kv*resdata.n1+q] = xstart + resdata.xq1[q]*(mesh.V[nd2-1][0] - mesh.V[nd1-1][0]);// / mesh.BN[edge][2];
            
        } // end of loop over quadrature points
        //kv++;
    } // end of loop through boundary edges
    
    
    // write to file
    // CL
    char CLfile_str[100]; sprintf(CLfile_str,"outputs/CL_p%d_bump%d.txt",resdata.p,bumplevel);
    FILE *CLfile = fopen(CLfile_str,"w+");
    fprintf(CLfile,"%e\n",CL);
    fclose(CLfile);
    // CD
    char CDfile_str[100]; sprintf(CDfile_str,"outputs/CD_p%d_bump%d.txt",resdata.p,bumplevel);
    FILE *CDfile = fopen(CDfile_str,"w+");
    fprintf(CDfile, "%e\n",CD);
    fclose(CDfile);
    // CP
    char CPfile_str[100]; sprintf(CPfile_str,"outputs/CP_p%d_bump%d.txt",resdata.p,bumplevel);
    FILE *CPfile = fopen(CPfile_str,"w+");
    for (int i=0; i<(jcb.Ncurvelem*resdata.n1); i++){
        fprintf(CPfile,"%e %e\n",xwall[i],CP[i]);
    }
    fclose(CPfile);
    
    free(phi);
    free(CP);
}





// function to calculate and write entropy output for domain
void getEntropyError(double **U, struct MeshData mesh, struct ResData resdata, struct JData jcb, struct Parameters bc, int bumplevel){
    int np = resdata.np;
    double state[4];
    // allocate basis function array
    double *phi = malloc(resdata.np * sizeof(double));
    
    // calculate the stagnation entropy
    double rhot = bc.Ptot / bc.R / bc.Ttot;
    double st = bc.Ptot / pow(rhot, bc.gamma);
    
    // set entropy error to zero
    double Es = 0.0;
    // initialize entropy variables
    double s = 0.0; double p = 0.0; double v = 0.0;
    // jacobian
    double J = 0.0;
    
    // loop through elements
    for (int k=0; k<mesh.Nelem; k++){
        // integrate within the element
        // loop over interior quadrature points
        for (int q=0; q<resdata.n2; q++){
            // calculate basis functions at the quadrature point
            calcPhi(phi, resdata.p, resdata.xq2[q][0], resdata.xq2[q][1]);
            // initialize temporary state variable to zero
            state[0] = 0.0; state[1] = 0.0; state[2] = 0.0; state[3] = 0.0;
            // loop over basis functions at each quad point to get the state at quad point
            for (int i=0; i<resdata.np; i++){
                // do for each state variable
                for (int var=0; var<4; var++){
                    state[var] += U[k*np+i][var]*phi[i];
                }
            }
            // get jacobian
            for (int kv=0; kv<jcb.Ncurvelem; kv++){
                // if element is curved, get jacobian for this quad point
                if (jcb.elemNum[kv] == (k+1)){
                    J = jcb.detJint[kv*jcb.n2+q];
                    break;
                // otherwise, set jacobian to the linear triangle determinant
                } else{
                    J = mesh.detJ[k];
                }
            }
            
            // calculate entropy error for this element
            v = (state[1]/state[0])*(state[1]/state[0]) + (state[2]/state[0])*(state[2]/state[0]);
            p = (bc.gamma-1)*(state[3] - 0.5*state[0]*v);
            s = p / pow(state[0], bc.gamma);
            Es += (s/st-1.0)*(s/st-1.0)*resdata.wq2[q]*J / mesh.Area[k];
        }
    }
    // square root of sum
//    printf("Entropy error = %e\n",sqrt(Es));
    
    // write to file
    char Esfile_str[100]; sprintf(Esfile_str,"outputs/Es_p%d_bump%d.txt",resdata.p,bumplevel);
    FILE *Esfile = fopen(Esfile_str,"w+");
    fprintf(Esfile,"%e\n",sqrt(Es));
    fclose(Esfile);
    
    // free phi array
    free(phi);
}
