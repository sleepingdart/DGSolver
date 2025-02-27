//
//  processmesh.c
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 3/12/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "processmesh.h"


struct MeshData processmesh(char *grifile){
    
    FILE *meshfile = fopen(grifile,"r");
    if(meshfile == 0) {
        perror("fopen");
        exit(1);
    }
    
    // read number of nodes, elements, (dimensions)
    int Nnodes, Nelem, dim;
    fscanf(meshfile,"%d %d %d",&Nnodes,&Nelem,&dim);
    // allocate V matrix
    double **V;
    V = (double **) malloc(Nnodes * sizeof(double));
    for (int i=0; i<Nnodes; i++){
        V[i] = (double *) malloc(2 * sizeof(double));
    }
    
    // read node coordinates
    float x,y;
    for (int i=0; i<Nnodes; i++){
        fscanf(meshfile,"%f %f",&x,&y);
        V[i][0] = x;
        V[i][1] = y;
    }
    
    // read number of boundary groups
    int Nbgroups;
    fscanf(meshfile,"%d",&Nbgroups);
    int Nbedges[Nbgroups];
    
    int ***B; // top level Boundaries pointer
    // allocate memory to top level pointer
    B = (int ***) malloc(Nbgroups * sizeof(int **));
    
    // read boundary edges
    for (int btype=0; btype<Nbgroups; btype++){
        // read num of edges in this group
        fscanf(meshfile,"%d %*d %*s",&Nbedges[btype]);
        
        // dyn allocate rows for this btype
        B[btype] = (int **) malloc(Nbedges[btype] * sizeof(int *));
        // dyn allocate cols for each row
        for (int i=0; i<Nbedges[btype]; i++){
            B[btype][i] = (int *) malloc(2 * sizeof(int));
        }
        
        // read node numbers defining each boundary edge
        int nd1, nd2;
        for (int j=0; j<Nbedges[btype]; j++){
            fscanf(meshfile,"%d %d",&nd1,&nd2);
            B[btype][j][0] = nd1;
            B[btype][j][1] = nd2;
        }
    }
    
    // allocate E matrix
    int **E;
    E = (int **) malloc(Nelem * sizeof(int *));
    for (int i=0; i<Nelem; i++){
        E[i] = (int *) malloc(3 * sizeof(int));
    }
    
    // read elements
    fscanf(meshfile,"%*d %*d %*s");
    for (int i=0; i<Nelem; i++){
        fscanf(meshfile,"%d %d %d",&E[i][0],&E[i][1],&E[i][2]);
    }
    
    // close file
    fclose(meshfile);
    
    
    /* Do analysis */
    // allocate memory to connectivity sparse matrix
    int **H;
    H = (int **) malloc(Nelem * sizeof(int *));
    for (int i=0; i<Nelem; i++){
        H[i] = (int *) calloc(Nelem, sizeof(int));
    }
    // allocate memory to I2E matrix (upper bound guess for size)
    int **I2E;
    I2E = (int **) malloc(Nelem*3 * sizeof(int *));
    for (int i=0; i<Nelem*3; i++){
        I2E[i] = (int *) malloc(4 * sizeof(int));
    }
    // allocate memory to B2E matrix (upper bound guess for size)
    int **B2E;
    B2E = (int **) malloc(Nelem*3 * sizeof(int *));
    for (int i=0; i<Nelem*3; i++){
        B2E[i] = (int *) malloc(3 * sizeof(int));
    }
    
    // identify edges through sparse connectivity table
    int niedge = 0;
    int nmin,nmax,oldelem;
    for (int elem=0; elem<Nelem; elem++){
        for (int edge=0; edge<3; edge++){
            
            if ( (E[elem][(edge+1)%3]-1) < (E[elem][(edge+2)%3]-1) ){
                nmin = E[elem][(edge+1)%3]-1;
                nmax = E[elem][(edge+2)%3]-1;
            } else {
                nmin = E[elem][(edge+2)%3]-1;
                nmax = E[elem][(edge+1)%3]-1;
            }
            
            if (H[nmax][nmin] == 0){
                // edge hit for first time
                H[nmin][nmax] = elem+1;
                H[nmax][nmin] = edge+1;
            } else {
                // edge hit for second time
                oldelem = H[nmin][nmax];
                if ((elem+1)<oldelem){
                    I2E[niedge][0] = elem+1;
                    I2E[niedge][1] = edge+1;
                    I2E[niedge][2] = oldelem;
                    I2E[niedge][3] = H[nmax][nmin];
                } else {
                    I2E[niedge][0] = oldelem;
                    I2E[niedge][1] = H[nmax][nmin];
                    I2E[niedge][2] = elem+1;
                    I2E[niedge][3] = edge+1;
                }
                H[nmin][nmax] = 0;
                H[nmax][nmin] = -1;
                niedge += 1;
            }
        }
    }
    // free unused rows of I2E
    for (int i=niedge; i<(Nelem*3); i++){
        free(I2E[i]);
    }
    
    // count total number of boundary edges
    int sumNbedges = 0;
    for (int i=0; i<Nbgroups; i++){
        sumNbedges += Nbedges[i];
    }
    
    // construct B2E matrix
    int nbedge = 0;
    int elemindex = 0;
    for (int btype=0; btype<Nbgroups; btype++){
        for (int edge=0; edge<Nbedges[btype]; edge++){
            // find element index of this edge
            for (int i=0; i<Nelem; i++){
                if ( (E[i][0]==B[btype][edge][0])||(E[i][1]==B[btype][edge][0])||(E[i][2]==B[btype][edge][0]) ){
                    if ( (E[i][0]==B[btype][edge][1])||(E[i][1]==B[btype][edge][1])||(E[i][2]==B[btype][edge][1]) ){
                        elemindex = i;
                    }
                }
            }
            B2E[nbedge][0] = elemindex+1;
            // find local edge index of this edge
            if ( ((E[elemindex][0]==B[btype][edge][0])&&(E[elemindex][1]==B[btype][edge][1])) || ((E[elemindex][0]==B[btype][edge][1])&&(E[elemindex][1]==B[btype][edge][0])) ){
                B2E[nbedge][1] = 3;
            } else if( ((E[elemindex][0]==B[btype][edge][0])&&(E[elemindex][2]==B[btype][edge][1])) || ((E[elemindex][0]==B[btype][edge][1])&&(E[elemindex][2]==B[btype][edge][0])) ) {
                B2E[nbedge][1] = 2;
            } else if( ((E[elemindex][2]==B[btype][edge][0])&&(E[elemindex][1]==B[btype][edge][1])) || ((E[elemindex][2]==B[btype][edge][1])&&(E[elemindex][1]==B[btype][edge][0])) ){
                B2E[nbedge][1] = 1;
            }
            // boundary type
            B2E[nbedge][2] = btype+1;
            
            nbedge += 1;
        }
    }
    // free unused rows of B2E
    for (int i=nbedge; i<(Nelem*3); i++){
        free(B2E[i]);
    }
    
    // calculate normal vectors
    // allocate normal vector matrices
    double **IN = malloc(niedge * sizeof(double *));
    for (int i=0; i<niedge; i++){
        IN[i] = (double *) malloc(3 * sizeof(double));
    }
    double **BN = malloc(nbedge * sizeof(double *));
    for (int i=0; i<nbedge; i++){
        BN[i] = (double *) malloc(3 * sizeof(double));
    }
    // allocate perimeter matrix
    double *Per = calloc(Nelem, sizeof(double));
    
    // loop over interior edges
    double x1,y1,x2,y2;
    int n1 = 0, n2 = 0;
    double edgelen;
    for (int edge=0; edge<niedge; edge++){
        // get nodes of edge
        if (I2E[edge][1]==1){
            n1 = E[I2E[edge][0]-1][1];
            n2 = E[I2E[edge][0]-1][2];
        } else if (I2E[edge][1]==2){
            n1 = E[I2E[edge][0]-1][2];
            n2 = E[I2E[edge][0]-1][0];
        } else if (I2E[edge][1]==3){
            n1 = E[I2E[edge][0]-1][0];
            n2 = E[I2E[edge][0]-1][1];
        }
        // get coordinates of point
        x1 = V[n1-1][0]; y1 = V[n1-1][1];
        x2 = V[n2-1][0]; y2 = V[n2-1][1];
        // calculate normal
        edgelen = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
        IN[edge][0] = (y2-y1)/edgelen;
        IN[edge][1] = (x1-x2)/edgelen;
        // store edge length
        IN[edge][2] = edgelen;
        // add edge to perimeter tally
        Per[I2E[edge][0]-1] += edgelen;
        Per[I2E[edge][2]-1] += edgelen;
    }
    
    // loop over boundary edges
    for (int edge=0; edge<nbedge; edge++){
        // get nodes of edge
        if (B2E[edge][1]==1){
            n1 = E[B2E[edge][0]-1][1];
            n2 = E[B2E[edge][0]-1][2];
        } else if (B2E[edge][1]==2){
            n1 = E[B2E[edge][0]-1][2];
            n2 = E[B2E[edge][0]-1][0];
        } else if (B2E[edge][1]==3){
            n1 = E[B2E[edge][0]-1][0];
            n2 = E[B2E[edge][0]-1][1];
        }
        // get coordinates of point
        x1 = V[n1-1][0]; y1 = V[n1-1][1];
        x2 = V[n2-1][0]; y2 = V[n2-1][1];
        // calculate normal
        edgelen = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
        BN[edge][0] = (y2-y1)/edgelen;
        BN[edge][1] = (x1-x2)/edgelen;
        // ensure normal points out of domain
        if (B2E[edge][2]==1){ // left boundary
            if (BN[edge][0]>0){
                BN[edge][0] = -BN[edge][0];
                BN[edge][1] = -BN[edge][1];
            }
        } else if (B2E[edge][2]==2){// right boundary
            if (BN[edge][0]<0){
                BN[edge][0] = -BN[edge][0];
                BN[edge][1] = -BN[edge][1];
            }
        } else if (B2E[edge][2]==3){// top boundary
            if (BN[edge][1]<0){
                BN[edge][0] = -BN[edge][0];
                BN[edge][1] = -BN[edge][1];
            }
        } else if (B2E[edge][2]==4){// bottom boundary
            if (BN[edge][1]>0){
                BN[edge][0] = -BN[edge][0];
                BN[edge][1] = -BN[edge][1];
            }
        }
        // store edge length
        BN[edge][2] = edgelen;
        // add edge to perimeter tally
        Per[B2E[edge][0]-1] += edgelen;
        
    }
    
    // compute area, Jacobian of each cell
    double x3,y3;
    // allocate memory
    double *Area = malloc(Nelem * sizeof(double));
    double ***J = malloc(Nelem * sizeof(double **));
    double ***Jinv = malloc(Nelem * sizeof(double **));
    double *detJ = malloc(Nelem * sizeof(double));
    for (int k=0; k<Nelem; k++){
        J[k] = malloc(2 * sizeof(double *));
        Jinv[k] = malloc(2 * sizeof(double *));
        for (int i=0; i<2; i++){
            J[k][i] = malloc(2 * sizeof(double));
            Jinv[k][i] = malloc(2 * sizeof(double));
        }
    }
    
    // loop through each element
    for (int k=0; k<Nelem; k++){
        // area
        x1 = V[E[k][0]-1][0]; y1 = V[E[k][0]-1][1];
        x2 = V[E[k][1]-1][0]; y2 = V[E[k][1]-1][1];
        x3 = V[E[k][2]-1][0]; y3 = V[E[k][2]-1][1];
        Area[k] = 0.5 * fabs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
        
        // Jacobian
        x1 = V[E[k][0]-1][0]; y1 = V[E[k][0]-1][1];
        x2 = V[E[k][1]-1][0]; y2 = V[E[k][1]-1][1];
        x3 = V[E[k][2]-1][0]; y3 = V[E[k][2]-1][1];
        J[k][0][0] = x2-x1; J[k][0][1] = x3-x1;
        J[k][1][0] = y2-y1; J[k][1][1] = y3-y1;
        // determinant
        detJ[k] = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
        // inverse jacobian matrix
        Jinv[k][0][0] = (y3-y1)/detJ[k]; Jinv[k][0][1] = (x1-x3)/detJ[k];
        Jinv[k][1][0] = (y1-y2)/detJ[k]; Jinv[k][1][1] = (x2-x1)/detJ[k];
    }
    
    
    
    // free connectivity matrix
    for (int i=0; i<Nelem; i++){
        free(H[i]);
    }
    free(H);
    
    // free Boundaries array
    for (int btype = 0; btype < Nbgroups; btype++){
        // free each row of B[btype]
        for (int i = 0; i < Nbedges[btype]; i++){
            free(B[btype][i]);
        }
        // free the boundary array pointer
        free(B[btype]);
    }
    // free the top level Boundaries pointer
    free(B);
    
    
    // assemble output struct
    struct MeshData mesh;
    mesh.I2E = I2E;
    mesh.IN = IN;
    mesh.B2E = B2E;
    mesh.BN = BN;
    mesh.E = E;
    mesh.V = V;
    mesh.Per = Per;
    mesh.Area = Area;
    mesh.J = J;
    mesh.Jinv = Jinv;
    mesh.detJ = detJ;
    mesh.Nelem = Nelem;
    mesh.Nnodes = Nnodes;
    mesh.niedge = niedge;
    mesh.nbedge = nbedge;
    mesh.Nbgroups = Nbgroups;
    
    return mesh;
}





void freeMesh(struct MeshData mesh){
    // free mesh arrays
    // E
    for (int i=0; i<mesh.Nelem; i++){
        free(mesh.E[i]);
    }
    free(mesh.E);
    // V
    for (int i=0; i<mesh.Nnodes; i++){
        free(mesh.V[i]);
    }
    free(mesh.V);
    // I2E, IN
    for (int i=0; i<mesh.niedge; i++){
        free(mesh.I2E[i]);
        free(mesh.IN[i]);
    }
    free(mesh.I2E); free(mesh.IN);
    // B2E, BN
    for (int i=0; i<mesh.nbedge; i++){
        free(mesh.B2E[i]);
        free(mesh.BN[i]);
    }
    free(mesh.B2E);
    free(mesh.BN);
    // Perimeter, Area
    free(mesh.Per);
    free(mesh.Area);
    // Jacobians
    for (int k=0; k<mesh.Nelem; k++){
        for (int i=0; i<2; i++){
            free(mesh.J[k][i]);
            free(mesh.Jinv[k][i]);
        }
        free(mesh.J[k]);
        free(mesh.Jinv[k]);
    }
    free(mesh.J);
    free(mesh.Jinv);
    free(mesh.detJ);
}


