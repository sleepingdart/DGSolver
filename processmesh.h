//
//  processmesh.h
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 3/12/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#ifndef processmesh_h
#define processmesh_h

#include <stdio.h>
#include <math.h>
#include "precompute.h"

struct ResData;

struct MeshData{
    int **I2E;      // interior edge matrix
    double **IN;    // interior edge normals (3rd column is edge length)
    int **B2E;      // boundary edge matrix
    double **BN;    // boundary edge normals (3rd column is edge length)
    int **E;        // element nodes
    double **V;     // node coordinates
    double *Per;    // perimeter of each element
    double *Area;   // area of each element
    double ***J;    // Jacobian matrix per linear element
    double ***Jinv; // inverse Jacobian matrix per linear element
    double *detJ;   // determinant of Jacobian matrix per linear element
    int Nelem;      // number of elements (total)
    int Nnodes;     // number of nodes
    int niedge;     // number of interior edges
    int nbedge;     // number of boundary edges
    int Nbgroups;   // number of boundary groups
};

// function to read gri file and process relevant data
struct MeshData processmesh(char *grifile);
// function to free the arrays created in processmesh
void freeMesh(struct MeshData);


#endif /* processmesh_h */
