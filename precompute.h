//
//  precompute.h
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 3/19/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#ifndef precompute_h
#define precompute_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calcRes.h"
#include "processmesh.h"

struct MeshData;

struct Parameters{
    double gamma;
    double Minf;
    double alpha;
    double R;
    double Pinf;
    double Ttot;
    double Ptot;
    double h;
    double tolerance;
    double freestate[4];
};

struct ResData{
    int p;          // order of approximation
    int np;         // number of basis functions = (p+1)(p+2)/2
    int Q;          // order of geometry = p + 1
    int nQ;         // number of geometry-order basis functions = (Q+1)(Q+2)/2
    int n1;         // number of quad pts in 1D
    int n2;         // number of quad pts in 2D
    double *xq1;    // quadrature points, 1D
    double **xq2;   // quadrature points, 2D
    double *wq1;    // quadrature weights, 1D
    double *wq2;    // quadrature weights, 2D
};

struct JData{
    int Ncurvelem;          // number of curved elements
    int n1;                 // number of edge quadrature points per edge
    int n2;                 // number of interior quadrature points per triangle
    int *elemNum;          // kth entry gives the element number corresponding to Jint[k]
    double ***Jint;         // interior quadrature point Jacobian for curved elements: (k+q)->(2)x(2) Jacobian
    double ***Jintinv;      // interior quadrature point Jacobian inverse for curved elements
    double *detJint;        // interior quadrature point determinant(J) for curved elements (k+q)
    int *bedgeNum;          // kth entry gives the edge index in B2E corresponding to Jedge[k]
    double ***Jedge;        // interior quadrature point Jacobian for curved boundary edges: (k+q)->(2)x(2) Jacobian
    double ***Jedgeinv;     // interior quadrature point Jacobian inverse for curved boundary edges
    double *detJedge;       // interior quadrature point determinant(J) for curved boundary edges (k+q)
};

// function to get the residual evaluation quadrature data for a specified order p
struct ResData getResData(int p);
// function to free arrays allocated to the resdata struct
void freeResData(struct ResData resdata);

// function to compute jacobian values at all quad points across all elements in the mesh
struct JData getJData(struct MeshData mesh, struct ResData resdata);
// function to free arrays allocated to the JData struct
void freeJData(struct JData jcb);

// function to return the y coordinate of the exact bump geometry
double exactBumpY(double x);

// function to find the maximum value in a double array (2D, ie the residual)
double getLinfNorm(double **Res, int nrows, int ncols);

// function to compute the mass matrix based on p
void getMassMatrix(double ***M, struct MeshData mesh, struct ResData resdata, struct JData jcb);

// function to invert the mass matrix
void getMMinverse(double ***M, double ***Minv, struct MeshData mesh, struct ResData resdata);
// supporting functions for the inverse mass matrix:
double Determinant(double **a,int n);
void CoFactor(double **a,int n,double **b);
void Transpose(double **a,int n);

// function to initialize the state vector U
//void initializeU(double **U, )

#endif /* precompute_h */
