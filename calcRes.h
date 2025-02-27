//
//  calcRes.h
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 3/14/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#ifndef calcRes_h
#define calcRes_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "processmesh.h"
#include "precompute.h"
#include "fluxfunctions.h"

struct MeshData;
struct JData;
struct Parameters;
struct ResData;

// function to calculate residual over entire mesh
void calcRes(double **Res, double *wavespeed, double **U, const struct MeshData mesh, const struct ResData resdata, struct JData jcb, const struct Parameters bc);

// function to calculate the basis functions phi for a given p and quadrature point coordinates (xi,eta)
void calcPhi(double *phi, int p, double xi, double eta);

// function to calculate the gradient of basis functions phi for a given p and quadrature point coordinates (xi,eta)
void calcGradPhi(double **gradPhi, int p, double xi, double eta);

/*
// function to calculate interior edge states U for 1D edge residual contribution
void calciEdgeState(double *uL, double *uR, double **U, int edge, int q, const struct MeshData mesh, const struct ResData resdata);

// function to calculate boundary edge states U for 1D edge residual contribution
void calcbEdgeState(double *uL, double **U, int edge, int q, const struct MeshData mesh, const struct ResData resdata);
 */
#endif /* calcRes_h */
