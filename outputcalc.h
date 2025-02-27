//
//  outputcalc.h
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 4/3/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#ifndef outputcalc_h
#define outputcalc_h

#include <stdio.h>
#include <math.h>
#include "precompute.h"
#include "processmesh.h"
#include "calcRes.h"

// function to calculate and write outputs along bottom wall
void getWallCoeffs(double **U, struct MeshData mesh, struct ResData resdata, struct JData jcb, struct Parameters bc, int bumplevel);

// function to calculate and write entropy output for domain
void getEntropyError(double **U, struct MeshData mesh, struct ResData resdata, struct JData jcb, struct Parameters bc, int bumplevel);

#endif /* outputcalc_h */
