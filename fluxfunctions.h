//
//  fluxfunctions.h
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 3/23/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#ifndef fluxfunctions_h
#define fluxfunctions_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "precompute.h"

struct Parameters;

// function to calculate numerical flux on an inviscid wall boundary
double inviscidWall(double *fhat, double *ui, double *nvector, struct Parameters bc);

// function to calculate numerical flux on subsonic inflow boundary
double inflow(double *fhat, double *ui, double *nvector, struct Parameters bc);

// function to calculate numerical flux on subsonic outflow boundary
double outflow(double *fhat, double *ui, double *nvector, struct Parameters bc);

// function to calculate the numerical Roe flux (with an entropy fix)
double roeflux(double *fhat, const double *uL, const double *uR, const double *nvector, double gamma);

#endif /* fluxfunctions_h */
