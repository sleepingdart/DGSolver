//
//  fluxfunctions.c
//  CFD2_Project3_Code
//
//  Created by Borren Moe on 3/23/19.
//  Copyright © 2019 Borren Moe. All rights reserved.
//

#include "fluxfunctions.h"


double inviscidWall(double *fhat, double *ui, double *nvector, struct Parameters bc){
    // calculate flux on an inviscid wall boundary
    
    // pull variable from struct
    double gamma = bc.gamma;
    
    double vi[2] = {(ui[1]/ui[0]),(ui[2]/ui[0])};
    double speedi = sqrt(vi[0]*vi[0] + vi[1]*vi[1]);
    double Pi = fabs((gamma-1)*(ui[3]-0.5*ui[0]*pow(speedi,2)));
    double ci = sqrt(gamma*Pi/ui[0]);
    
    // calculate flux
    fhat[0]=0; fhat[1]=Pi*nvector[0]; fhat[2]=Pi*nvector[1]; fhat[3]=0;
    
    // calculate wavespeed
    double un = vi[0]*nvector[0] + vi[1]*nvector[1];
    double wavespeed = un + ci;
    
    return wavespeed;
}



double inflow(double *fhat, double *ui, double *nvector, struct Parameters bc){
    // calculate flux on a subsonic inflow boundary
    
    // pull variables from struct
    double gamma = bc.gamma;
    double Ttot = bc.Ttot;
    double Ptot = bc.Ptot;
    double alpha = bc.alpha;
    double R = bc.R;
    
    double vi[2] = {(ui[1]/ui[0]),(ui[2]/ui[0])};
    double speedi = sqrt(pow(vi[0],2) + pow(vi[1],2));
    double Pi = fabs((gamma-1)*(ui[3]-0.5*ui[0]*pow(speedi,2)));
    double ci = sqrt(gamma*Pi/ui[0]);
    
    double Jplus = vi[0]*nvector[0]+vi[1]*nvector[1] + 2*ci/(gamma-1);
    double nin[2] = {cos(alpha),sin(alpha)};
    double dn = nin[0]*nvector[0]+nin[1]*nvector[1];
    
    // solve quadratic for boundary Mach number
    double A = gamma*R*Ttot*dn*dn - (gamma-1)/2*Jplus*Jplus;
    double B = 4*gamma*R*Ttot*dn / (gamma-1);
    double C = 4*gamma*R*Ttot/pow(gamma-1,2) - Jplus*Jplus;
    
    double Mbsoln[2] = {(-B+sqrt(B*B-4*A*C))/(2*A), (-B-sqrt(B*B-4*A*C))/(2*A)};
    double Mb;
    
    if ((Mbsoln[0]>=0) && (Mbsoln[1]>=0)) {
        // take the minimum value
        Mb = (Mbsoln[0] > Mbsoln[1]) ? Mbsoln[1] : Mbsoln[0];
        // otherwise, use the positive solution
    } else if (Mbsoln[0]>0){
        Mb = Mbsoln[0];
    } else if (Mbsoln[1]>0){
        Mb = Mbsoln[1];
    } else{
        // both solns negative, take closest to zero to prevent errors (should never occur)
        Mb = (Mbsoln[0] < Mbsoln[1]) ? Mbsoln[1] : Mbsoln[0];
    }
    
    // calculate exterior state
    double Tb = Ttot / (1+0.5*(gamma-1)*Mb*Mb);
    double Pb = Ptot*pow((Tb/Ttot),(gamma/(gamma-1)));
    double rhob = Pb / (R*Tb);
    double cb = sqrt(gamma*Pb/rhob);
    double vb[2] = {Mb*cb*nin[0], Mb*cb*nin[1]};
    double rhoEb = Pb/(gamma-1) + 0.5*rhob*(vb[0]*vb[0]+vb[1]*vb[1]);
    
    //assemble flux
    double Hb = (rhoEb+Pb)/rhob;
    double Fb[4][2] = {
        {rhob*vb[0], rhob*vb[1]},
        {rhob*vb[0]*vb[0]+Pb, rhob*vb[1]*vb[0]},
        {rhob*vb[0]*vb[1], rhob*vb[1]*vb[1]+Pb},
        {rhob*vb[0]*Hb, rhob*vb[1]*Hb}
    };
    
    fhat[0] = Fb[0][0]*nvector[0] + Fb[0][1]*nvector[1];
    fhat[1] = Fb[1][0]*nvector[0] + Fb[1][1]*nvector[1];
    fhat[2] = Fb[2][0]*nvector[0] + Fb[2][1]*nvector[1];
    fhat[3] = Fb[3][0]*nvector[0] + Fb[3][1]*nvector[1];
    
    // calculate wavespeed
    double un = vi[0]*nvector[0] + vi[1]*nvector[1];
    double wavespeed = un + ci;
    
    return wavespeed;
}



double outflow(double *fhat, double *ui, double *nvector, struct Parameters bc){
    // calculate flux on a subsonic outflow boundary
    
    // pull variables from struct
    double gamma = bc.gamma;
    double Pb = bc.Pinf;
    
    double vi[2] = {(ui[1]/ui[0]),(ui[2]/ui[0])};
    double speedi = sqrt(pow(vi[0],2) + pow(vi[1],2));
    double Pi = fabs((gamma-1)*(ui[3]-0.5*ui[0]*pow(speedi,2)));
    double ci = sqrt(gamma*Pi/ui[0]);
    
    double Jplus = vi[0]*nvector[0]+vi[1]*nvector[1] + 2*ci/(gamma-1);
    double Si = Pi/pow(ui[0],gamma); // interior entropy
    
    double rhob = pow((Pb/Si),(1/gamma));
    
    double cb = sqrt(gamma*Pb/rhob);
    double ubn = Jplus - 2*cb/(gamma-1);
    
    double vb[2] = {
        vi[0] - (vi[0]*nvector[0]+vi[1]*nvector[1])*nvector[0] + ubn*nvector[0],
        vi[1] - (vi[0]*nvector[0]+vi[1]*nvector[1])*nvector[1] + ubn*nvector[1]
    };
    
    
    double rhoEb = Pb/(gamma-1) + 0.5*rhob*(vb[0]*vb[0]+vb[1]*vb[1]);
    
    // assemble flux
    double Hb = (rhoEb+Pb)/rhob;
    
    double Fb[4][2] = {
        {rhob*vb[0], rhob*vb[1]},
        {rhob*vb[0]*vb[0]+Pb, rhob*vb[1]*vb[0]},
        {rhob*vb[0]*vb[1], rhob*vb[1]*vb[1]+Pb},
        {rhob*vb[0]*Hb, rhob*vb[1]*Hb}
    };
    
    fhat[0] = Fb[0][0]*nvector[0] + Fb[0][1]*nvector[1];
    fhat[1] = Fb[1][0]*nvector[0] + Fb[1][1]*nvector[1];
    fhat[2] = Fb[2][0]*nvector[0] + Fb[2][1]*nvector[1];
    fhat[3] = Fb[3][0]*nvector[0] + Fb[3][1]*nvector[1];
    
    // calculate wavespeed
    double un = vi[0]*nvector[0] + vi[1]*nvector[1];
    double wavespeed = un + ci;
    
    return wavespeed;
}




double roeflux(double *fhat, const double *uL, const double *uR, const double *nvector, double gamma){
    // pass in pointer to fhat vector, returns wavespeed value
    
    // compute velocity vectors
    double VL[2] = {uL[1]/uL[0], uL[2]/uL[0]};
    double VR[2] = {uR[1]/uR[0], uR[2]/uR[0]};
    
    // compute pressures
    double PL = (gamma-1)*(uL[3]-0.5*uL[0]*(VL[0]*VL[0]+VL[1]*VL[1]));
    double PR = (gamma-1)*(uR[3]-0.5*uR[0]*(VR[0]*VR[0]+VR[1]*VR[1]));
    
    // compute enthalpies
    double HL = (uL[3]+PL)/uL[0];
    double HR = (uR[3]+PR)/uR[0];
    
    // compute flux vectors
    double FL[4][2] = {
        {uL[1], uL[2]},
        {uL[1]*uL[1]/uL[0]+PL, uL[1]*uL[2]/uL[0]},
        {uL[1]*uL[2]/uL[0], uL[2]*uL[2]/uL[0]+PL},
        {uL[1]*uL[3]/uL[0]+uL[1]*PL/uL[0], uL[2]*uL[3]/uL[0]+uL[2]*PL/uL[0]}
    };
    
    double FR[4][2] = {
        {uR[1], uR[2]},
        {uR[1]*uR[1]/uR[0]+PR, uR[1]*uR[2]/uR[0]},
        {uR[1]*uR[2]/uR[0], uR[2]*uR[2]/uR[0]+PR},
        {uR[1]*uR[3]/uR[0]+uR[1]*PR/uR[0], uR[2]*uR[3]/uR[0]+uR[2]*PR/uR[0]}
    };
    
    // compute Roe state variables
    double Vroe[2] = {
        (sqrt(uL[0])*VL[0] + sqrt(uR[0])*VR[0]) / (sqrt(uL[0])+sqrt(uR[0])),
        (sqrt(uL[0])*VL[1] + sqrt(uR[0])*VR[1]) / (sqrt(uL[0])+sqrt(uR[0]))
    };
    double Hroe = (sqrt(uL[0])*HL + sqrt(uR[0])*HR) / (sqrt(uL[0])+sqrt(uR[0]));
    
    double q = sqrt(Vroe[0]*Vroe[0]+Vroe[1]*Vroe[1]);
    double c = sqrt((gamma-1)*(Hroe-0.5*q*q));
    
    // apply entropy fix to eigenvalues
    double epsilon = 0.1*c;
    double lambd[4] = {
        (Vroe[0]*nvector[0]+Vroe[1]*nvector[1])+c,
        (Vroe[0]*nvector[0]+Vroe[1]*nvector[1])-c,
        (Vroe[0]*nvector[0]+Vroe[1]*nvector[1]),
        (Vroe[0]*nvector[0]+Vroe[1]*nvector[1])
    };
    
    for(int i=0; i<4; i++){
        if (fabs(lambd[i]) < epsilon){
            lambd[i] = (epsilon*epsilon+lambd[i]*lambd[i])/(2*epsilon);
        }
    };
    
    // assemble analytical flux vector
    double s1 = 0.5*(fabs(lambd[0])+fabs(lambd[1]));
    double s2 = 0.5*(fabs(lambd[0])-fabs(lambd[1]));
    
    double G1 = (gamma-1)*(q*q/2*(uR[0]-uL[0]) - (Vroe[0]*(uR[1]-uL[1])+Vroe[1]*(uR[2]-uL[2])) + (uR[3]-uL[3]));
    double G2 = -(Vroe[0]*nvector[0]+Vroe[1]*nvector[1])*(uR[0]-uL[0]) + ((uR[1]-uL[1])*nvector[0]+(uR[2]-uL[2])*nvector[1]);
    
    double C1 = G1/(c*c)*(s1-fabs(lambd[2])) + G2/c*s2;
    double C2 = G1/c*s2 + (s1-fabs(lambd[2]))*G2;
    
    double term2[4] = {
        fabs(lambd[2])*(uR[0]-uL[0]) + C1,
        fabs(lambd[2])*(uR[1]-uL[1]) + C1*Vroe[0] + C2*nvector[0],
        fabs(lambd[2])*(uR[2]-uL[2]) + C1*Vroe[1] + C2*nvector[1],
        fabs(lambd[2])*(uR[3]-uL[3]) + C1*Hroe + C2*(Vroe[0]*nvector[0]+Vroe[1]*nvector[1])
    };
    
    // generate output
    fhat[0] = 0.5*( (FL[0][0]*nvector[0]+FL[0][1]*nvector[1]) + (FR[0][0]*nvector[0]+FR[0][1]*nvector[1]) - term2[0] );
    fhat[1] = 0.5*( (FL[1][0]*nvector[0]+FL[1][1]*nvector[1]) + (FR[1][0]*nvector[0]+FR[1][1]*nvector[1]) - term2[1] );
    fhat[2] = 0.5*( (FL[2][0]*nvector[0]+FL[2][1]*nvector[1]) + (FR[2][0]*nvector[0]+FR[2][1]*nvector[1]) - term2[2] );
    fhat[3] = 0.5*( (FL[3][0]*nvector[0]+FL[3][1]*nvector[1]) + (FR[3][0]*nvector[0]+FR[3][1]*nvector[1]) - term2[3] );
    
    // return wavespeed
    return lambd[0];// * nvector[2];
}
