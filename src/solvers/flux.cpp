#include "build.h"
#include <cmath>
#include <array>

void SEMO_Solvers_Builder::calcFluxHAUS(
	double& PL, double& UL, double& VL, double& WL, double& TL, 
	vector<double>& YL, double& rhoL, double& CL, double& HtL, 
	double& PR, double& UR, double& VR, double& WR, double& TR, 
	vector<double>& YR, double& rhoR, double& CR, double& HtR,
	vector<double>& nvec,
	vector<double>& flux
	){
	

	
    // properties of Left
    double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
    // properties of Right
    double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
	
	double unhat = (UnL+UnR)/2.0;
	double rhohat = 0.5*(rhoL+rhoR);
    double chat= 0.5*(CL+CR);
	
	double ML = UnL/chat; 
	double MR = UnR/chat;
	double KLR = sqrt(0.5*(UL*UL+VL*VL+WL*WL+UR*UR+VR*VR+WR*WR));
	double MLP  = M_func(ML,1.0,0.0); 
	double MRM = M_func(MR,-1.0,0.0);
	double preP = pre_func(ML,1.0,0.0); 
	double preM = pre_func(MR,-1.0,0.0);
	
    double Mbar = ( rhoL*abs(ML)+rhoR*abs(MR) ) / ( rhoL + rhoR );

	// SLAU
	double Mcy = min(1.0,KLR/chat);
	double phi_c = pow((1.0-Mcy),2.0);
	double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );

    double D_L = ML+(1.0-g_c)*abs(ML);
    double D_R = MR-(1.0-g_c)*abs(MR);
    double D_rho = Mbar*g_c;

	double preLs = abs(PL) + 0.1 * rhoL*CL*CL;
	double preRs = abs(PR) + 0.1 * rhoR*CR*CR;
	double w = 1.0 - pow( min(preLs/preRs,preRs/preLs),2.0);

	double ps = preP*PL+preM*PR;
	double pll;
	if( 3.0/4.0 <= min(PL/PR,PR/PL) && 1.0 > min(PL/PR,PR/PL) ){
		pll=4.0*min(PL/PR,PR/PL)-3.0;
	}
	else{
		pll=0.0;
	}
	double fL;
	if( abs(ML) <= 1.0 ) {
		fL = (PL/ps-1.0)*pll*abs(MLP)*min(1.0,pow( (abs(UnL)/chat),0.25 ));
	}
	else{
		fL = 0.0;
	}
		
	double fR;
	if( abs(MR) <= 1.0 ) {
		fR = (PL/ps-1.0)*pll*abs(MRM)*min(1.0,pow( (abs(UnR)/chat),0.25 ));
	}
	else{
		fR = 0.0;
	}

	double MPP = MLP+MRM;
	double MLP_AUSM = 0.5*(MPP+abs(MPP));
	double MRM_AUSM = 0.5*(MPP-abs(MPP));

	double MLP_SLAU = 0.5*(D_L+D_rho);
	double MRM_SLAU = 0.5*(D_R-D_rho);

	double fa1 = w;

	double MLPL = fa1*MLP_AUSM + (1.0-fa1)*MLP_SLAU;
	double MRMR = fa1*MRM_AUSM + (1.0-fa1)*MRM_SLAU;

	double mdot = rhoL*chat*MLPL + rhoR*chat*MRMR
	   - 0.5
		 *(1.0-0.5*(1.0-cos(3.141592*min(1.0,max(abs(ML),abs(MR))))))
		 *(1.0-0.5*(1.0+cos(3.141592*min(abs(PL/PR),abs(PR/PL)))))
		 *(PR-PL)
		 /chat;

	MLP = 0.5*(ML+abs(ML));
	MRM = 0.5*(MR-abs(MR));
	double f1L;
	double f1R;
    if( mdot >= 0.0 ) {
      f1L = mdot - (rhoL*chat*MRM)*( w*(1.0+fR)-fR+fL );
      f1R = (rhoR*chat*MRM)*( w*(1.0+fR) );
	}
    else{
      f1L = (rhoL*chat*MLP)*( w*(1.0+fL) );
      f1R = mdot - (rhoR*chat*MLP)*( w*(1.0+fL)-fL+fR );
    }
	
	// Pressure Term
	// if( cosface <= 1.d0 .and. cosface >= -0.5d0 ) then
	  // gam3 = 1.d0 - 1.d0*(2.d0/3.d0*cosface+1.d0/3.d0)
	// else; gam3 = 1.d0; endif
	
	double fa2 = w;
	// gam2 = 0.5d0*(1.d0-tanh(15.d0*(min(pwL,pwR)+0.0d0))) !1.d0-dmin1(dabs(pL/pR),dabs(pR/pL))**3.d0
	// gam = dmax1( 0.d0, gam2 )
	
	// double gam2 = 0.5*(1.0-tanh(15.5*cpi2));
	// double gam_HR = 0.6;
	// double gam = dmax1( 0.0, gam_HR, gam2 );
	double gam = 0.6;

	double PLR = 0.5*(PL+PR) 
				- fa2*(KLR/chat)*0.5*preP*preM*0.5*(PL+PR)/chat*(UnR-UnL)
				+ max(0.2,gam)*(KLR/chat)*0.5*(PL+PR)*(preP+preM-1.0)
				- 0.5*(preP-preM)*(PR-PL);

	flux.clear();
	flux.push_back(f1L+f1R);
	flux.push_back(f1L*UL + f1R*UR + PLR*nvec[0]);
	flux.push_back(f1L*VL + f1R*VR + PLR*nvec[1]);
	flux.push_back(f1L*WL + f1R*WR + PLR*nvec[2]);
	flux.push_back(f1L*HtL+ f1R*HtR);
	
	for(int i=0; i<YL.size()-1; ++i){
		flux.push_back(f1L*YL[i]+ f1R*YR[i]);
	}
	

}


double SEMO_Solvers_Builder::M_func(double M, double op, double alp){

	double mu;

	if( abs(M) > 1.0 ) {
		mu = 0.5*(M+op*abs(M));
	}
	else{
		mu = op*0.25*pow((M + op),2.0) + op*alp*pow((M*M-1.0),2.0);
	}
	
	return mu;
}
double SEMO_Solvers_Builder::pre_func(double M, double op, double alp){

	double mu;

	if( abs(M) > 1.0 ) {
		mu = 0.5*(1.0 + op * ( M>0.0 ? 1.0 : -1.0 ) );
	}
	else{
		mu = 0.25*pow((M + op),2.0)*(2.0-op*M) + op*alp*M*pow((M*M-1.0),2.0);
	}
	
	return mu;
}