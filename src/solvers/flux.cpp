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
	
	// if( 
	// abs(ML) < std::numeric_limits<double>::epsilon() &&
	// abs(MR) < std::numeric_limits<double>::epsilon()){
		// f1L = 0.0;
		// f1R = 0.0;
	// }
	
	
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

	// if(MPP>=0.0){
		// f1L = MPP;
		// f1R = 0.0;
	// }
	// else{
		// f1L = 0.0;
		// f1R = MPP;
	// }
	// double PLR = preP*PL+preM*PR;
	
	
	
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



void SEMO_Solvers_Builder::calcFluxSLAU2(
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

	double MLPL = 0.5*(D_L+D_rho);
	double MRMR = 0.5*(D_R-D_rho);

	double mdot = rhoL*chat*MLPL + rhoR*chat*MRMR - 0.5*phi_c/chat*(PR-PL);

	double f1L;
	double f1R;
    if( mdot >= 0.0 ) {
		f1L = mdot; f1R = 0.0;
	}
    else{
		f1L = 0.0; f1R = mdot;
    }
	
	double PLR = 0.5*(PL+PR) - 
				 0.5*Mcy*preP*preM*0.5*(PL+PR)/chat*(UnR-UnL) + 
				 Mcy*0.5*(PL+PR)*(preP+preM-1.0) - 
				 0.5*(preP-preM)*(PR-PL);

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









void SEMO_Solvers_Builder::calcFluxRoe(
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
    double Uhat= 0.5*(UL+UR);
    double Vhat= 0.5*(VL+VR);
    double What= 0.5*(WL+WR);
    double Hthat= 0.5*(HtL+HtR);
	vector<double> Yihat;
	for(int i=0; i<YL.size()-1; ++i){
		Yihat.push_back( 0.5*(YL[i]+YR[i]) );
	}
	
	// left advection flux
	vector<double> fluxAdvL;
	fluxAdvL.push_back(rhoL*UnL);
	fluxAdvL.push_back(rhoL*UnL*UL + PL*nvec[0]);
	fluxAdvL.push_back(rhoL*UnL*VL + PL*nvec[1]);
	fluxAdvL.push_back(rhoL*UnL*WL + PL*nvec[2]);
	fluxAdvL.push_back(rhoL*UnL*HtL);
	for(int i=0; i<YL.size()-1; ++i){
		fluxAdvL.push_back(rhoL*UnL*YL[i]);
	}
	
	// right advection flux
	vector<double> fluxAdvR;
	fluxAdvR.push_back(rhoR*UnR);
	fluxAdvR.push_back(rhoR*UnR*UR + PR*nvec[0]);
	fluxAdvR.push_back(rhoR*UnR*VR + PR*nvec[1]);
	fluxAdvR.push_back(rhoR*UnR*WR + PR*nvec[2]);
	fluxAdvR.push_back(rhoR*UnR*HtR);
	for(int i=0; i<YR.size()-1; ++i){
		fluxAdvR.push_back(rhoR*UnR*YR[i]);
	}
	
	
	// dissipation flux
	vector<double> fluxDissDW;
	double tmpUn = abs(unhat);
	fluxDissDW.push_back( tmpUn*(rhoR-rhoL) );
	fluxDissDW.push_back( tmpUn*(rhoR*UR-rhoL*UL) );
	fluxDissDW.push_back( tmpUn*(rhoR*VR-rhoL*VL) );
	fluxDissDW.push_back( tmpUn*(rhoR*WR-rhoL*WL) );
	fluxDissDW.push_back( tmpUn*(rhoR*(HtR-PR/rhoR)-rhoL*(HtL-PL/rhoL)) );
	for(int i=0; i<YL.size()-1; ++i){
		fluxDissDW.push_back( tmpUn*(rhoR*YR[i]-rhoL*YL[i]) );
	}
	
	vector<double> fluxDissDP;
	double tmpDP = unhat/chat*(PR-PL);
	fluxDissDP.push_back(0.0);
	fluxDissDP.push_back(tmpDP*nvec[0]);
	fluxDissDP.push_back(tmpDP*nvec[1]);
	fluxDissDP.push_back(tmpDP*nvec[2]);
	fluxDissDP.push_back(tmpDP*unhat);
	for(int i=0; i<YL.size()-1; ++i){
		fluxDissDP.push_back(0.0);
	}
	
	vector<double> fluxDissDU;
	double tmpDU = (chat-abs(unhat))/chat/chat/rhohat*(PR-PL) + unhat/chat*(UnR-UnL);
	fluxDissDU.push_back(tmpDU*rhohat);
	fluxDissDU.push_back(tmpDU*rhohat*Uhat);
	fluxDissDU.push_back(tmpDU*rhohat*Vhat);
	fluxDissDU.push_back(tmpDU*rhohat*What);
	fluxDissDU.push_back(tmpDU*rhohat*Hthat);
	for(int i=0; i<YL.size()-1; ++i){
		fluxDissDU.push_back(tmpDU*rhohat*Yihat[i]);
	}
	
	// flux
	flux.clear();
	for(int i=0; i<5+(YL.size()-1); ++i){
		// flux.push_back( 
			// 0.5*(fluxAdvL[i]+fluxAdvR[i]) 
			// - 0.5*(fluxDissDW[i]+fluxDissDP[i]+fluxDissDU[i]) );
			
		flux.push_back( 
			0.5*(fluxAdvL[i]+fluxAdvR[i]) );
	}
	

}








void SEMO_Solvers_Builder::calcFluxAPRoe(
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
	 
    double RT = sqrt(rhoR/rhoL); 
	double uhat = (UL+RT*UR)/(1.0+RT); 
	double vhat = (VL+RT*VR)/(1.0+RT);
	double what = (WL+RT*WR)/(1.0+RT); 
	double Hthat = (HtL+RT*HtR)/(1.0+RT); 
	double rhohat = RT*rhoL;
	vector<double> Yihat;
	for(int i=0; i<YL.size()-1; ++i){
		Yihat.push_back( (YL[i]+RT*YR[i])/(1.0+RT) );
	}
	
    double chat= (CL+RT*CR)/(1.0+RT); 
	double unhat = (UnL+RT*UnR)/(1.0+RT); 
	
    //======= carefully, sensitive ===========
    double preLs = PL + 0.1 * min(rhoL*CL*CL,rhoR*CR*CR);
    double preRs = PR + 0.1 * min(rhoL*CL*CL,rhoR*CR*CR);
	//========================================
    double cpi = min(preLs/preRs,preRs/preLs);
    double cpi3 = cpi;
	
	
	// left advection flux
	vector<double> fluxAdvL;
	fluxAdvL.push_back(rhoL*UnL);
	fluxAdvL.push_back(rhoL*UnL*UL + PL*nvec[0]);
	fluxAdvL.push_back(rhoL*UnL*VL + PL*nvec[1]);
	fluxAdvL.push_back(rhoL*UnL*WL + PL*nvec[2]);
	fluxAdvL.push_back(rhoL*UnL*HtL);
	for(int i=0; i<YL.size()-1; ++i){
		fluxAdvL.push_back(rhoL*UnL*YL[i]);
	}
	
	// right advection flux
	vector<double> fluxAdvR;
	fluxAdvR.push_back(rhoR*UnR);
	fluxAdvR.push_back(rhoR*UnR*UR + PR*nvec[0]);
	fluxAdvR.push_back(rhoR*UnR*VR + PR*nvec[1]);
	fluxAdvR.push_back(rhoR*UnR*WR + PR*nvec[2]);
	fluxAdvR.push_back(rhoR*UnR*HtR);
	for(int i=0; i<YR.size()-1; ++i){
		fluxAdvR.push_back(rhoR*UnR*YR[i]);
	}
	
	
	// dissipation flux
	vector<double> DW;
	DW.push_back( (rhoR-rhoL) );
	DW.push_back( (rhoR*UR-rhoL*UL) );
	DW.push_back( (rhoR*VR-rhoL*VL) );
	DW.push_back( (rhoR*WR-rhoL*WL) );
	DW.push_back( (rhoR*(HtR-PR/rhoR)-rhoL*(HtL-PL/rhoL)) );
	for(int i=0; i<YL.size()-1; ++i){
		DW.push_back( (rhoR*YR[i]-rhoL*YL[i]) );
	}
	
	double xi = abs(unhat);
	
	double Mcy = sqrt(0.5*(UL*UL+VL*VL+WL*WL+UR*UR+VR*VR+WR*WR))/chat;
	
	double theta = min(max(Mcy*Mcy,0.01*0.01),1.0);
	double Um = unhat*min(max(Mcy,0.01),1.0);
	double Cm = chat*min(max(Mcy,0.01),1.0);
	
	double thetad = min(Mcy*Mcy,1.0);
	double Ud = 0.5*(1.0+thetad);
	double Cd = 0.5*sqrt(4.0*chat*chat*thetad+pow((1.0-thetad),2.0)*unhat*unhat);
	
	vector<double> DU(5+YL.size()-1,0.0);
	DU[0] = (Cm - 0.5*(1.0-theta)*unhat*Um/Cm - theta*abs(unhat))*(PR-PL)/rhohat/theta/chat/chat
			+ Ud/Cm*(UnR-UnL);
	DU[0] = DU[0]*rhohat;
	DU[1] = DU[0]*uhat; 
	DU[2] = DU[0]*vhat;	
	DU[3] = DU[0]*what;
	DU[4] = DU[0]*Hthat;
	for(int i=0; i<YL.size()-1; ++i){
		DU[5+i] = DU[0]*Yihat[i];
	}
	
    //> X. Li et al., 2017
	vector<double> DP(5+YL.size()-1,0.0);	
	DP[0] = (Cd - abs(unhat) + 0.5*(1.0-thetad)*unhat*Ud/Cm)*rhohat*(UnR-UnL) + Ud*(PR-PL)/Cm;
	DP[1] = DP[0]*nvec[0]; 
	DP[2] = DP[0]*nvec[1];
	DP[4] = DP[0]*nvec[2]; 
	DP[0] = 0.0;
	
	
	
	// flux
	flux.clear();
	for(int i=0; i<5+(YL.size()-1); ++i){
		flux.push_back( 
			0.5*(fluxAdvL[i]+fluxAdvR[i]) - 0.5*(xi*DW[i] + DP[i] + DU[i]) );
	}
	

}








/*==============================================================================
!> @author yyl
!> @date 2018-03-20
!> @brief comp. RoeM Riemann solver for multiphase
!> Ref : Preconditioning Applied to Variable and Constant Density Flows, 1995
!> Ref : Role of Momentum Interpolation Mechanism of the Roe
!> Scheme in Shock Instability, 2017
*/
void SEMO_Solvers_Builder::calcFluxRoeAMTP(
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
	 
    double RT = sqrt(rhoR/rhoL); 
	double uhat = (UL+RT*UR)/(1.0+RT); 
	double vhat = (VL+RT*VR)/(1.0+RT);
	double what = (WL+RT*WR)/(1.0+RT); 
	double Hthat = (HtL+RT*HtR)/(1.0+RT); 
	double rhohat = RT*rhoL;
	vector<double> Yihat;
	for(int i=0; i<YL.size()-1; ++i){
		Yihat.push_back( (YL[i]+RT*YR[i])/(1.0+RT) );
	}
	
    double chat= (CL+RT*CR)/(1.0+RT); 
	double unhat = (UnL+RT*UnR)/(1.0+RT); 
	
    //======= carefully, sensitive ===========
    double preLs = PL + 0.1 * min(rhoL*CL*CL,rhoR*CR*CR);
    double preRs = PR + 0.1 * min(rhoL*CL*CL,rhoR*CR*CR);
	//========================================
    double cpi = min(preLs/preRs,preRs/preLs);
    double cpi3 = cpi;
	
	
	// left advection flux
	vector<double> fluxAdvL;
	fluxAdvL.push_back(rhoL*UnL);
	fluxAdvL.push_back(rhoL*UnL*UL + PL*nvec[0]);
	fluxAdvL.push_back(rhoL*UnL*VL + PL*nvec[1]);
	fluxAdvL.push_back(rhoL*UnL*WL + PL*nvec[2]);
	fluxAdvL.push_back(rhoL*UnL*HtL);
	for(int i=0; i<YL.size()-1; ++i){
		fluxAdvL.push_back(rhoL*UnL*YL[i]);
	}
	
	// right advection flux
	vector<double> fluxAdvR;
	fluxAdvR.push_back(rhoR*UnR);
	fluxAdvR.push_back(rhoR*UnR*UR + PR*nvec[0]);
	fluxAdvR.push_back(rhoR*UnR*VR + PR*nvec[1]);
	fluxAdvR.push_back(rhoR*UnR*WR + PR*nvec[2]);
	fluxAdvR.push_back(rhoR*UnR*HtR);
	for(int i=0; i<YR.size()-1; ++i){
		fluxAdvR.push_back(rhoR*UnR*YR[i]);
	}
	
	
	// dissipation flux
	vector<double> delW;
	delW.push_back( (rhoR-rhoL) );
	delW.push_back( (rhoR*UR-rhoL*UL) );
	delW.push_back( (rhoR*VR-rhoL*VL) );
	delW.push_back( (rhoR*WR-rhoL*WL) );
	delW.push_back( (rhoR*(HtR-PR/rhoR)-rhoL*(HtL-PL/rhoL)) );
	for(int i=0; i<YL.size()-1; ++i){
		delW.push_back( (rhoR*YR[i]-rhoL*YL[i]) );
	}
	
	double Morg = 0.5*(sqrt(UL*UL+VL*VL+WL*WL)/CL+sqrt(UR*UR+VR*VR+WR*WR)/CR);
	double Mbar = abs(unhat)/chat;
	double fef = 0.05*chat;
	
	double ksqrt = sqrt(pow((UR-UL),2.0)+pow((VR-VL),2.0)+pow((WR-WL),2.0));
	double n1x, n1y, n1z;
	if( ksqrt < 1.e-5 ) {
	  n1x = nvec[0]; n1y = nvec[1]; n1z = nvec[2];
	}
	else{
	  n1x = (UR-UL)/ksqrt; n1y = (VR-VL)/ksqrt; n1z = (WR-WL)/ksqrt;
	}
	double n2x = - n1y*(n1x*nvec[1] - n1y*nvec[0]) - n1z*(n1x*nvec[2] - n1z*nvec[0]);
	double n2y = n1x*(n1x*nvec[1] - n1y*nvec[0]) - n1z*(n1y*nvec[2] - n1z*nvec[1]);
	double n2z = n1x*(n1x*nvec[2] - n1z*nvec[0]) + n1y*(n1y*nvec[2] - n1z*nvec[1]);
	
	double a1=n1x*nvec[0]+n1y*nvec[1]+n1z*nvec[2];	
	double a2=n2x*nvec[0]+n2y*nvec[1]+n2z*nvec[2];
	double U1=n1x*uhat+n1y*vhat+n1z*what;	
	double U2=n2x*uhat+n2y*vhat+n2z*what;
	double frr = abs(a1*U1)+abs(a2*U2);
			 
	double xi = max( abs(0.5*(UnR+UnL))+0.5*(UnR-UnL) , 
	               (1.0-pow( f_func(Mbar),8.0 ))* pow( f_func(Morg),8.0 )*min(fef,frr) );
	double Und = abs(unhat) - pow( f_func(Mbar),8.0 )*(
			( (unhat+chat > 0.0) ? 1.0 : -1.0 )*max(0.0,UnR-UnL)-
			( (unhat-chat > 0.0) ? 1.0 : -1.0 )*max(0.0,UnR-UnL))/4.0;
											
	double theta = 1.0;
    double Umul = 0.5*(1.0+theta)*unhat;
    double Cmul = 0.5*sqrt(4.0*chat*chat*theta+pow((1.0-theta),2.0)*unhat*unhat);
	cpi3 = min( abs(PL/PR), abs(PR/PL) );
    double s1 = pow(f_func(cpi3),8.0);
	
	vector<double> dUp(5+YL.size()-1,0.0);
	dUp[0] = s1*(1.0-pow( f_func(Morg),8.0 ))*(max(0.0,Cmul-Und)+
	         (1.0-theta)*(Und-unhat*0.5/Cmul* 
			 ( (Umul > 0.0) ? 1.0 : -1.0 )*min(abs(Umul),Cmul)))*
			 (PR-PL)/(theta*chat*chat);
	dUp[1] = dUp[0]*uhat; 
	dUp[2] = dUp[0]*vhat;	
	dUp[3] = dUp[0]*what;
	dUp[4] = dUp[0]*Hthat;
	for(int i=0; i<YL.size()-1; ++i){
		dUp[5+i] = dUp[0]*Yihat[i];
	}
	
	vector<double> dUuLR(5+YL.size()-1,0.0);
    double xi_LR = 0.5*( (unhat > 0.0) ? 1.0 : -1.0 )*min(Und,chat)*(UnR-UnL)/chat;
    dUuLR[0] = xi_LR*(rhoL+rhoR);
	dUuLR[1] = xi_LR*(rhoL*UL+rhoR*UR);
	dUuLR[2] = xi_LR*(rhoL*VL+rhoR*VR);
	dUuLR[3] = xi_LR*(rhoL*WL+rhoR*WR);
	dUuLR[4] = xi_LR*(rhoL*HtL+rhoR*HtR);
	for(int i=0; i<YL.size()-1; ++i){
		dUuLR[5+i] = xi_LR*(rhoL*YL[i]+rhoR*YR[i]);
	}
	
    //> X. Li et al., 2017
	vector<double> dPuPp(5+YL.size()-1,0.0);	
	double Mcy = sqrt(0.5*(UL*UL+VL*VL+WL*WL+UR*UR+VR*VR+WR*WR))/chat;
	double phi_rer = 1.0-s1*(1.0-f_func(Mcy));
	dPuPp[0] = phi_rer*max(0.0,chat - Und)*rhohat*(UnR-UnL) + 
	           ( (unhat > 0.0) ? 1.0 : -1.0 )*min(Und,chat)*(PR-PL)/chat;
	dPuPp[1] = dPuPp[0]*nvec[0]; 
	dPuPp[2] = dPuPp[0]*nvec[1];
	dPuPp[4] = dPuPp[0]*nvec[2]; 
	dPuPp[0] = 0.0;
	
	
	
	// flux
	flux.clear();
	for(int i=0; i<5+(YL.size()-1); ++i){
		flux.push_back( 
			0.5*(fluxAdvL[i]+fluxAdvR[i]) 
			- 0.5*(xi*delW[i]+dPuPp[i]+dUp[i]+dUuLR[i]) );
	}
	
		
}


double SEMO_Solvers_Builder::f_func(double M){

	double mu;

	mu = min( M*sqrt(4.0+pow((1.0-M*M),2.0))/(1.0+M*M), 1.0 );
	
	return mu;
}






void SEMO_Solvers_Builder::calcFluxUpwind(
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
	
	double f1L;
	double f1R;
	if(UnL+UnR>=0.0){
		f1L=rhoL*0.5*(UnL+UnR);
		f1R=0.0;
	}
	else{
		f1L=0.0;
		f1R=rhoR*0.5*(UnL+UnR);
	}

	double PLR = 0.5*(PL+PR);
	
	// flux
	flux.clear();
	flux.push_back( f1L    + f1R );
	flux.push_back( f1L*UL + f1R*UR   + PLR*nvec[0] );
	flux.push_back( f1L*VL + f1R*VR   + PLR*nvec[1] );
	flux.push_back( f1L*WL + f1R*WR   + PLR*nvec[2] );
	flux.push_back( f1L*HtL+ f1R*HtR );
	for(int i=0; i<YL.size()-1; ++i){
		flux.push_back( f1L*YL[i] + f1R*YR[i] );
	}
	
}






/*==============================================================================
!> AUSMPW+_N Riemann solver for multiphase
!> Ref : Computational Methods for Homogeneous
!> Multi-phase Real Fluid Flows at All Speeds, 2017
!> Ref : Accurate and efficient computations of phase-changing flows in thermal
!> vapor compressors, 2018
!> Ref : Methods for Accurate Computations of Homogeneous
!> Multi-phase Real Fluid Flows at All Speeds, 2017
!> Ref : Computations of Homogeneous Multiphase Real
!> Fluid Flows at All Speeds
*/
void SEMO_Solvers_Builder::calcFluxAUSMPWP_N(
	double& PL, double& UL, double& VL, double& WL, double& TL, 
	vector<double>& YL, double& rhoL, double& CL, double& HtL, 
	double& PR, double& UR, double& VR, double& WR, double& TR, 
	vector<double>& YR, double& rhoR, double& CR, double& HtR,
	vector<double>& nvec,
	double& w1 , double& w2 ,  double& w3 , double& fC, double Uco, double dtstep, double Lch,
	vector<double>& flux){
	
	// properties of Left
	double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
	// properties of Right
	double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];

	double RT = sqrt(rhoR/rhoL);
	double rhohat = (rhoL+RT*rhoR)/(1.0+RT); //0.5*(rhoL+rhoR);
	double uhat = (UL+RT*UR)/(1.0+RT); //0.5*(uL+uR); 
	double vhat = (VL+RT*VR)/(1.0+RT); //0.5*(vL+vR); 
	double what = (WL+RT*WR)/(1.0+RT); //0.5*(wL+wR);
	//> Low-Diffusion Flux-Splitting Methods for Real Fluid Flows at All Speeds, 1999
	double chat = fC;
	double Mco = Uco/chat;
	double Mkhat = sqrt(uhat*uhat+vhat*vhat+what*what)/chat;
	double Mun = Lch/3.141592/dtstep/chat;
	//======= Scaling ============
	double thetaTemp = max(Mkhat,Mco);
	double thetaP = min(1.0,max(thetaTemp,Mun));
	double thetaU = min(1.0,thetaTemp);
	double phiP = thetaP*(2.0-thetaP);
	double phiU = thetaU*(2.0-thetaU); 
	double aSDST = 0.1875*(-4.0+5.0*phiU*phiU);
	//========================================
	double ML = UnL/chat; 
	double MR = UnR/chat;
	//========================================
	// calculate M+ and P+ for left state
	double MLP  = M_func(ML,1.0,0.125); 
	double preP = pre_func(ML,1.0,aSDST);
	// calculate M- and P- for left state
	double MRM  = M_func(MR,-1.0,0.125); 
	double preM = pre_func(MR,-1.0,aSDST);
	// ======= carefully, sensitive ===========
	double rhoLR = rhoL;
	double MLR = MLP + MRM;
	if( MLR < 0.0 ) rhoLR = rhoR;
	w1 = 1.0 - w1*w1*w1;
	w2 = 1.0 - w2*w2;
	double w = max(w1,w2);
	double fL = PL/(rhohat*chat*chat)*(1.0-w)*rhohat/rhoLR   /phiP;
	double fR = PR/(rhohat*chat*chat)*(1.0-w)*rhohat/rhoLR   /phiP;
	// =======================================
	double MBLP, MBLM;
	if(MLR >= 0.0) {
	  MBLP = MLP + MRM*((1.0-w)*(1.0+fR)-fL); 
	  MBLM = MRM * w * (1.0+fR);
	}
	else{
	  MBLP = MLP * w * (1.0+fL); 
	  MBLM = MRM + MLP*((1.0-w)*(1.0+fL)-fR);
	}
	double f1L = chat*rhoL*MBLP; 
	double f1R = chat*rhoR*MBLM;
	double Ku = 0.5;
	double pu = -2.0*Ku*preP*preM*rhohat*chat*(UnR-UnL);
	double PLR = preP*PL + preM*PR +  pu * phiU;
	
	//> comp. convective flux
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




/*==============================================================================
!> RoeM_N Riemann solver for multiphase
!> Ref : Computations of Homogeneous Multiphase Real
!> Fluid Flows at All Speeds, 2018
*/
void SEMO_Solvers_Builder::calcFluxRoeM_N(
	double& PL, double& UL, double& VL, double& WL, double& TL, 
	vector<double>& YL, double& rhoL, double& CL, double& HtL, 
	double& PR, double& UR, double& VR, double& WR, double& TR, 
	vector<double>& YR, double& rhoR, double& CR, double& HtR,
	vector<double>& nvec,
	double& w1 , double& w2 , double& w3 , double& fC, double Uco, double dtstep, double Lch,
	vector<double>& flux){

	// properties of Left
	double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
	// properties of Right
	double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
	
	double RT = sqrt(rhoR/rhoL);
	double rhohat = (rhoL+RT*rhoR)/(1.0+RT); //0.5*(rhoL+rhoR);
	double uhat = (UL+RT*UR)/(1.0+RT); //0.5*(uL+uR); 
	double vhat = (VL+RT*VR)/(1.0+RT); //0.5*(vL+vR); 
	double what = (WL+RT*WR)/(1.0+RT); //0.5*(wL+wR);
	double Hthat = (HtL+RT*HtR)/(1.0+RT); //0.5*(wL+wR);
	double Unhat = uhat*nvec[0] + vhat*nvec[1] + what*nvec[2];
	vector<double> Yihat;
	for(int i=0; i<YL.size()-1; ++i){
		Yihat.push_back((YL[i]+RT*YR[i])/(1.0+RT)); 
	} 
	double chat = fC;
	
	// double rhohat = 0.5*(rhoL+rhoR);
	// double uhat = 0.5*(UL+UR); 
	// double vhat = 0.5*(VL+VR); 
	// double what = 0.5*(WL+WR);
	// double Hthat = 0.5*(HtL+HtR);
	// double Unhat = uhat*nvec[0] + vhat*nvec[1] + what*nvec[2];
	// vector<double> Yihat;
	// for(int i=0; i<YL.size()-1; ++i){
		// Yihat.push_back(0.5*(YL[i]+YR[i])); 
	// } 
	// double chat = 0.5*(CL+CR);
	
	double Mco = Uco/chat;
	double Mkhat = sqrt(uhat*uhat+vhat*vhat+what*what)/chat;
	double Mun = Lch/3.141592/dtstep/chat;
	//======= Scaling ============
	double thetaTemp = max(Mkhat,Mco);
	double thetaP = min(1.0,max(thetaTemp,Mun));
	double thetaU = min(1.0,thetaTemp);
	double phiP = thetaP*(2.0-thetaP);
	double phiU = thetaU*(2.0-thetaU); 
	double aSDST = 0.1875*(-4.0+5.0*phiU*phiU);
	//========================================
	double Mhat = Unhat/chat;
	//========================================
	double h = 1.0 - w3;
	double f = 1.0;
    if( uhat*uhat+vhat*vhat+what*what != 0.0 ) f=pow(abs(Mhat),h);
	double g = 1.0;
    if( Mhat != 0.0 ) g = pow(abs(Mhat),1.0-w1);
	
	vector<double> fluxcL, fluxcR, delWdash, BdelW1dash, BdelW2dash;
	
    fluxcL.push_back(rhoL*UnL); 
	fluxcL.push_back(rhoL*UnL*UL + PL*nvec[0]);
	fluxcL.push_back(rhoL*UnL*VL + PL*nvec[1]);
	fluxcL.push_back(rhoL*UnL*WL + PL*nvec[2]);
	fluxcL.push_back(rhoL*UnL*HtL);
	for(int i=0; i<YL.size()-1; ++i){
		fluxcL.push_back(rhoL*UnL*YL[i]);
	}
	
    fluxcR.push_back(rhoR*UnR); 
	fluxcR.push_back(rhoR*UnR*UR + PR*nvec[0]);
	fluxcR.push_back(rhoR*UnR*VR + PR*nvec[1]);
	fluxcR.push_back(rhoR*UnR*WR + PR*nvec[2]);
	fluxcR.push_back(rhoR*UnR*HtR);
	for(int i=0; i<YL.size()-1; ++i){
		fluxcR.push_back(rhoR*UnR*YR[i]);
	}

    delWdash.push_back(rhoR - rhoL); 
    delWdash.push_back(rhoR*UR - rhoL*UL); 
    delWdash.push_back(rhoR*VR - rhoL*VL); 
    delWdash.push_back(rhoR*WR - rhoL*WL); 
    delWdash.push_back(rhoR*HtR - rhoL*HtL);
	for(int i=0; i<YL.size()-1; ++i){
		delWdash.push_back(rhoR*YR[i] - rhoL*YL[i]);
	} 
	
    BdelW1dash.push_back(rhoR - rhoL - f*(PR-PL)/chat/chat  /phiP); 
    BdelW1dash.push_back(BdelW1dash[0]*uhat + rhohat*(UR-UL)); 
    BdelW1dash.push_back(BdelW1dash[0]*vhat + rhohat*(VR-VL)); 
    BdelW1dash.push_back(BdelW1dash[0]*what + rhohat*(WR-WL)); 
    BdelW1dash.push_back(BdelW1dash[0]*Hthat + rhohat*(HtR-HtL)); 
	for(int i=0; i<YL.size()-1; ++i){
		BdelW1dash.push_back(BdelW1dash[0]*Yihat[i] + rhohat*(YR[i]-YL[i])); 
	} 

    BdelW2dash.push_back(0.0); 
    BdelW2dash.push_back(rhohat*(UnR-UnL)*(-nvec[0])); 
    BdelW2dash.push_back(rhohat*(UnR-UnL)*(-nvec[1])); 
    BdelW2dash.push_back(rhohat*(UnR-UnL)*(-nvec[2])); 
    BdelW2dash.push_back(0.0); 
	for(int i=0; i<YL.size()-1; ++i){
		BdelW2dash.push_back(0.0); 
	} 

    double b1 = max(Unhat+chat,UnR+chat); b1 = max(b1,0.0);
    double b2 = min(Unhat-chat,UnL-chat); b2 = min(b2,0.0);
	double Mtilde = ( Mhat>0.0 ? 1.0 : -1.0 ) * min(1.0,abs(Mhat));
	double b1v = max(Unhat+chat*phiU,UnR+chat*phiU); b1v = max(b1v,0.0);
	double b2v = min(Unhat-chat*phiU,UnL-chat*phiU); b2v = min(b2v,0.0);
	double Mvtilde = ( Unhat/chat/phiU>0.0 ? 1.0 : -1.0 )*min(1.0,abs(Unhat/chat/phiU));

	//> comp. convective flux
	flux.clear();
	// cout << b1 << " " << b2 << endl;
	for(int i=0; i<5+YL.size()-1; ++i){
		double fluxTemp = (b1*fluxcL[i] - b2*fluxcR[i])/(b1-b2)
			+ b1*b2/(b1-b2)*(delWdash[i]-g/(1.0+abs(Mtilde))*BdelW1dash[i])
			- b1v*b2v/(b1v-b2v)*g/(1.0+abs(Mvtilde))*BdelW2dash[i];
		flux.push_back(fluxTemp);
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




