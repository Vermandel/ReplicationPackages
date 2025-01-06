
% load steady state functions
addpath('../core_files')

//---------------------------------------------------------------------
// 1. Variable declaration
//---------------------------------------------------------------------

@#ifndef OPTIMAL
	@#define OPTIMAL = 0
@#endif

@#ifndef EZ
	@#define EZ = 1
@#endif

@#ifndef DIET
	@#define DIET = 1
@#endif

@#ifndef TaxSS
	@#define TaxSS = 0
@#endif

@#ifndef HAB
	@#define HAB = 1
@#endif

@#ifndef FULL_INFO
	@#define FULL_INFO = 0
@#endif

@#if EZ==1
	var  ce; 
@#endif

var sdf, c, kT, x, kY, kH, iT, lambda,   q,  phi,  varrho, gamma,  y,  iH, d,  aY,  e,   w,  rF,  pE, realEPquat,   gy,  gc,  giT, cyrat, rFquat, kHkYrat, SCC,c_x,z, multH, u;  

varexo eaY   ${std(\eta_{At})}$ (long_name='Productivity shock');                          


//---------------------------------------------------------------------
// 2. Parameter declaration and calibration
//---------------------------------------------------------------------

parameters nu, EIS, xbar, iotaE, iotaX, scale_y, scale_ce,  chi, estar, kappa, 
			thetaEZ  	${\theta}$ (long_name='Utility curvature'),
			mhab  		${\varpi}$ (long_name='Habit smoothing'), 
			eps 		${\epsilon}$ (long_name='Investment cost'),
			squig 		${\xi}$ (long_name='Capital in R\&D'),
			beta 		${\beta}$ (long_name='Discount factor'),
			gamma_sts 	${g}$ (long_name='Steady state growth'),
			sigma 		${\sigma}$ (long_name='Risk aversion'),
			rhoY 		${\rho_A}$ (long_name='Persistence productivity'),
			deltaK 		${\delta_K}$ (long_name='Capital depreciation'),
			Sg, nY_sts, nH_sts, th1, th2, zetaY, zetaY_tilda, deltaX,  alpha, g, std_y, xfactor, scale_u, vartheta, vtax;

// addons for Bayesian estimation
@#if FULL_INFO==1
	var aD aC aE aG ge; 
	varexo  eC  ${std(\eta_{Bt})}$ (long_name='Preference shock'),  
			eD  ${std(\eta_{Dt})}$ (long_name='Dividend shock'), 
			eE  ${std(\eta_{Et})}$ (long_name='Emission shock'), 
			eG  ${std(\eta_{Gt})}$ (long_name='Spending shock');                          
	parameters 	rhoD  ${\rho_{D}}$ (long_name='Dividend presistence'), 
				rhoC  ${\rho_{C}}$ (long_name='Preference presistence'),
				rhoE  ${\rho_{E}}$ (long_name='Emissions presistence'), 
				rhoG  ${\rho_{G}}$ (long_name='Sprending presistence'); 
	rhoD = 0; rhoC = 0.8; rhoE = 0.8; rhoG = 0.8;
	std_y 		= 1;
@#else 
	parameters  eD eC eE eG aD aC aE aG; eD = 1; eC = 1; eE = 1; eG = 1; aD = 1; aC = 1; aE = 1; aG = 1; 
	std_y 		= 0.0123;
@#endif

// addons for non estimated models
@#if DIET==0
	var U gq gd gE X pB;
	var 
	ln_SCC ln_c ln_y  Esdf muSDF varSDF skewSDF ln_pd asset_return EC1 pB  EretB bpann ln_c100 ESCC vE_irfs EPquat ln_z ln_EP ln_rF
	;
	
@#endif

//*************************

// Calibrated
gamma_sts  = .0060; 
// Optimized parameters
rhoY       = 0.9859;
eps        = 4.09595;
mhab       = 0.94081; 
sigma      = 5.98712;
thetaEZ    = 0.68716;
deltaK     = 0.00820;
squig      = 0.11649;    
beta       = 0.98444;
nu         = 0.16070;
//************************
EIS 		= 1/(1-thetaEZ);
scale_y   	= 4.3;
iotaX     	= 3/11;
Sg          = 0.2;  
alpha       = 0.33;   
deltaX 		= 0.0021;
xbar   		= 0;
iotaE  		= 0.45;
chi    		= 3; 
zetaY_tilda = 0.85;  
nH_sts   	= .029;
nY_sts  	= 0.2-nH_sts;           
vartheta    = 1;
vtax		= 0.0000;
//***********************
parameters eta1, eta2;
eta1		= 0.02;
eta2		= 2.6;
@#if OPTIMAL
	var mu vE vX;
@#else
	parameters mu vE;
@#endif
@#if HAB==0
	mhab=1;
@#endif
// scaling variable
xfactor 	= 100; 



//---------------------------------------------------------------------
// 3. Model declaration
//---------------------------------------------------------------------


model;  

// PREFERENCES
@#if EZ==1
u*scale_u       = ((1-beta)*((aC*c/((x(-1)*xfactor)^nu)-z(-1))^(thetaEZ))+beta*(gamma^thetaEZ)*((scale_ce*ce)^(thetaEZ/(1-sigma))))^(1/thetaEZ);   
scale_ce*ce      = ((u(+1)*scale_u)^(1-sigma));
lambda  = ((u*scale_u)^(1-thetaEZ))*(1/((x(-1)*xfactor)^nu))*(((1-beta)*((aC*c/((x(-1)*xfactor)^nu)-z(-1))^(thetaEZ-1))-multH*(1-mhab)));
multH   = beta*(gamma^(thetaEZ-1))*((scale_ce*ce)^(((thetaEZ-(1-sigma))/(1-sigma))))*((u(+1)*scale_u)^(1-sigma-thetaEZ))*((1-beta)*((aC(+1)*c(+1)/((x*xfactor)^nu)-z)^(thetaEZ-1))+multH(+1)*mhab);
sdf     = beta*(gamma(-1)^(thetaEZ-1))*(((u*scale_u)/((ce(-1)*scale_ce)^(1/(1-sigma))))^(1-sigma-thetaEZ))*(((x(-1)*xfactor)/(x(-2)*xfactor))^-nu)*(((1-beta)*((aC*c/((x(-1)*xfactor)^nu)-z(-1))^(thetaEZ-1))-multH*(1-mhab))/((1-beta)*((aC(-1)*c(-1)/((x(-2)*xfactor)^nu)-z(-2))^(thetaEZ-1))-multH(-1)*(1-mhab)));
@#else
u       = ((1-sigma)^-1) * ((aC*c/((x(-1)*xfactor)^nu)-z(-1))^(1-sigma))+beta*(gamma^(1-sigma))*u(+1);   
sdf     = beta*((gamma(-1))^-sigma)*lambda/lambda(-1);
lambda   = ( ((aC/((x(-1)*xfactor)^nu))^(1-sigma))*(c^-sigma) - aC/((x(-1)*xfactor)^nu)*multH*(1-mhab) ); 
multH   = beta*((gamma(-1))^-sigma)* ( multH(+1)*mhab + ((aC(+1)*c(+1)/((x*xfactor)^nu)-z)^-sigma) ) ;
@#endif

@#if HAB==1
	gamma*z = mhab*z(-1) + (1-mhab)*aC*c/((x(-1)*xfactor)^nu);
@#else
	z = 0;
@#endif


// ENDOGENOUS GROWTH
iH                  = (exp(-zetaY*((x(-1)*xfactor)-xbar)))*kappa*aY*(kH^squig)*(nH_sts^(1-squig));
phi                 =  aD(+1)*sdf(+1)*phi(+1)*(1 + (1-squig)*iH(+1)) + sdf(+1)*varrho(+1)*(1-alpha)*y(+1) ;  
gamma               = 1 + iH ; 
varrho*alpha*(y/kY) = phi*squig*(iH/kH) ;

// DIVIDENDS
w         = ((1-alpha)*y/nY_sts);
d         = ( 1 - eta1*(mu^eta2) ) * y - w*nY_sts - iT - vE*e;
varrho    = 1 - eta1*(mu^eta2) - vE*(1-mu)*iotaE;

// PRODUCTION
1         = q*th1*((iT/kT(-1))^-eps);
q         = aD(+1)*sdf(+1)*q(+1)*((1-deltaK)+(th1/(1-eps))*(iT(+1)/kT)^(1-eps)-th1*(iT(+1)/kT)^(1-eps)+th2) + sdf(+1)*varrho(+1)*alpha*(y(+1)/kY(+1)) ;
gamma*kT  = (1-deltaK)*kT(-1)+ ((th1/(1-eps))*((iT/kT(-1))^(1-eps))+th2)*kT(-1);
y         = scale_y*vartheta*(exp(-zetaY*((x(-1)*xfactor)-xbar)))*aY*(kY^alpha)*(nY_sts^(1-alpha));
kT(-1)    = kH+kY;

// EMISSIONS
x*xfactor  	= (1-deltaX)*x(-1)*xfactor + iotaX*(e+estar);
e  			= aE*(1-mu)*iotaE*y;

// MARKET CLEARING
y* (1-eta1*(mu^eta2)) =  c + iT + g*aG;

// EXOGENOUS PROCESSES
log(aY)   = rhoY*log(aY(-1)) + std_y*eaY;
@#if FULL_INFO==1                      
log(aD)   = rhoD*log(aD(-1)) + eD;
log(aC)   = rhoC*log(aC(-1)) + eC;
log(aE)   = rhoE*log(aE(-1)) + eE;
log(aG)   = rhoG*log(aG(-1)) + eG;
@#endif

// ASSET PRICING
pE       = sdf(+1)*gamma*aD(+1)*(d(+1)+pE(+1));
1/(1+rF) = sdf(+1);
realEPquat = (aD*gamma(-1)*(pE+d)/pE(-1)-(1+rF))*100;
rFquat = rF*100;

// MACRO
cyrat 	= 100*c/y;
kHkYrat = 100*kH/kY;
gy      =  (log(gamma(-1))+ log(y)-log(y(-1)))*100;
gc   	=  (log(gamma(-1))+ log(c)-log(c(-1)))*1;
giT  	=  (log(gamma(-1))+ log(iT)-log(iT(-1)))*100;
c_x 	=  1000*nu*c/(x(-1)*xfactor);
@#if FULL_INFO==1                      
ge   	=  (log(e)-log(e(-1)))*100;
@#endif

// CLIMATE
@#if TaxSS==0
	SCC        = gamma*sdf(+1)*( iotaX*1000*(vtax+nu*c(+1)/(xfactor*x) + zetaY*varrho(+1)*y(+1) +  zetaY*phi(+1)*iH(+1)) + SCC(+1)*(1-deltaX) ) ;
@#else
	SCC        = steady_state(SCC);
@#endif
@#if OPTIMAL
	vE*iotaE  = eta2*eta1*(mu^(eta2-1)); 
	vX*iotaX  = vE; 
	@#if TaxSS==0
		vX        = gamma*sdf(+1)*( vtax+nu*c(+1)/(xfactor*x) + zetaY*varrho(+1)*y(+1) +  zetaY*phi(+1)*iH(+1) + vX(+1)*(1-deltaX) ) ;
	@#else
		vX        = steady_state(vX);
	@#endif
@#endif

@#if DIET==0                      
U  				= u*scale_u;
gq   			= (log(q)-log(q(-1)))*100; 
gd   			= (log(gamma(-1))+ log(d)-log(d(-1)))*100;
gE   			= (log(e)-log(e(-1)))*100;
X 				= x*xfactor;

asset_return 	= (aD*gamma(-1)*(pE+d)/pE(-1))-1;
EPquat 	= (aD(+1)*gamma*(pE(+1)+d(+1))/pE-(1+rF))*100;
Esdf	 = sdf(+1);
muSDF	 = sdf(+1);
varSDF	 = (sdf(+1)-muSDF)^2;
skewSDF	 = (sdf(+1)-muSDF)^3;
ln_pd	 = log(pE/d);
EC1		= (SCC/c-steady_state(SCC)/steady_state(c))/(steady_state(SCC)/steady_state(c));
pB   	= sdf(+1)/gamma*(1+pB(+1));
EretB 	= (1+pB(+1))/pB;
bpann 	= (EretB-(1+rF))*100;
ESCC  =  SCC(+1);
ln_SCC  = log(SCC/steady_state(SCC));
ln_c  	= log(c/steady_state(c));
ln_y  	= log(y/steady_state(y));
ln_c100	= ln_c*100; 
ln_EP   = log((aD(+1)*gamma*(pE(+1)+d(+1))/pE)/(1+rF));
ln_rF   = log(1+rF);
@#if HAB==1
	ln_z    = log(z/steady_state(z));
@#else
	ln_z = 0;
@#endif

@#if OPTIMAL
	vE_irfs = log(vE);
@#else
	vE_irfs = 0;
@#endif
@#endif

end;


//---------------------------------------------------------------------
// 4. Initial values and steady state
//---------------------------------------------------------------------


steady_state_model; 
@#if EZ == 0
	thetaEZ = 1-sigma;
	scale_u     = 1;
	scale_ce    = 1;
@#else
	@#if HAB==1
		scale_u     = 1;
	@#else
		scale_u     = (1-beta);
	@#endif 
@#endif  
gamma_lf 	= 1+gamma_sts;
btilda_lf   = beta*(gamma_lf^(thetaEZ-1));
q      		= 1;
varrho_lf 	= 1;
kYyratio_lf = varrho_lf*btilda_lf*alpha/(q*((1-btilda_lf*(1-deltaK))));
y_lf      	= ((scale_y)^(1/(1-alpha)))*((zetaY_tilda)^(1/(1-alpha)))*((kYyratio_lf)^(alpha/(1-alpha)))*nY_sts; 
yxrat_lf    = deltaX/( iotaX * (1+chi)*iotaE  ) ;
x_lf     	= (1/yxrat_lf)*y_lf/xfactor;
%zetaY   	= -log(zetaY_tilda)/((xfactor*x_lf)-xbar);
zetaY   	= -log(zetaY_tilda)*0.0021/( iotaX * (1+chi)*iotaE * y_lf );

%zetaY_tilda = exp(-zetaY*( iotaX * (1+chi)*iotaE * y_lf )/deltaX);
e_lf      	= iotaE*y_lf;
phi_lf    	= varrho_lf*(btilda_lf*(1-alpha)*y_lf)/(1-btilda_lf*(1+(1-squig)*(gamma_lf-1)));
kHiHrat_lf  = (phi_lf/varrho_lf)*(squig/alpha)*kYyratio_lf; 
iH_lf     	= gamma_lf-1;
kH_lf     	= kHiHrat_lf*iH_lf;
kappa      	= iH_lf/((zetaY_tilda)*(kH_lf^squig)*(nH_sts^(1-squig)));
@#if OPTIMAL
	[gamma,mu,zetaY_tilda_now]   =  model_file_solve(deltaX,eta1,eta2,beta,sigma,squig,alpha,nY_sts,nH_sts,deltaK,zetaY,kappa,scale_y*vartheta,iotaX,iotaE,thetaEZ,gamma_lf,y_lf,x_lf,phi_lf,iH_lf,Sg,nu,xfactor,chi,vtax);
	varrho      = 1*(1-(eta1)*(mu^eta2))-eta1*eta2*(mu^(eta2-1))*(1-mu);
	btilda    	= beta*(gamma^(thetaEZ-1));
	kYyratio   	= varrho*btilda*alpha/(q*((1-btilda*(1-deltaK))));
	y      		= ((scale_y*vartheta)^(1/(1-alpha)))*((zetaY_tilda_now)^(1/(1-alpha)))*((kYyratio)^(alpha/(1-alpha)))*nY_sts; 
	iH     		= gamma-1;
	phi    		= varrho*(btilda*(1-alpha)*y)/(1-btilda*(1+(1-squig)*iH));
	e  			= (1-mu)*iotaE*y;
    estar      	= chi*e;
	x  			= iotaX/(xfactor*deltaX)*((1-mu)*iotaE*y+estar);
	vE  		= eta2*eta1*(mu^(eta2-1))/iotaE;
	vX		  	= vE/iotaX; 
	kHiHrat    	= (phi/varrho)*(squig/alpha)*kYyratio; 
	kH     		= kHiHrat*iH;
@#else
	mu 		= 0;
	vE 		= 0;
	varrho	= varrho_lf;
	%% applies welfare cost calculation if vartheta ~=1
	[y, x, zetaY_tilda_lf, phi, kH, gamma] = model_ss_lf_welfcost(scale_y,zetaY,alpha,xfactor,deltaX,vartheta,kYyratio_lf,kappa,squig,nY_sts,nH_sts,beta,thetaEZ,y_lf, x_lf, zetaY_tilda, phi_lf, kH_lf, gamma_lf,iotaX,iotaE,chi,deltaK);	
	btilda    	= beta*(gamma^(thetaEZ-1));
	kYyratio   	= varrho*btilda*alpha/(q*((1-btilda*(1-deltaK))));
	iH     		= gamma-1;
	e  			= iotaE*y;
    estar      	= chi*e;
@#endif
kY     		= kYyratio*y;
kT     		= kY+kH;
iT     		= (gamma-(1-deltaK))*kT;
g      		= Sg*y; 
w      		= (1-alpha)*(y/nY_sts);
d      		= (1-eta1*(mu^eta2))*y-w*(nY_sts)-iT- vE*e;
c      		= y*(1-eta1*(mu^eta2))-iT-g;
@#if HAB==1
	z      	= ((1-mhab)/(gamma-mhab))*(c/((xfactor*x)^nu));
@#else
	z = 0;
@#endif
@#if EZ==1
	multH  		= (btilda*(1-beta)/(1-btilda*mhab))*((c/((xfactor*x)^nu)-z)^(thetaEZ-1));
	u      		= scale_u^-1*(((1-beta)/(1-beta*(gamma^thetaEZ)))^(1/thetaEZ))*(c/((xfactor*x)^nu)-z);
	ce          = 1;scale_ce     		= (u*scale_u)^(1-sigma);
	lambda 		= ((u*scale_u)^(1-thetaEZ))*(1/((xfactor*x)^nu))*(((1-beta)*((c/((xfactor*x)^nu)-z)^(thetaEZ-1))-multH*(1-mhab)));
	sdf    		= btilda;
@#else
	sdf     	= beta*((gamma)^-sigma);
	u       	= ((1-sigma)^-1)*(((c/((x*xfactor)^nu)-z)^(1-sigma)))/(1-beta*(gamma^(1-sigma)));   
	multH  		= ((c/((x*xfactor)^nu)-z)^-sigma)/(1/btilda-mhab) ;
	lambda   	= ( ((1/((x*xfactor)^nu))^(1-sigma))*(c^-sigma) - 1/((x*xfactor)^nu)*multH*(1-mhab) ); 
@#endif
rF     		= 1/(btilda)-1;
pE     		= (btilda*gamma/(1-btilda*gamma))*d;
th1        	= (gamma-1+deltaK)^(eps);                   
th2        	= -(gamma-1+deltaK)*(eps/(1-eps)); 
aY  		= 1;
EretE 		= gamma*(pE+d)/pE;
EPquat 		= (gamma*(pE+d)/pE-(1+rF))*100;
EPann 		= (EretE-(1+rF))*100;
rFann 		= rF*100;
pB   		= 1/(gamma/sdf-1);
EretB 		= (1+pB)/pB;
bpann 		= (EretB-(1+rF))*100;
rFquat 		= rF*100;
realEPquat 	= (gamma*(pE+d)/pE-(1+rF))*100;
vX      	= ( vtax+nu*c/(xfactor*x) + zetaY*varrho*y +  zetaY*phi*iH  )/(1/(gamma*sdf)-(1-deltaX));
gy   		= (log(gamma))*100;
gc   		= (log(gamma))*1;
giT  		= (log(gamma))*100;
cyrat    	= 100*c/y;
kHkYrat 	= 100*kH/kY;
SCC			= 1000*iotaX*(vX);
c_x 		= 1000*nu*c/(x*xfactor);
@#if FULL_INFO==1                      
	aD = 1; aC = 1; aE = 1; aG = 1;
	ge = 0;
@#endif

@#if DIET==0                      
U  = u*scale_u; gE = 0; gq = 0; gd=log(gamma)*100; X 		= x*xfactor;
ESCC     		= SCC;
Esdf	 		= sdf;
muSDF			= Esdf;
varSDF			= (Esdf-muSDF)^2; skewSDF	 = (Esdf-muSDF)^3;
ln_pd			= log(pE/d);
asset_return	= (aD*gamma*(pE+d)/pE-1);
ln_SCC = 0; ln_c = 0; ln_y = 0; ln_c100	= ln_c*100;  ESCC  =  SCC;
ln_z    = 0;
ln_EP   = log((aD*gamma*(pE+d)/pE)/(1+rF));
ln_rF   = log(1+rF);
@#if OPTIMAL
	vE_irfs = log(vE);
@#else
	vE_irfs = 0;
@#endif
@#endif

EC1			= 0;


end;

resid;

/*
	set initial state variables for simulations
*/
histval;
	x(0)  = 60/xfactor;
end;

steady;

//---------------------------------------------------------------------
// 5. Shock declaration                       
//---------------------------------------------------------------------

shocks;
var eaY; stderr 1;        
end;

stoch_simul(order=3,pruning,noprint,nomoments,nofunctions,irf=0,ar=0,nocorr);
