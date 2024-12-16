
addpath('estimation_v6')


%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------



% MA (GtC)
% GtC = (3/11) × GtCO2
% emission are GtCO2
% 1ppm CO2 = 2.13 GtC
% 1961 -> 318 ppm CO2 = 2.13*318
close all;


@#ifndef SIMULATIONS
	@#define SIMULATIONS = 0
@#endif



%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

 
var Z   		${Z_t}$ (long_name='TFP'),  
	L   		${L_t}$ (long_name='population (billion)'),
	gZ  		${g_{z,t}}$ (long_name='TFP growth'),  
	gL  		${g_{l,t}}$ (long_name='population growth'),  
	SIG 		${\Sigma_t}$ (long_name='Decoupling rate '),
	gSIG 		${g_{\sigma,t}}$ (long_name='Decoupling rate '),
	THETA1 		${\theta_{1,t}}$ (long_name='Abatement trend'),
	delthet  	${\delta_{p,t}}$ (long_name='Abatement depreciation factor'),
 	dy  		${g_{y,t}}$ (long_name='GDP growth'),  
	de  		${g_{y,t}}$ (long_name='Emissions growth'),
	pi  		${\pi_t}$ (long_name='Headline inflation'),  
	r 			${r_t}$ (long_name='Interest rate (pp)'),  
	y 			${\hat{y}_t}$ (long_name='Detrended output'), 
	x 			${{x}_t}$ (long_name='cons to gdp'), 
    M  			${m_t}$ (long_name='Carbon stock (Gt) '),
	mc 			${mc_t}$ (long_name='marginal cost'),
	y_n 		${\hat{y}^n_t}$ (long_name='Natural GDP'), 
	x_n 		${\hat{\lambda}^n_t}$ (long_name='Natural MUC'),  
	rr_n 		${\hat{r}^n_t}$ (long_name='Natural rate'),  
	pi_bar 		${r^{*}_{t}}$ (long_name='Inflation target'), 
	D 			${D_{t}}$ (long_name='Transfert policy'), 	 
	s_r 		${\epsilon^R_{t}}$ (long_name='Monetary policy shock'), 	 
	s_p 		${\epsilon^P_{t}}$ (long_name='Inflation shock')
	s_b 		${\epsilon^B_{t}}$ (long_name='Preference shock'),
	s_e			${\epsilon^E_{t}}$ (long_name='Emissions shock')	;


varexo 	e_b		${\sigma_b}$ 		(long_name='Std demand'),
		e_p		${\sigma_p}$ 		(long_name='Std price'), 
		e_r		${\sigma_{r}}$ 		(long_name='Std MPR'), 
		e_e		${\sigma_e}$ 		(long_name='Std emissions')
			; 

parameters   
	M0			$m_{t_{0}}$ 				(long_name='Initial stock of carbon (GtC)'),  
	xi 			$\xi_{m}$ 					(long_name='Marginal atmospheric retention ratio'),  
	LT  		$L_{\infty}$ 				(long_name='Terminal population (billion)'),  
	GZ1  		${delta_{z}}$ 				(long_name='Decay TFP'),  
	Z0 			${Z_{t_{0}}}$ 				(long_name='Initial TFP'),  
	lg			${l_{g}}$ 					(long_name='Population growth'),  
	rho			${\rho}$ 					(long_name='MPR smoothing'), 
	phi_pi		${(\phi_{\pi}-1)}$ 			(long_name='Inflation stance'), 
	sigmaC		${\sigma_{c}}$ 				(long_name='Risk aversion'),
	sigmaL		${\sigma_{h}}$ 				(long_name='Labor disutility'),
	beta 		${\beta}$ 					(long_name='Discount factor'),
	alpha 		${\alpha}$ 					(long_name='Labor intensity'), 
	chi  		${\chi}$ 					(long_name='Labor disutility'),
	PI 			${\pi^{\star}_{\infty}}$ 	(long_name='Long-term inflation target'),
	kappa		${\kappa}$ 					(long_name='Rotemberg Cost'),
	b			${(\beta^{-1}-1)\times 100}$ 	(long_name='Discount rate'),
	varsigma	${\varsigma}$ 				(long_name='Goods substitution elasticity'),
	omega 		${\omega}$ 					(long_name='Share low productive workers'),
	Dc 			${d/c}$ 					(long_name='Low productivity worker payoff-to-consumption'),
	E0  		${e_{t_{0}}}$ 				(long_name='Initial emissions (GtCO$_{2}$)'),  
	Y0  		${y_{t_{0}}}$ 				(long_name='Initial GDP (trillion USD PPP)'),  
	L0  		${l_{t_{0}}}$ 				(long_name='Initial population (billion)'),   
	rr0  		${r_{t_{0}}}$ 				(long_name='Initial Interest rate'),   
	tau0   		$\tau_{t_{0}}$ 				(long_name='Initial tax'),  
	delthet0  	$\delta_{\theta,t_{0}}$ 	(long_name='Initial decay'),   
	h0  		$h_{t_{0}}$ 				(long_name='Initial hours worked'),    
	mu0   		$\mu_{t_{0}}$ 				(long_name='Initial abatement share'), 
	THETA0   	$\theta_{1,t_{0}}$ 			(long_name='Initial abatement cost-to-gdp'),  
	pi0   		$\pi_{t_{0}}$ 				(long_name='Initial inflation'),  
	c0   		$c_{t_{0}}$ 				(long_name='Initial consumption'),  
	y0   		$y_{1,t_{0}}$ 				(long_name='Initial output'),  
	gamma   	$\gamma$ 					(long_name='Climate damage elasticity'),   
	theta2   	$\theta_{2}$ 				(long_name='Abatement cost curvature'),    
	pb			${p_{p}}$ 					(long_name='Cost abatement'),   
	deltapb		${\delta_{pb}}$ 			(long_name='Decay abatement cost'),  
	M_1750		${m_{1750}}$ 				(long_name='Pre-industrial stock of carbon (GtC)'), 
	delta_M		${\delta_{m}}$ 				(long_name='CO$_{2}$ rate of transfer to deep oceans'), 
	gZ0   		$g_{z,t_{0}}$ 				(long_name='Initial TFP growth'),  
	gS0   		$g_{\sigma,t_{0}}$ 			(long_name='Initial decoupling rate'), 
	SIG0 
	GS1 
	rr0_400		${r_{t_{0}} \times 400}$ 	(long_name='Initial interest rate'), 
	gZ0_400		${g_{z,t_{0}} \times 400}$ 			(long_name='Initial TFP growth'), 
	gS0_400		${g_{\sigma,t_{0}}}$ 		(long_name='Decay rate decoupling'), 
	GS1_400		${\delta_{\sigma}}$ 		(long_name='Decay rate emission intensity'), 
	PI_400		${\pi_^{\star}{\infty} \times 400}$ 			(long_name='Inflation ss'), 
	rho_b		${\rho_b}$ 					(long_name='AR demand'),
	rho_p		${\rho_b}$ 					(long_name='AR price'),
	rho_e		${\rho_e}$ 					(long_name='AR emissions'),
	phi_y 		${\phi_y}$ 					(long_name='MPR GDP stance'),
	phi_dy 		${\phi_y}$ 					(long_name='MPR GDP stance'),
	rho_r 		${\rho_r}$ 					(long_name='AR MPR'), 
	MP_R		${\phi_{r}}$ 				(long_name='Penalty R OSR'),
	MP_PI		${\phi_{\pi}}$ 				(long_name='Penalty inflation'),
	tf 
	nu 			${\nu}$		 				(long_name='Firm exit shock'),
	Et			${\varphi}$ 				(long_name='Mitigation policy belief'),
	PI0 		
	PI0_400		${\pi^{\star}_{t_{0}} \times 400}$ 	 (long_name='Initial inflation trend (annualized)'),
	delta_pi_star 		${\rho_{\pi^{\star}} }$ 	(long_name='Decline rate inflation target'),  
    phi_star	${\phi_{*}}$ 				(long_name='Trend stance')
	MP_y_n
	MP_r_n
	Et_rn
	yy0,
	tau0_USD		${\tau_{t_{0}}}$  (long_name='Initial carbon price (\$/ton)')
	mc_no_damage ygap MP_transfer
	
	;

 
@#if SIMULATIONS == 1
	varexo e_tau;
	var 	
	mu 					${\mu_t}$ (long_name='Abatement share'),
	c 					${\hat{c}_t}$ (long_name='Consumption (detrended)'), 
	lb 					${\lambda_t}$ (long_name='Marginal utility consumption'),
	rotemberg_to_gdp 	${nc_t}$ (long_name='Nominal costs'),
	abatement_to_gdp 	${ac_t}$ (long_name='Abatement costs'),
	tau_USD  			${\tau_t}$ (long_name='carbon tax (USD/ton)'), 
	h 					${h_t}$ (long_name='hours worked') ,
	d  					${d_t}$ (long_name='damage factor'),
	tau 				${\tau_t}$ (long_name='carbon tax')
	E  					${e_t}$ (long_name='emissions (Gt)'),
	rr_real 			${rr_{t}}$ (long_name='Real rate'),  	 
	pi100 				${\pi_{t}}$ (long_name='Inflation (\%)'),  	 
	rr100 				${r_{t}}$ (long_name='Interest rate (\%)'),  	 
	rr_real100,
 	output_gap ${log(y_{t}/y^{n}_{t})}$ (long_name='Output Gap')   
	lny_n 	${\hat{y}^n}$ 				(long_name='Natural detrended GDP')
	lb_n  
	c_n
	w
	WELF 		${\W_t}$ (long_name='welfare'),
	WELF_n 		${\W^n_t}$ (long_name='welfare')
	SCC 		${SCC_t}$ (long_name='Social Cost of Carbon'), 
	SCC_n 		${SCC^{\star}_t}$ (long_name='Social Cost of Carbon (natural)'),
	IS 			${IS_t}$ (long_name='Consumption'),
	rrn100 		${r^{\star}_t}$ (long_name='Natural Rate')
	;
	var lny 	${\hat{y}}$ 				(long_name='Detrended GDP'),
        lnc 	${\hat{c}}$ 				(long_name='Detrended cons'),%		
		mc_hat mc_hat_s mc_hat_m mc_hat_g
		pi_tot 		${\hat{\pi}}$ 				(long_name='Inflation total'),
		pi_hat_s 	${\hat{\pi}^s}$ 			(long_name='Inflation standard'),
		pi_hat_m 	${\hat{\pi}^m}$ 			(long_name='Climateflation'),
		pi_hat_g 	${\hat{\pi}^g}$ 			(long_name='Greenflation'),
		pi_hat_x 	${\hat{\pi}^x}$ 			(long_name='Inflation exogenous')
		;
	parameters MP_taylor_target_m MP_taylor_target_g ;
	MP_taylor_target_m = 0;  MP_taylor_target_g = 0; 

@#else
	varexo e_tau;
@#endif

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
% shock process
rho_p   	= 0.95; 				% productivity 
rho_b   	= 0.95; 				% productivity 
rho_e		= 0;
rho_r		= 0.2;

%% Trends
L0			= 4.85;%3.8;4.38;4.85; 	
LT			= 10.48;
lg			= 0.025/4;
gZ0_400 	= 1.70;
%GZ1_400 	= 0.3;
GZ1			= 0.0072/4;

gS0_400 	= 0.8;
GS1_400		= 0;
rr0_400		= 5;
%% initial values
Y0			= 30/4; % 1970: 18  ; 1980: 26.4 ; 1985: 30
E0			= 20.30/4; % 1980:19.48/4
mu0			= 0.0001;
M0			= 338*2.13;% 1980:338*2.13 1985:346*2.13;
h0			= 1;
b			= 0.85;
alpha 		= .7;
sigmaL		= 1/2;
kappa		= 28;
varsigma	= 4;
rho 		= 0.75;
phi_pi		= 1.5;
phi_y		= 0.09;
% Climate parameters
gamma		= 1.5*2.379e-5;

% abatement sector
%deltapb		= 1-(1-0.025)^(1/(4*5));
%pb2019		= 544/0.9287;
%pb			= pb2019/(1-deltapb)^(4*(2016-1961));	
deltapb		= 1-(1-0.017)^(1/4);
delthet0	= 1;
theta2		= 2.6;

delta_M		= 0;
xi			= 3/11;
M_1750  	= 545;


kappa		= 102.72410808;
rho			= 0.63370700;
phi_pi		= .80776810;
phi_y		= 0.20414580;
Dc			= 0.97;
b			= .76613436;
omega		= 0.02;
y0          = 1.0;
rho_b 		= 0.99170991 ;
rho_p 		= 0.86668108 ;
gZ0_400 	= 1.79151844;
sigmaC 		= 1.64088274  ;
sigmaL 		= .33875988 ;
tf			= 1;
nu			= 0.05;
MP_R		= 0;
MP_PI		= 0;
Et			= 0.75;
phi_dy 		= 0;
PI0_400     = 10;
PI_400		= 1;
phi_star    = 0;
MP_y_n    	= 1;
MP_r_n    	= 0;
Et_rn		= 0;
sigmaC 		= 2.5;
delta_pi_star = 7.5/400;
MP_transfer = 0;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model;
 	%% TREND BLOCK
	pi_bar = (1-delta_pi_star)*pi_bar(-1) + delta_pi_star*steady_state(pi_bar);
	[name='productivity trend']
	Z 		= Z(-1)*(1+gZ(-1));
	gZ      = gZ(-1)*(1-GZ1);
	[name='population trend']
	L 		= L(-1)^(1-lg)*LT^lg;	
	gL      = L/L(-1)-1;
	[name='emissions trend']
	SIG 	= SIG(-1) * (1-gSIG(-1));
	gSIG    = (1-GS1)*gSIG(-1);
	[name='cost of abatement - level']
	THETA1 	= max(pb/1000/theta2*delthet*SIG,0);
	[name='efficiency trend']
	delthet = delthet(-1)*(1-deltapb);
	%% SHOCKS
	[name='shocks']
	s_b = 1-rho_b + rho_b*s_b(-1) + e_b;
	s_p = 1-rho_p + rho_p*s_p(-1) + e_p;
	s_e = 1-rho_e + rho_e*s_e(-1) + e_e;
	s_r = 1-rho_r + rho_r*s_r(-1) + e_r;
	%% ENDOGENOUS VARIABLES
 	[name='IS']
	((x*y/(1-omega))-omega*D/(1-omega))^(-sigmaC) = s_b(+1)/s_b*beta*(
		   (1-omega)*((x(+1)*y(+1)/(1-omega))-omega*D(+1)/(1-omega))^(-sigmaC)
		 + omega*D(+1)^-sigmaC  
		    ) * ((1+r))/(1+pi(+1));
	[name='Aux:']
	x = (1-(1-nu)*kappa/2*(pi-pi_bar)^2-tf*THETA1*(tau0+Et*e_tau)^(theta2/(theta2-1))-nu*(1-s_p*mc));
	[name='PC']	
	(1+pi)*(pi-pi_bar) = (1-varsigma)/kappa + s_p*varsigma/kappa*mc +  (1-nu)*beta*(1+gZ(+1))*(y(+1)/y)*(pi(+1)-pi_bar(+1))*(1+pi(+1));
	[name='CC']
	M - M_1750 = (1-delta_M)*(M(-1) - M_1750) + xi*(1-(tau0+Et*e_tau)^(1/(theta2-1)))*SIG*y*Z*L*s_e;
	[name='Aux:mc']
	mc = (chi/(alpha*(1-omega)^sigmaL))*
		 ((x*y/(1-omega))-omega*D/(1-omega))^(sigmaC) * y^((1+sigmaL)/alpha-1) / (exp(-gamma*(M(-1)-M_1750)))^((1+sigmaL)/alpha)
		 + tf*THETA1*(tau0+Et*e_tau)*(theta2+(1-theta2)*(tau0+Et*e_tau)^(1/(theta2-1))) ;
 	[name='MP']
    (1+r)	/	 (   (1-MP_r_n)*(1+steady_state(r)) + MP_r_n*(1+rr_n)*(1+pi_bar) ) = 
				((1+r(-1))		/	(   (1-MP_r_n)*(1+steady_state(r)) + MP_r_n*(1+rr_n(-1))*(1+pi_bar(-1)) )    )^rho   
				*( (((1+pi)/((1+pi_bar)))^((1+phi_pi)) * (y/((1-MP_y_n)*y0+MP_y_n*y_n))^(phi_y)))^(1-rho) *    ((1+pi_bar)/(1+steady_state(pi_bar)))^(1-rho) * ((1+pi_bar)/(1+pi_bar(-1)))^phi_star *  s_r*ygap;
	[name='Transfert policy']
	D 		= exp(-gamma*(M(-1)-M_1750))*Dc*steady_state(x)*steady_state(y) + MP_transfer * ( tau0 + Et*e_tau)*tf*THETA1/(SIG+1e-8) * ((1-(tau0+Et*e_tau)^(1/(theta2-1)))*SIG*s_e+1e-8);
	[name='Output growth']
	dy = log((1+gZ)*(1+gL)*y/y(-1));
	[name='Emissions growth']
	de = log((1-(tau0+Et*e_tau)^(1/(theta2-1)))*SIG*y*Z*L*s_e+1e-8) - log((1-(tau0+Et*e_tau(-1))^(1/(theta2-1)))*SIG(-1)*y(-1)*Z(-1)*L(-1)*s_e(-1)+1e-8);
	

	% natural block
	[name='Natural GDP']
	steady_state(mc) = 
		(chi/(alpha*(1-omega)^sigmaL))*
		((x_n*y_n/(1-omega))-omega*D/(1-omega))^(sigmaC) * y_n^((1+sigmaL)/alpha-1) / (exp(-gamma*(M(-1)-M_1750)))^((1+sigmaL)/alpha)
		+tf*THETA1*(tau0+Et*e_tau)*(theta2+(1-theta2)*(tau0+Et*e_tau)^(1/(theta2-1))) 
		;
	[name='Cons-to-GDP natural economy']
	x_n = (1-tf*THETA1*(tau0+Et*e_tau)^(theta2/(theta2-1))-nu*(1-steady_state(mc)));
	((x_n*y_n/(1-omega))-omega*D/(1-omega))^(-sigmaC) = beta*(
		   (1-omega)*((x_n(+1)*y_n(+1)/(1-omega))-omega*D(+1)/(1-omega))^(-sigmaC)
		 + omega*D(+1)^-sigmaC  
		    ) * ((1+rr_n));


	@#if SIMULATIONS == 1
		rotemberg_to_gdp 	= 100*((1-nu)*kappa/2*(pi-pi_bar)^2+nu*(1-s_p*mc));
		IS					= ln(c/y0)*100;
		abatement_to_gdp 	= 100*tf*THETA1*mu^theta2;
		rr_real 			= (1+r)/(1+pi(+1))-1;
		
		tau_USD		= tau*tf*THETA1/(SIG+1e-8)*1000;
		pi100		= 100*pi;
		rr100 		= 100*r;
		rrn100 		= 100*rr_n;
		rr_real100 	= 100*rr_real;
		lny_n       = 100*ln(y_n/y0) ;

		[name='Optimal abatement']
		tau = tau0 + Et*e_tau; % where: tau = tau_level*SIG/THETA1
		[name='technology']
		y 	= d*h^alpha;
		[name='Emissions']	
		d = exp(-gamma*(M(-1)-M_1750));
		[name='Optimal abatement']
		mu 	= ((tau0 + Et*e_tau)*s_e)^(1/(theta2-1));
		[name='Emissions']	
		E = (1-mu)*SIG*y*Z*L*s_e;
		[name='Ressources Constraint']
		y*x = c;
		y_n*x_n= c_n;

		[name='Welfare']
		WELF   	= s_b*(Z)^(1-sigmaC)*(L)*(1/((1-sigmaC))*c  ^(1-sigmaC)-chi/(1+sigmaL)*(                h/(1-omega))^(1+sigmaL)) - MP_R*(r/r(-1)-1)^2  + (Z)^(1-sigmaC)*(L)*(  1*MP_PI*(y-y_n)^2 + 0*MP_PI*(pi-pi_bar)^2 )  /*+ (Z)^(1-sigmaC)*(L)*MP_PI*(pi-pi_bar)^2*/ + beta*WELF(+1);
		WELF_n 	=     (Z)^(1-sigmaC)*(L)*(1/((1-sigmaC))*c_n^(1-sigmaC)-chi/(1+sigmaL)*((y_n/d)^(1/alpha)/(1-omega))^(1+sigmaL)) - MP_R*(r/r(-1)-1)^2                                                                         /*+ (Z)^(1-sigmaC)*(L)*MP_PI*(pi-pi_bar)^2*/ + beta*WELF_n(+1);
	
		
		mc_hat   = 100*(mc-mc_no_damage)/mc_no_damage;
		mc_hat_s = 100*((chi/(alpha*(1-omega)^sigmaL))*((x*y/(1-omega))-omega*Dc*steady_state(c)/(1-omega))^(sigmaC) * y^((1+sigmaL)/alpha-1) - mc_no_damage)/mc_no_damage;
		mc_hat_g = 100*(tf*THETA1*(tau0+Et*e_tau)*(theta2+(1-theta2)*(tau0+Et*e_tau)^(1/(theta2-1))))/mc_no_damage;
		mc_hat_m = 100*((chi/(alpha*(1-omega)^sigmaL))*((x*y/(1-omega))-omega*D/(1-omega))^(sigmaC) * y^((1+sigmaL)/alpha-1) / (exp(-gamma*(M(-1)-M_1750)))^((1+sigmaL)/alpha) - (chi/(alpha*(1-omega)^sigmaL))*((x*y/(1-omega))-omega*Dc*steady_state(c)/(1-omega))^(sigmaC) * y^((1+sigmaL)/alpha-1) )/mc_no_damage;

		w 		 = (chi*((y/exp(-gamma*(M(-1)-M_1750)))^(1/alpha)/(1-omega))^sigmaL*((x*y/(1-omega))-omega*(D)/(1-omega))^(sigmaC));
	
		pi_tot   	=  pi_hat_x + pi_hat_s + pi_hat_m + pi_hat_g;
		pi_hat_m 	= (varsigma-1)/kappa*mc_hat_m 								+ (1-nu)*beta*(1+gZ(+1))*y(+1)/y*pi_hat_m(+1);
		pi_hat_s 	= 100*(pi-pi_bar)  - pi_hat_x - pi_hat_m - pi_hat_g;
		pi_hat_g 	= (varsigma-1)/kappa*mc_hat_g 								+ (1-nu)*beta*(1+gZ(+1))*y(+1)/y*pi_hat_g(+1);
		pi_hat_x 	= 100*(varsigma-1)/kappa*(s_p*mc-mc)/mc_no_damage    	+ (1-nu)*beta*(1+gZ(+1))*y(+1)/y*pi_hat_x(+1);
		
		lnc 		= 100*ln(c/steady_state(c)) 	;
		lny 		= 100*ln(y/y0) 	;
		output_gap  = 100*ln(y/y_n) 	;
		lb  		= s_b*((c/(1-omega))-omega*(D)/(1-omega))^(-sigmaC);
		lb_n  		= ((c_n/(1-omega))-omega*(D)/(1-omega))^(-sigmaC);
			
		[name='Social Cost of Carbon']
		SCC   		= beta*((1+gZ)^-sigmaC*lb(+1)/lb*SCC(+1)     +  1000*3/11*gamma/alpha*chi*Z(+1)*L(+1)*(  (y(+1)/d(+1))^(1/alpha)/(1-omega))^(1+sigmaL) * (c)^sigmaC );
		SCC_n 		= beta*((1+gZ)^-sigmaC*lb_n(+1)/lb_n*SCC_n(+1) +  1000*3/11*gamma/alpha*chi*Z(+1)*L(+1)*((y_n(+1)/d(+1))^(1/alpha)/(1-omega))^(1+sigmaL) * (c_n)^sigmaC );
	@#endif

end;


%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
endval;
	Z 		= Z0;
	L		= L0;
	y		= yy0;
	gZ		= gZ0;
	SIG		= SIG0;
	gSIG	= gS0;
	M 		= M0;
	r		= rr0;
	pi		= PI;
	delthet	= delthet0;
	THETA1	= THETA0;
	y_n    = y0;
	pi_bar = PI0;
	@#if SIMULATIONS
%	c		= c0;
	@#endif
end;


steady_state_model;	
	% <---- Parameters Rescaling -----> %
	beta		= 1/(1+b/100);
	gZ0			= gZ0_400/400;
	gS0			= gS0_400/400;
	GS1			= GS1_400/400;
	PI			= PI_400/400;
	PI0			= PI0_400/400;

	Hss=1;
	c_yss 	 	= 1-nu*(1-(varsigma-1)/varsigma);
	chi			=  (alpha*(1-omega)^sigmaL)*(varsigma-1)/varsigma/(((c_yss*Hss^alpha/(1-omega))-omega*Dc*c_yss*Hss^alpha/(1-omega))^(-sigmaC) * (Hss^alpha)^((1+sigmaL)/alpha-1)) ;
	
	% <----   initial state   ----> %
	tau0		= mu0^(theta2-1);
	SIG0		= E0/((1-mu0)*Y0);
	pb 			= 1000*theta2*0.109/(SIG0*(1-gS0)^(4*(2020-1985))*(1-deltapb)^(4*(2020-1985)));
	THETA0 		= pb/1000/theta2*delthet0*SIG0;
	d0 			= exp(-gamma*(M0-M_1750));
    pi_bar0     = PI0; 
	pi0 		= pi_bar0;
	mc0 		= kappa/varsigma*((1+pi0)*(pi0-pi_bar0) - (1-varsigma)/kappa - (1-nu)*beta*(1+gZ0)*(pi0-pi_bar0)*(1+pi0));
	gL0 		= (LT/L0)^lg-1;	
	c_y0 	 	= 1-THETA0*mu0^theta2-(1-nu)*kappa/2*(pi0-pi_bar0)^2-nu*(1-mc0);
	yy0 		= d0*h0^alpha;
	c0			= yy0*c_y0;
	lb0  		= ((c0/(1-omega)-omega*(Dc*d0*c0))/(1-omega))^-sigmaC;
	w0			= alpha*yy0/h0*(mc0-THETA0*mu0^theta2-THETA0*tau0*(1-mu0));
	Z0			= Y0/(L0*yy0);
    rr0 		= rr0_400/400; % set different set 
    
	% <------ terminal state -----> %
	SIG		= 0;
	gSIG	= 0;
	THETA1	= 0;
	delthet	= 0;
	tau		= tau0;
	pi		= PI;
	L		= LT;
	gZ      = 0;
	gL      = 0;
	mc 		= kappa/varsigma*((1+pi)*(pi-PI) -(1-varsigma)/kappa - (1-nu)*beta*(1+gZ)*(pi-PI)*(1+pi));
	[Z,M]   = get_Z(Z0,gZ0,GZ1,L0,LT,lg,gS0,SIG0,GS1,deltapb,tau0,theta2,pb,chi,omega,Dc,Et,M_1750,nu,mc,sigmaC,xi,M0,tf,gamma,alpha,sigmaL);
	d		= exp(-gamma*(M-M_1750));
%	Z 		= get_Z(Z0,gZ0,GZ1);
	x 		= (1-tf*THETA1*(tau0+Et*e_tau)^(theta2/(theta2-1))-nu*(1-mc));
	junk  	= 0;
	dy 		= log((1+gZ)*(1+gL));
	de		= 0;
	mu		= (tau)^(1/(theta2-1));
	h 		= (d^(1-sigmaC)*(x*(1-omega*d*Dc)/(1-omega))^-sigmaC*alpha*(1-omega)^sigmaL/(chi)*(mc-THETA1*mu^theta2 - THETA1*tau*(1-mu)))^(1/(1+sigmaL-alpha*(1-sigmaC)));
	y 		= d*h^alpha;
	w 		= alpha*y/h*(mc-THETA1*mu^theta2 - THETA1*tau*(1-mu)); 
	c		= y*x;
	D 		= d*Dc*c;
	lb  	= (c/(1-omega)-omega*D/(1-omega))^-sigmaC;
	r		= (1+pi)/(((1-omega) + omega*(D^-sigmaC)/lb)*beta) -1;
	s_p 	= 1; s_b = 1; s_e = 1;
	E		= 0;
	s_r  	= 1 ;
	pi_core	= PI;
	y_n		= y;
	x_n		= x;
	rr_n	= (1+r)/(1+pi)-1;
	pi_bar 	= PI;
	tau0_USD = 0*tau0*tf*THETA0/SIG0*1000;
	lny = 100*ln(y/y0) 	;
    ygap 	= 1/((((1+pi)/((1+pi_bar)))^((1+phi_pi)) * (y/((1-MP_y_n)*y0+MP_y_n*y_n))^(phi_y)));
	
	
	@#if SIMULATIONS == 1
		rotemberg_to_gdp = 100*((1-nu)*kappa/2*(pi-pi_bar)^2+nu*(1-s_p*mc));
		abatement_to_gdp = 100*THETA1*mu^theta2;
		tau_USD			= tau*THETA1/(SIG+1e-8)*1000;
		SDF 			= beta;
		C 				= c*Z*L;
		IS				= ln(c/y0)*100;
		
		mc_no_damage 	= (chi/(alpha*(1-omega)^sigmaL))*((x*y/(1-omega))-omega*Dc*c/(1-omega))^(sigmaC) * y^((1+sigmaL)/alpha-1) ;
		mc_hat   		= 100*(mc-mc_no_damage)/mc_no_damage; 
		mc_hat_s 		= 100*((chi/(alpha*(1-omega)^sigmaL))*((x*y/(1-omega))-omega*Dc*c/(1-omega))^(sigmaC) * y^((1+sigmaL)/alpha-1) - mc_no_damage)/mc_no_damage;
		mc_hat_m 		= 100*((chi/(alpha*(1-omega)^sigmaL))*((x*y/(1-omega))-omega*D/(1-omega))^(sigmaC) * y^((1+sigmaL)/alpha-1) / (exp(-gamma*(M-M_1750)))^((1+sigmaL)/alpha) - (chi/(alpha*(1-omega)^sigmaL))*((x*y/(1-omega))-omega*Dc*c/(1-omega))^(sigmaC) * y^((1+sigmaL)/alpha-1) )/mc_no_damage ;
		mc_hat_g 		= 0;
		lb_n			= lb;
		rr_n			= (1+r)/(1+pi)-1;
		c_n				= c;
		y_n				= y;
		rr_real			= rr_n;
		pi_hat_x 		= 0; 
		pi_hat_g = 0;
		pi_hat_m 		= (varsigma-1)/kappa*mc_hat_m/(1-(1-nu)*beta*(1+gZ));
		pi_hat_s 		= 100*(pi-pi_bar) - pi_hat_x - pi_hat_m - pi_hat_g;
		pi_tot   		=  pi_hat_x + pi_hat_s + pi_hat_m + pi_hat_g;
		ln_gap  		= omega*log((D^-sigmaC)/lb_n) - omega*Dc^-sigmaC;
		lny_n 			= 100*ln(y_n/y0) 	;
		pi100			= 100*pi;
		rr100 			= 100*r;
		rr_real100 		= 100*rr_real;
		rrn100 			= 100*rr_n;

		w = (chi*((y/exp(-gamma*(M-M_1750)))^(1/alpha)/(1-omega))^sigmaL*((x*y/(1-omega))-omega*(D)/(1-omega))^(sigmaC));

		
		WELF   = (Z)^(1-sigmaC)*(L)*(1/((1-sigmaC))*    c    ^(1-sigmaC)-chi/(1+sigmaL)*(              h/(1-omega))^(1+sigmaL))/(1-beta);
		WELF_n = (Z)^(1-sigmaC)*(L)*(1/((1-sigmaC))*(x_n*y_n)^(1-sigmaC)-chi/(1+sigmaL)*((y_n/d)^(1/alpha)/(1-omega))^(1+sigmaL))/(1-beta);
		SCC    = (1000*3/11*gamma/alpha*chi*Z*L*c^sigmaC*((y/d)^(1/alpha)/(1-omega))^(1+sigmaL))/(1/beta-1);
		SCC_n  = (1000*3/11*gamma/alpha*chi*Z*L*c^sigmaC*((y/d)^(1/alpha)/(1-omega))^(1+sigmaL))/(1/beta-1);
		
		@#endif
	@#if SIMULATIONS == 2
	@#endif
end;

resid;
steady;


gen_histval_func(M_);

options_.initial_guess_path = 'SSV_sims0';
options_.expectation_window = 100;
options_.forward_path 		= 3000-options_.expectation_window;
options_.ep.Tdrop 			= 1; 
