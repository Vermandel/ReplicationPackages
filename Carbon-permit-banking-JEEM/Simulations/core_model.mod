//*************************************************************
//
// A General Equilibrium Approach to Carbon Permit Banking
//
// L. Dubois, J.-G. Sahuc, G. Vermandel  - August 2024
//
// E-DSGE Model
//
//*************************************************************


@#define variable_phi1 = 1

var y c n e mu w pe thet b mc lb k e_a e_pe  e_e lambda q i dy_obs pe_obs  de_obs mu_cost;

varexo eta_a eta_e eta_pe;

parameters  alpha sigma phi_1 phi_2 beta varphi psi nu cap chi delta E Y A N chi_I h L rho_a rho_pe rho_e MU gamma rss p ;



%-------------------------------
% 1. Calibration
%-------------------------------       

alpha 	= 1/3;        						% Technology parameter
varphi 	= 1;         						% Inverse of the frisch elasticity of labor supply
sigma 	= 2;          						% Intertemporal elasticity of consumption
rss     = (1+1/100)^(1/12);					% Steady-state risk-free rate
nu 	= 0.2;          					% Elasticity of e wrt y (heutel2012: 1-nu = .696)
delta	= 0.01/3;						% Capital depreciation rate
phi_1	= 0.1;							% Steady-state backstop technology
phi_2	= 2.6;							% Parameter of the abatement cost function
E 	= 1; 							% 1Gt quarterly
Y	= 14.5/4;
N	= 38*4;
MU    	= .3;
chi_I  	= 4;
h	= 0.7;
rho_a 	= 0.408227792033937;					% Autocorrelation parameter - productivity shock
rho_pe 	= 0.881793185217435; 					% Autocorrelation parameter - carbon price shock
rho_e	= 0.662108179477432;					% Autocorrelation parameter - carbon emission shock
L       = 400;							% EU population
gamma	= 1.0016;						% GDP growth - monthly


@#if  variable_phi1 
	var  		deltaphi;
	parameters 	phi1min;
	phi1min     = phi_1*(1-1/1200)^1525;
@#else
	parameters deltaphi phi1min;
	deltaphi	= 1;
	phi1min     = 0.1;
@#endif


load('DSV_mode.mat')

for iParam = 4 : size(e_obj.thet_ids,2)
    set_param_value(deblank(e_obj.theta.names(iParam,:)), theta(iParam))
end


%-------------------------------
% 2. Steady state
%-------------------------------  

steady_state_model;
    e_a = 1; e_pe = 1; e_e = 1;
	beta    = (gamma^(sigma-1))/rss;
	betahat = beta*gamma^-sigma;
	betehat = beta*gamma^(1-sigma);
	p 	= 1;
	q 	= 1;
	% pre-state no shocks (i.e. perturbations state)
	cap 	= E;
	psi 	= E/((1-MU)*Y^(1-nu));
	pe0 	= 1000*phi_1*phi_2*MU^(phi_2-1)*Y^nu/psi;
	mc0	= (p-phi_1*MU^(phi_2)-pe0/1000*(1-nu)*E/Y);
	k0	= alpha*Y*mc0/(1/betahat-(1-delta));
	c0 	= Y - phi_1*MU^(phi_2)*Y - pe0/1000*E - (gamma - (1-delta))*k0;
	A 	= Y/((L*N)^(1-alpha)*k0^alpha);
	chi	= (c0-h/gamma*c0)^-sigma*(1-alpha)*Y/N*mc0/(N^varphi);
    	thet 	= e_e*cap;
	b 	= 0;
	deltaphi = phi1min/phi_1;
%    deltaphi = 1;
	[mc,k,c,n,y,pe,mu] = policy_anticipations_ss(e_pe,e_a,deltaphi*phi_1,phi_2,psi,nu,thet,p,A,beta,alpha,h,L,varphi,sigma,delta,chi,mc0,k0,c0,N,Y,pe0,MU,gamma);	
    	e = psi*(1-mu)*y^(1-nu);                                                             
	w 	= (1-alpha)*y/n*mc  ;
	i 	= (gamma - (1-delta))*k;
    lb 	= (c-h/gamma*c)^-sigma;
	lambda= pe*(1-betehat);
 	dy_obs = log(gamma)*100; pe_obs = 0; de_obs = 0;
    	mu_cost = e_pe*deltaphi*phi_1*(mu^phi_2)*100 ;
	
end;

%-------------------------------
% 3. Model
%-------------------------------  

model;%(bytecode);

    %***********
    % Households
    %***********

    lb*w/p  = chi*(n^varphi);                 									                                        % FOC-n
    lb      = (c-h/gamma*c(-1))^-sigma;										                                            % FOC-c
    q       = p + chi_I*(i/i(-1)*gamma-gamma) - beta*gamma^-sigma*lb(+1)/lb*0.5*chi_I*((i(+1)/i*gamma)^2-gamma^2);   	% FOC-i
    beta*gamma^-sigma*lb(+1)*(alpha*y(+1)/k*mc(+1)+q(+1)*(1-delta)) = q*lb ;					                        % FOC-k
    i       = gamma*k - (1-delta)*k(-1);										                                        % capital accumulation

    %***********
    % Firms
    %***********

    y = e_a*A*(L*n)^(1-alpha)*k(-1)^alpha;                 						% Production function
    e = psi*(1-mu)*y^(1-nu);                               						% Emission process
    b = b(-1) + thet  - e;                                 						% Bank of allowances
    (1-alpha)*y/n*mc = w ;   										            % FOC-n
    mc = p-e_pe*deltaphi*phi_1*(mu^phi_2)-pe/1000*(1-nu)*e/y;					% FOC-y
    e_pe*deltaphi*phi_1*phi_2*(mu)^(phi_2-1)*y = pe/1000*psi*y^(1-nu);          % FOC-mu                   
    pe = beta*gamma^(1-sigma)*lb(+1)/lb*pe(+1) + lambda;                        % FOC-b                                       
    mu_cost = e_pe*deltaphi*phi_1*(mu^phi_2)*100 ;     							% Abatement cost/y

    %**********************
    % Regulatory authority
    %**********************

    @#ifndef msr
    thet = e_e*cap;
    @#else 
    thet = max(e_e*cap - (b>(833*12/1529))*(b<(1096*12/1529))*(b-(833*12/1529))/12 -(b>(1096*12/1529))*(0.24/12)*b,0.029758011772404) ;
    @#endif


    %*************
    % Equilibrium
    %*************
                                                 
    y = c + e_pe*deltaphi*phi_1*mu^phi_2*y + pe/1000*thet + (i+i(-1)/gamma*(chi_I/2)*(i/i(-1)*gamma-gamma)^2);  

    %********
    % Shocks
    %********

    e_a  = 1-rho_a + rho_a*e_a(-1) + eta_a;        	% Productivity shock
    e_pe = 1-rho_pe + rho_pe*e_pe(-1) + eta_pe;    	% Carbon price shock
    e_e  = 1 + eta_e;                              	% Carbon emission shock

    %*************
    % Observables
    %*************

    dy_obs   = 100*log(gamma*y/y(-1));			% Output growth
    pe_obs  =  (pe - steady_state(pe));			% Carbon price
    de_obs   = 100*log(e/e(-1));			    % Emissions growth

    [ mcp = 'lambda > 0' ]
    b = 0;
	
	@#if  variable_phi1 
		deltaphi = max(phi1min/phi_1,deltaphi(-1)*(1-1/1200));
	@#endif 

end;  

steady;

