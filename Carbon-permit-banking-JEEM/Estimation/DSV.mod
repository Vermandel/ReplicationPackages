//*************************************************************
//
// A General Equilibrium Approach to Carbon Permit Banking
//
// L. Dubois, J.-G. Sahuc, G. Vermandel  - August 2024
//
// E-DSGE Model - Estimation
//
//*************************************************************

% Define constraints to bind
@#define c1_bind = 0
close all;

addpath('estimation');


var y c n e mu w pe thet b mc lb k e_a e_pe  e_e lambda q i  dy_obs pe_obs de_obs dc di;

varexo eta_a eta_e  eta_pe ;

parameters  alpha sigma phi_1 phi_2 beta varphi psi nu cap chi delta gamma rss E Y A N chi_I h L p rho_a rho_pe  rho_e MU;

%-------------------------------
% 1. Calibration
%-------------------------------  

@#if c1_bind == 0

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
MU    	= .25;
chi_I  	= 4;
h	= 0;
rho_a 	= 0.408227792033937;					% Autocorrelation parameter - productivity shock
rho_pe 	= 0.881793185217435; 					% Autocorrelation parameter - carbon price shock
rho_e	= 0.662108179477432;					% Autocorrelation parameter - carbon emission shock
L       = 400;							% EU population
gamma	= 1.0016;	

%-------------------------------
% 2. Steady state
%-------------------------------  

steady_state_model;
	e_a = 1; e_pe = 1; e_e = 1;
	beta    = (gamma^(sigma-1))/rss;
	betahat = beta*gamma^-sigma;
	betehat = beta*gamma^(1-sigma);
	p 		= 1;
	q 		= 1;
	n 		= N;
	y 		= Y;
	mu 		= MU; 
	e 		= E; 
    cap 	= E;
    psi     = e/((1-mu)*y^(1-nu)); 
    thet    = cap;    
	pe 		= 1000*phi_1*phi_2*(mu)^(phi_2-1)*y^nu/psi;
	lambda	= pe*(1-betehat);
	mc 		= (p-e_pe*phi_1*(mu)^(phi_2)-pe/1000*(1-nu)*e/y);
	k		= alpha*y*mc/(1/betahat-(1-delta));
	A 		= y/(e_a*(L*n)^(1-alpha)*k^alpha);
	i 		= (gamma-1+delta)*k;
	c 		= y - phi_1*mu^(phi_2)*y - pe/1000*thet - i;
	w 		= (1-alpha)*y/n*mc  ;
    lb 		= (c-h/gamma*c)^-sigma;
	chi		= lb*w/(n^varphi);
	b 		= 0; 
    dy_obs  = log(gamma)*100; 
    pe_obs  = 0*pe; 
    de_obs  = 0;
	dc      = 100*log(gamma*c/c);
    di      = 100*log(gamma*i/i);
end;
@#endif

%--------------------------------------------------------------------------
% 3. Model
%--------------------------------------------------------------------------

model;

    %***********
    % Households
    %***********

    lb*w/p = chi*(n^varphi);                 									                                % FOC-n
    lb = (c-h/gamma*c(-1))^-sigma;										                                        % FOC-c
    q = p + chi_I*(i/i(-1)*gamma-gamma) - beta*gamma^-sigma*lb(+1)/lb*0.5*chi_I*((i(+1)/i*gamma)^2-gamma^2);   	% FOC-i
    beta*gamma^-sigma*lb(+1)*(alpha*y(+1)/k*mc(+1)+q(+1)*(1-delta)) = q*lb ;					                % FOC-k
    i = gamma*k - (1-delta)*k(-1);										                                        % capital accumulation



    %***********
    % Firms
    %***********

    y = e_a*A*(L*n)^(1-alpha)*k(-1)^alpha;                          % Production technology
    e = psi*(1-mu)*y^(1-nu);                                        % Emissions
    b = b(-1) + thet  - e;                                          % Bank of allowances
    (1-alpha)*y/n*mc = w ;                                          % FOC 1
    mc = p-e_pe*phi_1*(mu^phi_2)-pe/1000*(1-nu)*e/y;		        % FOC-y
    e_pe*phi_1*phi_2*(mu)^(phi_2-1)*y = pe/1000*psi*y^(1-nu);       % FOC-mu                   
    pe = beta*gamma^(1-sigma)*lb(+1)/lb*pe(+1) + lambda;            % FOC-b                                       


    %**********************
    % Regulatory authority
    %**********************

    thet = e_e*cap;

    %*************
    % Equilibrium
    %*************
          
    y = c + e_pe*phi_1*(mu^phi_2)*y + pe/1000*thet +  (i+i(-1)/gamma*(chi_I/2)*(i/i(-1)*gamma-gamma)^2);                                                    
     
    %********
    % Shocks
    %********

    e_a  = 1-rho_a + rho_a*e_a(-1) + eta_a;             % Productivity shock
    e_pe = 1-rho_pe + rho_pe*e_pe(-1) + eta_pe;         % Carbon price shock
    e_e  = 1-rho_e + rho_e*e_e(-1) + eta_e;             % Carbon emission shock

    %*************
    % Observables
    %*************

    dy_obs   = 100*log(gamma*y/y(-1));      % Output growth
    pe_obs  =  (pe - steady_state(pe));     % Carbon price
    de_obs   = 100*log(e/e(-1));            % Emissions growth
	dc      = 100*log(gamma*c/c(-1));
    di      = 100*log(gamma*i/i(-1));

	@#if c1_bind
		lambda = 0;
	@#else
		b + 0.000*b(-1) = 0;
	@#endif
	
	%(phi_1 = phi_1(-1)*(1-1/1200); % 1% decline by year
end;     


resid(1);

%-------------------------------
% 4. Estimation
%-------------------------------  

varobs dy_obs  de_obs pe_obs; 

estimated_params;
%%	PARAM NAME,		INITVAL,		LB,		UB,		PRIOR_SHAPE,		PRIOR_P1,		PRIOR_P2,		PRIOR_P3,		PRIOR_P4,		JSCALE
stderr eta_a,       0.00115163,         0,      50;%,     inv_gamma2_pdf,     .1,             200;
stderr eta_pe,      0.14701702,         0,      50;%,     inv_gamma2_pdf,     .1,             200;
stderr eta_e,       0.02997834,         0,      50;%,     inv_gamma2_pdf,     .1,             200;
rho_a,				0.97104583,         0,		0.99;%,		beta_pdf,			.5,             .2;
rho_pe,             0.92342373,         0,		0.99;%,		beta_pdf,           .5,             .2;
rho_e,				0.94678561,         0,		0.99;%,		beta_pdf,			.5,             .2;
nu,					0.20211098,         0,		1;%,      beta_pdf,  , ,   0.2,              0.05; %beta tight
gamma,				1.00155400,            1,		1.01;%,		gamma_pdf,          2.5,            .5;
MU,					0.30344853,            0.2,		0.5;%,		gamma_pdf,          2.5,            .5;
varphi,				1.10805133,            0,     5;%,		normal_pdf,          1,              .5;
sigma,				2.00,            0,     5;%,		normal_pdf,          2,              .25;
h,				    0.70382285,            0,     0.999;%,		beta_pdf,          0.7,              .1;
chi_I,				6.23488246,            0,     10;%,		normal_pdf,          4,              1.5;
end;

%try
estimation(datafile='Data_Estimation_062009_122019.xlsx',prefilter=0,xls_sheet=Feuil1,xls_range=A1:F128,  order=1, mode_compute=0, mh_replic=0, mh_jscale=0.9, mh_nblocks=8, plot_priors=0);
%catch
%end
%% solve the model
stoch_simul(order=1,irf=0,noprint,nocorr,nograph,nodisplay);

options_.reload_last		= 1;
options_.mode_compute 		= 7;
%options_.mode_file 			= 'DSV_mode.mat';

%% build the piece-wise solution method
bayestopt=bayestopt_;estim_params=estim_params_;dataset=dataset_;options=options_;
pw = build_pw_struct(oo_,M_,options_);
bayestopt_=bayestopt;estim_params_=estim_params;dataset_=dataset;

% impose constraints
pw.c1_bind  = 'lambda <= 0';

%gg=read_bind(pw.c1_relax,M_)
%newStr = regexprep(str,expression,replace)
% selection matrix per regime
% one row = one regime in pw.combinations
% one column is one corresponding shock/variable = 1 activated
pw.shock_per_regime = ones(size(pw.combinations,1),size(estim_params_.var_exo,1));
pw.obs_per_regime    = ones(size(pw.combinations,1),size(dataset_.data,2));
% now let's disable useless obs
for i1=1:size(pw.combinations,2)
	% shocks
	regname = ['c' num2str(i1) '_shock'];
	if isfield(pw,regname) % if user gave a restriction
		idx = strmatch(eval(['pw.' regname]),M_.exo_names,'exact');
		if ~isempty(idx)
			% search when the regime binds
			
			idx2 = find(pw.combinations(:,i1)==1);
			for i2 = 1:length(idx2)
				% disable shock
				pw.shock_per_regime(idx2(i2),idx) = 0;
			end
		else
			warning([ regname ' ignored'])
		end
	end
	% observables
	regname = ['c' num2str(i1) '_obs'];
	if isfield(pw,regname) % if user gave a restriction
		idx = strmatch(eval(['pw.' regname]),dataset_.name,'exact');
		if ~isempty(idx)
			% search when the regime binds
			idx2 = find(pw.combinations(:,i1)==1);
			for i2 = 1:length(idx2)
				% disable shock
				pw.obs_per_regime(idx2(i2),idx) = 0;
			end
		else
			warning([ regname ' ignored'])
		end
	end
end

%% ESTIMATION
e_obj.pw = pw;

options_.order = 1;
options_       = options;
init_IF;


%% SMOOTH VARIABLES
smooth_pw;
var_names = {'y','c','i','pe','b','mu','e','thet'};
ny = size(var_names,2);
figure;
for i1=1:ny
	subplot(ceil(sqrt(ny)),ceil(sqrt(ny)),i1)
	plot(y_pw(:,strcmp(var_names{i1},M_.endo_names)))
    hold on;
    plot(y_lin(:,strcmp(var_names{i1},M_.endo_names)))
    hold off;
    title(var_names{i1})
end



%pw = build_pw_struct(oo_,M_,options_);

%% IRFs
Tirf 		= 50;
sign_shock  = ones(M_.exo_nbr,1);
sign_shock(1) = 1;
sign_shock(3) = -1;
% 
IRF_mat     = nan(Tirf+1, M_.endo_nbr, M_.exo_nbr);
IRF_mat_lin = nan(Tirf+1, M_.endo_nbr, M_.exo_nbr);
for ix = 1: M_.exo_nbr
	% shock(t=0,t=T)
	zz = zeros(M_.exo_nbr,Tirf);
	% impulse D
	D = zeros(M_.exo_nbr,1);
	D(ix) = sign_shock(ix);
	% initial impulsion
	zz(:,2) = (chol(M_.Sigma_e)')*D;
	% simulate the non-linear model
	[IRF_mat(:,:,ix),IRF_mat_lin(:,:,ix)] =  simul_pw_mat(pw,zz);
end


var_names = {'y','pe','b','mu','lambda'};
ny = size(var_names,2);
for ix = 1: M_.exo_nbr
	figure('name',M_.exo_names{ix});
	for i1=1:ny
		subplot(ceil(sqrt(ny)),ceil(sqrt(ny)),i1)
		plot(IRF_mat(:,strcmp(var_names{i1},M_.endo_names),ix))
		hold on;
		plot(IRF_mat_lin(:,strcmp(var_names{i1},M_.endo_names),ix))
		hold off;
		title(var_names{i1})
	end
end
