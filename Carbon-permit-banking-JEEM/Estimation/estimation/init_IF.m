
options_bak = options_;

% housekeeping of persistent variables 
clear IF_EP dsge_llk;

% solve the model up to second order
% with minimal setting
options_.irf=0;
options_.noprint=1;
options_.nocorr=1;
options_.ar=0;
options_.nograph=1;
options_.nodisplay=1;

% dynare version
options_.dynare_versx = split(dynare_version,'.');
options_.dynare_vers    = nan(1,2);
options_.dynare_vers(1) = str2double(options_.dynare_versx{1});
options_.dynare_vers(2) = str2double(options_.dynare_versx{2});

if options_.dynare_vers(1) ~= 5
	error('Dynare version must be 5.x')
end


if ~isfield(options_,'penalized_function')
	options_.penalized_function = 1e7;
end


% check names across dynare versions
if options_.dynare_vers(2) < 6 && options_.dynare_vers(1)==4
	M_.exo_names2 = M_.exo_names;
	M_.endo_names2 = M_.endo_names;
	M_.param_names2 = M_.param_names;
	M_.exo_names_long2 = M_.exo_names_long;
	M_.exo_names_tex2 = M_.exo_names_tex;
	M_.endo_names_long2 = M_.endo_names_long;
	M_.endo_names_tex2 = M_.endo_names_tex;
	M_.param_names_tex2 = M_.param_names_tex;
	M_.param_names_long2 = M_.param_names_long;
else
	M_.exo_names2 = char(M_.exo_names);
	M_.endo_names2 = char(M_.endo_names);
	M_.param_names2 = char(M_.param_names);
	M_.exo_names_long2 = char(M_.exo_names_long);
	M_.exo_names_tex2 = char(M_.exo_names_tex);
	M_.endo_names_long2 = char(M_.endo_names_long);
	M_.endo_names_tex2 = char(M_.endo_names_tex);
	M_.param_names_tex2 = char(M_.param_names_tex);
	M_.param_names_long2 = char(M_.param_names_long);
end


% impose pruning state space representation
if options_.order > 1
	options_.pruning=1;
end
% activate ksolver
if options_.order==3
	options_.k_order_solver=1;
end

% run stochastic stimulations
if options_.dynare_vers(2) < 6 && options_.dynare_vers(1) == 4 
	stoch_simul(M_.endo_names(1,:));
else
	var_list_ = {};
	[info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_);
end

%% First get infos
%dataset_.data
%dataset_.name
%bayestopt_.name
e_obj.id.shocks		= [];
e_obj.id.params		= [];
e_obj.idu.params	= [];
e_obj.idu.shocks	= [];
for i1 = 1:size(bayestopt_.name,1)
	idx = strmatch(bayestopt_.name{i1},M_.param_names,'exact');
	if ~isempty(idx)
		e_obj.id.params = [e_obj.id.params; idx];
		e_obj.idu.params = [e_obj.idu.params; i1];
	else
		idx = strmatch(bayestopt_.name{i1},M_.exo_names,'exact');
		if ~isempty(idx)
			e_obj.id.shocks = [e_obj.id.shocks; idx];
			e_obj.idu.shocks = [e_obj.idu.shocks; i1];
		end
	end
end

%% ESTIMATED PARAMETERS
% START WITH SHOCKS STD
e_obj.theta.lb = [estim_params_.var_exo(:,3);estim_params_.param_vals(:,3)];
e_obj.theta.ub = [estim_params_.var_exo(:,4);estim_params_.param_vals(:,4)];
e_obj.theta.st = [estim_params_.var_exo(:,2);estim_params_.param_vals(:,2)];
e_obj.theta.id = [estim_params_.var_exo(:,1);estim_params_.param_vals(:,1)];
e_obj.theta.is_sd = 1:size(estim_params_.var_exo,1);
e_obj.theta.is_pm = (e_obj.theta.is_sd(end)+1):(e_obj.theta.is_sd(end)+size(estim_params_.param_vals,1));
e_obj.theta.names = char(M_.exo_names2(estim_params_.var_exo(:,1),:),M_.param_names2(estim_params_.param_vals(:,1),:));
theta0				= e_obj.theta.st;
% get dynare extremums
e_obj.theta.lb = max([bayestopt_.p3 e_obj.theta.lb],[],2);
e_obj.theta.ub = min([bayestopt_.p4 e_obj.theta.ub],[],2);

% check empty priors starting value
for i1=1:length(theta0)
	% then replace by prior mean
	if isnan(theta0(i1))
		theta0(i1) = bayestopt_.p1(i1);
	end
end

% set theta into bayesopt ids
e_obj.thet_ids = [];
for i1 = 1:length(theta0)
	e_obj.thet_ids = [e_obj.thet_ids strmatch(deblank(e_obj.theta.names(i1,:)),bayestopt_.name,'exact')];
end

% look for observation equations
e_obj.id.obs =[];
for i1 = 1:size(dataset_.name,1)
	idx = strmatch(dataset_.name{i1},M_.endo_names,'exact');
	if ~isempty(idx)
		e_obj.id.obs = [e_obj.id.obs;idx];
	end
end

% set estimated std as first elements of theta 
e_obj.theta.M_sd = estim_params_.var_exo(:,1);


% create the selection matrix
e_obj.Q=zeros(length(e_obj.theta.is_sd),M_.endo_nbr);
for i = 1:length(e_obj.id.obs)
	e_obj.Q(i,e_obj.id.obs(i))=1;
end

Tstart = 1;%+options_.first_obs;%+options_.presample;
obs = dataset_.data(Tstart:end,:);
if options_.prefilter
	obs=obs-mean(obs);
end

%% if reload
if isfield(options_,'reload_last') > 0 && options_.reload_last 
    load([ M_.fname '_mle_estimates_temp.mat'],'params')
    load('ve_names.mat')
    for i1=1:length(ve_names)
       idx = find(strcmp(bayestopt_.name,ve_names{i1}));
       if ~isempty(idx)
           theta0(idx) = params(i1);
       end
    end
end

%%%
disp('STARTING HIGHER ORDER ESTIMATION')
if isempty(options_.mode_file)
	
	% save variable names
	ve_names = bayestopt_.name;
	save ve_names ve_names;
	% initialize priors persistent variables
	lnprior = priordens(theta0(e_obj.thet_ids),bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4,1);        
	
	[llk] = dsge_llk(theta0,obs,e_obj,oo_,M_,options_,bayestopt_);
	disp(['Initial value of Log-likelihood is ' num2str(-llk) '.']);

	if options_.mode_compute ==1	 %FMINCON

		% Set default optimization options for fmincon.
		%optim_options = optimset('Algorithm','interior-point','display','iter', 'LargeScale','off', 'MaxFunEvals',100000, 'TolFun',1e-8, 'TolX',1e-8,'TolFun',1e-8,'DiffMinChange',1e-8,'UseParallel',false,'OutputFcn', @(x,optimValues,state)fminOut(x,optimValues,state));
		optim_options = optimset('Algorithm','active-set','display','iter', 'LargeScale','off', 'MaxFunEvals',100000, 'TolFun',1e-8, 'TolX',1e-8,'TolFun',1e-8,'DiffMinChange',1e-8,'UseParallel',true,'OutputFcn', @(x,optimValues,state)fminOut(x,optimValues,state,e_obj));
		
		%% IF USE PARALLEL COMP
		%% PERSISTENT VARIABLES ARE TROUBLESOME FOR PRIORDENS
		%% NEED TO REWRITE IT LATER!
		e_obj.theta.lb(find(e_obj.theta.lb==0))=.0000001;
		e_obj.theta.lb(isinf(e_obj.theta.lb==0))=-10^6;
		e_obj.theta.ub(isinf(e_obj.theta.ub==0))=10^6;
		[theta,fval,exitflag,output,lamdba,grad,hessian_mat] = fmincon(@(theta) dsge_llk(theta,obs,e_obj,oo_,M_,options_,bayestopt_),theta0,[],[],[],[],e_obj.theta.lb,e_obj.theta.ub,[],optim_options);
	elseif options_.mode_compute == 9 % CMAES
		%set CMAES options
		% Set defaults
		H0 = (e_obj.theta.ub-e_obj.theta.lb)*0.2;
		H0(~isfinite(H0)) = 0.01;
		while max(H0)/min(H0)>1e6 %make sure initial search volume (SIGMA) is not badly conditioned
			H0(H0==max(H0))=0.9*H0(H0==max(H0));
		end

		cmaesOptions.SaveVariables	= 'no';
		cmaesOptions.DispFinal	= 'on';
		cmaesOptions.WarnOnEqualFunctionValues	= 'no';
		cmaesOptions.DispModulo	= '10';
		cmaesOptions.LogModulo	= '0';
		cmaesOptions.LogTime	= '0';
		cmaesOptions.TolFun		= 1e-4;
		cmaesOptions.TolX		= 1e-6;
		cmaesOptions.Resume		= 0;
		cmaesOptions.LBounds	= e_obj.theta.lb;
		cmaesOptions.UBounds	= e_obj.theta.ub;
		cmaesOptions.MaxFunEvals= 10000;
		%H0 = ones(size(H0))*.01;
		%H0=mean([ones(size(H0))*.01,H0],2);
		H0 = 0.01+(H0-0.01)/15;
		[theta, loss_end, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = cmaes_estobin('dsge_llk',theta0,H0,cmaesOptions,obs,e_obj,oo_,M_,options_,bayestopt_);
	elseif options_.mode_compute == 10

    simplexOptions = options_.simplex;
    if options_.silent_optimizer
        simplexOptions.verbosity = 0;
    end
    [opt_par_values,fval,exitflag] = simplex_optimization_routine('dsge_llk',theta0,simplexOptions,cellstr(e_obj.theta.names),obs,e_obj,oo_,M_,options_,bayestopt_);

	elseif options_.mode_compute == 6
		hh=[];
		[opt_par_values, hessian_mat, Scale, fval] = gmhmaxlik('dsge_llk', theta0, ...
                                                      hh, options_.mh_jscale, [e_obj.theta.lb e_obj.theta.ub], bayestopt_.p2, options_.gmhmaxlik, options_.optim_opt,...
													  obs,e_obj,oo_,M_,options_,bayestopt_);
	
	elseif options_.mode_compute==7
		e_obj.theta.lb(find(e_obj.theta.lb==0))		= 1e-6;
		e_obj.theta.ub(isinf(e_obj.theta.ub==0))	= 10^6;
		tolerance = 1e-05;
	  options = optimset('Display','on','TolFun',tolerance,'TolX',tolerance,'TolFun',tolerance,...
		'DiffMinChange',tolerance,'Algorithm','Interior-Point','OutputFcn', @(x,optimValues,state)fminOut(x,optimValues,state,e_obj), ...
		'MaxFunEvals',25000,'MaxIter',25000);
	  [theta,fval]=fminsearchbnd(@(theta) dsge_llk(theta,obs,e_obj,oo_,M_,options_,bayestopt_),theta0,e_obj.theta.lb,e_obj.theta.ub,options);
    elseif options_.mode_compute ==0
        % do nothing
        theta=theta0;
    else
		error('optimizer unknown.')
	  
	end

	% once estimated compute shocks
	[llk,ee,res,yy] = dsge_llk(theta,obs,e_obj,oo_,M_,options_,bayestopt_);
	% save result
	eval(['save ' M_.fname '_mode.mat e_obj theta llk ee res']);
else
	eval(['load ' options_.mode_file ' theta;'])
	
end

%% COMPUTE COVARIANCE
if (options_.mh_replic > 0 && options_.load_mh_file == 0)
	
	% first reload last optimization
	if ~isempty(options_.mode_file)
		eval(['load ' options_.mode_file ' theta;']);
		disp('Mode from previous estimation successfully loaded.');
	else
		eval(['load ' M_.fname '_mode theta;']);
	end
	% compute initial llk
	[llk,ee,res,yy] = dsge_llk(theta,obs,e_obj,oo_,M_,options_,bayestopt_);
	
	
	
	%% FIRST GET COVARIANCE ACROSS ESTIMATED PARAMETERS
	disp('Compute covariances across parameters...')
	
	try
		%% DYNARE FASHION: COMPUTE HESSIAN TO GET COVARIANCES
	
%		hhh = hessian('dsge_llk',theta, options_.gstep,obs,e_obj,oo_,M_,options_,bayestopt_);
		hhh = hessian2('dsge_llk',theta, options_.gstep,obs,e_obj,oo_,M_,options_,bayestopt_);

		hessian_mat = reshape(hhh,[length(theta) length(theta)]);
		check = chol(inv(hessian_mat));
		disp('Hessian computed!')
	catch
		disp('hessian not finite...')
		[idx,idx2]=find(hessian_mat==min(min(hessian_mat)));
		disp(['Maximum hessian is for: ' e_obj.theta.names(idx(1),:) ' and ' e_obj.theta.names(idx(2),:)])
		% use priors to generate the initial varcov matrix of parameters
%		hh = diag(1./(bayestopt_.p2.^2));
%		hsd = sqrt(diag(hh));
%		invhess = inv(hh./(hsd*hsd'))./(hsd*hsd');
%		hessian_priors  = inv(invhess);
%		save hessian_priors;
        load([M_.fname '_mle_estimates_temp.mat'])
        [ hessian_reg stdh_reg hessian_fmin stdh_fmin ] = compute_hessian(xstory',fstory,30);
        hessian_mat = diag(diag(hessian_fmin));
        disp('Hessian computed!')
        disp('')
 

	end

    %% save hessian
    eval(['save ' M_.fname '_hessian_mat hessian_mat']);
		

end

% run MH using covariances computed in the previous conditional statement
if options_.mh_replic > 0 
	% first reload last optimization
	if ~isempty(options_.mode_file)
		eval(['load ' options_.mode_file ' theta;']);
	else
		eval(['load ' M_.fname '_mode theta;']);
	end
	eval(['load ' M_.fname '_hessian_mat']);

	disp('Starting Metropolis-Hasting Iterations')
	
	% Start iterations
	%try
	% try parallel setting
	mhalgo_init('dsge_llk',theta,hessian_mat,options_.mh_nblck,options_.mh_replic,options_.load_mh_file,options_.mh_jscale,obs,e_obj,oo_,M_,options_,bayestopt_);
	%catch
	%mhalgo_init_simple('dsge_llk',theta,hessian_mat,options_.mh_nblcks,options_.mh_replic,options_.load_mh_file,options_.mh_jscale,obs,e_obj,oo_,M_,options_,bayestopt_);
	%end
	check_metropolis;
	
end

% try to load previous estimation
if options_.load_mh_file
	check_metropolis;
	theta = mean_metropolis;
elseif sum(bayestopt_.pshape)==0
    % compute llk uncertainty
        eval(['load ' M_.fname '_mle_estimates_temp xstory fstory'])
        [ hessian_reg stdh_reg hessian_fmin stdh_fmin ] = compute_hessian(xstory',fstory,200);
        hessian_mat = diag(inv(hessian_fmin));
        disp('Hessian computed!')
        disp('')
 
end

disp('ESTIMATED PARAMETERS:')
for i1=1:size(e_obj.theta.names,1)
	fprintf('%s %.8f \n',e_obj.theta.names(i1,:),theta(i1))
end
[llk,ee,res,yy] = dsge_llk(theta,obs,e_obj,oo_,M_,options_,bayestopt_);
disp(['VALUE OF THE LOG-POSTERIOR: ' num2str(-llk)])

%% update the model's
M_.params(e_obj.theta.id(e_obj.theta.is_pm)) = theta(e_obj.theta.is_pm);
if length(e_obj.idu.shocks) == M_.exo_nbr
    M_.Sigma_e = zeros(M_.exo_nbr);
    M_.Sigma_e(1:length(e_obj.idu.shocks)+1:end)=theta(e_obj.theta.is_sd).^2;
    M_.Sigma_e(e_obj.theta.id(e_obj.theta.is_sd),e_obj.theta.id(e_obj.theta.is_sd)) = M_.Sigma_e;
else
   for i1 = 1:length(e_obj.idu.shocks)
        M_.Sigma_e(e_obj.idu.shocks(i1),e_obj.idu.shocks(i1)) = theta(e_obj.theta.is_sd(i1)).^2;
   end
end



% if perturbations methods
% update policy rule and steady state
if isfield(options_.ep,'estim') 
	if ~options_.ep.estim
		% update policy rule
		[oo_.dr,~,M_,~,oo_] = resol(0,M_,options_,oo_);
	end
end

if ~exist('pw') % (non-)linear solution from Dynare
	
	[llk,ee,res,yy_] = dsge_llk(theta,obs,e_obj,oo_,M_,options_,bayestopt_);
	if res > 0.00001
		warning('Inversion did not worked at mean, get mode.')
		load([M_.fname '_mode'],'theta')
		[llk,ee,res,yy_] = dsge_llk(theta,obs,e_obj,oo_,M_,options_,bayestopt_);	
	end
	
	if 	~options_.ep.estim
		yy_=simult_(M_,options_,oo_.dr.ys,oo_.dr,ee,options_.order)
		xx=e_obj.Q*yy_;
	else
		xx=e_obj.Q*yy_;
	end
	
	
	% try to get the x-axis (dates)
	try
		load(options_.datafile,'Tq');
		Tnum=Tq;
        Tnum=Tnum(options_.first_obs:options_.nobs)
	catch
    end
    if ~exist('Tnum')
        	Tnum = 1:length(obs(:,1));
    end;
    
    ls = size(obs,1);

    if ls ~= length(Tnum)
        if length(Tnum) > ls
            Tnum = Tnum(1:ls);
        else
            freq = mean(diff(Tnum));
            gap = ls-length(Tnum) ;
            Tnum = [-flip(1:gap)*freq+Tnum(1);Tnum]
        end
    end

	% model vs data
	figure;
	for i1=1:length(e_obj.id.obs)
		subplot(ceil(sqrt(length(e_obj.id.obs))),ceil(sqrt(length(e_obj.id.obs))),i1)
		idd = strmatch(deblank(dataset_.name{i1}),M_.endo_names2,'exact');
		plot([Tnum(1) Tnum(end)], [oo_.dr.ys(idd) oo_.dr.ys(idd)],':','Color',[0.1 0.1 0.1],'LineWidth',1.2,'HandleVisibility','off')
		hold on;
		plot(Tnum,xx(i1,1+(1:length(Tnum))),'Color',[0 102 204]/255,'LineWidth',1.2)
		plot(Tnum,obs(:,i1),'or', 'LineStyle', 'none' )
		hold off
		if isfield(M_,'endo_names_tex') && ~isempty(M_.endo_names_tex2(idd,:)) && ~isempty(M_.endo_names_long2(idd,:))
				title([  M_.endo_names_long2(idd,:) ' ($' deblank(M_.endo_names_tex2(idd,:)) '$)'],'interpreter','latex')	
		else
			title(['Observable ' dataset_.name{i1}])
		end
		xlim([Tnum(1) Tnum(end)])
	end
	legend('Model','Data')
	
	% plot Smoothed shocks
	figure('Name','Smoothed shocks');
	for i1=1:length(e_obj.id.shocks)
		subplot(ceil(sqrt(length(e_obj.id.shocks))),ceil(sqrt(length(e_obj.id.shocks))),i1)
		idd = strmatch(deblank(M_.exo_names(e_obj.id.shocks(i1),:)),M_.exo_names2,'exact');
		plot([Tnum(1) Tnum(end)], [0 0],':','Color',[0.1 0.1 0.1],'LineWidth',1.2,'HandleVisibility','off')
		hold on;
		plot(Tnum,ee(1:length(Tnum),e_obj.id.shocks(i1)),'Color',[0 102 204]/255)
		hold off;
		if isfield(M_,'exo_names_tex') && ~isempty(M_.exo_names_tex2(idd,:)) && ~isempty(M_.exo_names_long2(e_obj.theta.id(idd),:))
			title([  M_.exo_names_long2(idd,:) ' ($' deblank(M_.exo_names_tex2(idd,:)) '_t$)'],'interpreter','latex')	
		else
			title([M_.exo_names(idd,:) ])					
		end
		xlim([Tnum(1) Tnum(end)])
	end

else
    smooth_pw;
end
%options_ = options_bak;


