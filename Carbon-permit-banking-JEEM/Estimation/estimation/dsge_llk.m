function [llk,ee,res,y_,lnprior,llks] = dsge_llk(theta,obs,e_obj,oo_,M_,options_,bayestopt_,varargin)
	%DSGE_INVERT_2O Summary of this function goes here
	%   Detailed explanation goes here
	
	
	persistent current_maxllk;
	persistent	best_exo;
	persistent	best_endogenous;

	oo = oo_;
	oo_.dr.ghx = [];
    
	% initialize
	llk = 0;
    if isempty(current_maxllk)
		% initialize
		current_maxllk 	= -inf;
		best_exo 		= zeros(size(obs,1),M_.exo_nbr);
		best_endogenous = 1;
		if exist('ee0') && size(ee0,1) == size(best_exo,1) && size(ee0,2) == M_.exo_nbr
			best_exo = ee0;
		end
	end
	
	% if theta is not a vector, vectorize it
    if size(theta,2) > 1
        theta = theta';
    end
	
	% solve the model/update matrices & steady state
    infos(1)=0;
	if sum(isnan(theta))==0

		% check parameters' interval
		theta2 = min(max(e_obj.theta.lb,theta),e_obj.theta.ub);
		% penalize objective function if bound not restricted
		llk = llk - 10^7*sum(abs(theta2-theta));

		theta=theta2;
		warning off;

		% try update policy functions
		%try

			% copy the structure
			M = M_;
			
			M.params(e_obj.theta.id(e_obj.theta.is_pm)) = theta(e_obj.theta.is_pm);
			M.Sigma_e = zeros(M.exo_nbr);
			M.Sigma_e(1:length(e_obj.idu.shocks)+1:end)=theta(e_obj.theta.is_sd).^2;
			M.Sigma_e(e_obj.theta.id(e_obj.theta.is_sd),e_obj.theta.id(e_obj.theta.is_sd)) = M.Sigma_e;

			if isfield(options_,'dynare_vers') && (options_.dynare_vers(2)<6 && options_.dynare_vers(1)==4)
				[oo.dr.ys, M.params,infoss] = feval([M.fname '_steadystate2'],oo.dr.ys, zeros(M_.exo_nbr,1), M.params);
			else
				[oo.dr.ys, M.params,infoss] = feval([M.fname '.steadystate'],oo.dr.ys, zeros(M_.exo_nbr,1), M.params);
			end
			
			if ~isfield(options_.ep,'estim')
				% perturbation methods
				[oo.dr,infos] = stochastic_solvers(oo.dr,0,M,options_,oo);
			else
				% deterministic simulations
				oo.steady_state=oo.dr.ys;
				% avoid error check for deterministic models
				oo.dr.ghx = 1;
			end
			
        %catch
       %     disp('dont compute PQ')
		%	infos(end+1)=1;
		%end
	
	else
		infos(end+1)=1;
	end	
	
	% check if solution looks ok
	if (sum(isnan(oo.dr.ys))+sum(isinf(oo.dr.ys))+ sum(abs(imag(oo.dr.ys)))) ~= 0 || sum(oo.dr.ys) > 10^8
		infos(end+1) = 1;
	end

	% if solution has good shape
	if sum(infos)==0
		
		% start likelihood estimation
		N 		= length(e_obj.id.obs);
		T		= size(obs,1)-options_.presample;
		T1		= 1+options_.presample;
		% initialize likelihood
		llk		= llk-N*T/2*log(2*pi)-T/2*log(det(M.Sigma_e));
		% use the inversion filter to extract the sequence of shocks
		if isfield(options_.ep,'estim')
			% extended path
			[llks,ee,res,y_,stop_dist] = IF_EP(obs,e_obj,oo,M,options_,best_exo,best_endogenous,varargin);		
		elseif options_.order==1 && ~isfield(e_obj,'pw')
			% linearized model
			[llks,ee,res,y_,stop_dist] = IF1(obs,e_obj,oo,M,options_,varargin);
		elseif options_.order==1 && isfield(e_obj,'pw')
			% piece-wise linear model
			[llks,ee,res,y_,stop_dist] = IF1_pw(obs,e_obj,oo,M,options_,varargin);		
		elseif options_.order==2
			% second order model
			[llks,ee,res,y_,stop_dist] = IF2(obs,e_obj,oo,M,options_,varargin);
		elseif options_.order==3
			% third order model
			[llks,ee,res,y_,stop_dist] = IF3(obs,e_obj,oo,M,options_,varargin);		
		end
        
		% check residuals
		if sum(res) < 0.0001% && stop_dist == 0
			res = 0;
        end
        
		% compute likelihood 
		llk = llk + sum(llks(T1:end)) - (res*options_.penalized_function+stop_dist*options_.penalized_function*10);
		
		% compute priors
        if sum(bayestopt_.pshape)>0
		    lnprior = priordens(theta(e_obj.thet_ids),bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);        
        else
            % likelihood estimation
            lnprior = 0;
        end

		%% check priors
		if lnprior == 0 && sum(bayestopt_.pshape)>0
			disp('Null prior')
		end
		
		%% add priors
		if isnan(lnprior) || isinf(lnprior) || isinf(llk) || isnan(llk)
			llk = -10^12;
		else
			llk = llk + lnprior;
		end
	
	else
		llk = -10^12;
		ee=nan;res=nan;y_=nan;
	end
	
	if current_maxllk < llk
		best_exo 		= ee;
		best_endogenous	= y_;
		current_maxllk  = llk;
	end
	
	% minimize the minus
	llk = - llk;
end

