function [endogenous_variables_paths,err,c] = solve_one_period(spfm_exo_simul,endogenous_variables_paths,M,oo,options_,ny,nrc,iyp,iyf,dynamic_model,iz,is,isf1,isf,nrs,isp,icf,dynamicmodel)
%SOLVE_ONE_PERIOD Summary of this function goes here
%   Detailed explanation goes here


    Tmax = size(spfm_exo_simul,1)-1;
    
	% Initialize while-loop index.
	t = 1;
   
	% Set period index.
	t = t+1;
				



	%% solve each period
	iter = 1;
	error_growth = 0;
	info.iterations.error = nan(options_.simul.maxit,1);
	for iter = 1:options_.simul.maxit
				
		% initial c
		c = zeros(ny*(Tmax),nrc);
				
		% this is the current date 
		it_ = t;
		it_init = t;
				
		%% FIRST PERIOD
		z = [ endogenous_variables_paths(iyp,it_-1) ; endogenous_variables_paths(:,it_) ; endogenous_variables_paths(iyf,it_+1) ];
		%[d1,jacobian] = feval(dynamic_model,z,spfm_exo_simul, M.params, oo.steady_state, 2+it_-t);
        [d1, jacobian] = dynamicmodel(z, spfm_exo_simul, M.params, oo.steady_state, 2+it_-t);
	

 
		% [Jac_y(t-1) Jac_y(t) Jac_y(t+1) -res]
		jacobian = [jacobian(:,iz) , -d1];

		ic = 1:ny;
		icp = iyp;
				
		%          Jac_now      \ [ Jac_frwd res ]
		c(ic,:) = jacobian(:,is)\jacobian(:,isf1) ;
				
				
		%% NEXT PERIODS
		for it_ = (1+t+(0:(Tmax-2)))
			z = [ endogenous_variables_paths(iyp,it_-1) ; endogenous_variables_paths(:,it_) ; endogenous_variables_paths(iyf,it_+1)];
%			[d1,jacobian] = feval(dynamic_model,z,spfm_exo_simul, M.params, oo.steady_state, 2+it_-t);
            [d1, jacobian] = dynamicmodel(z, spfm_exo_simul, M.params, oo.steady_state, 2+it_-t);
 			jacobian = [jacobian(:,iz) , -d1];
			jacobian(:,[isf nrs]) = jacobian(:,[isf nrs])-jacobian(:,isp)*c(icp,:);
			ic = ic + ny;
			icp = icp + ny;
			c(ic,:) = jacobian(:,is)\jacobian(:,isf1);
        end

		% Terminal condition is Y_{T}=Y^{\star}
        c = back_subst_lbj(c, ny, iyf, Tmax);
		
		% F(Y^k) = dF/dY*(Y(k+1)-Y(k))
		% set expectation paths
		endogenous_variables_paths(:,it_init+(0:Tmax-1)) = endogenous_variables_paths(:,it_init+(0:Tmax-1))+c;

				
		err = max(max(abs(c)));
				
		info.iterations.error(iter) = err;
		if iter>1
			error_growth = error_growth + (info.iterations.error(iter)>info.iterations.error(iter-1));
		end
		if isnan(err) || error_growth>3
			last_line = iter;
			break
		end
		if err < options_.dynatol.f
			stop = 1;
			break
		end
    end
    
end

