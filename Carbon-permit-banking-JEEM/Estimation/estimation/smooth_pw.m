[llk,ee,res,yy] = dsge_llk(theta,obs,e_obj,oo_,M_,options_,bayestopt_);
%% update the model's
M_.params(e_obj.theta.id(e_obj.theta.is_pm)) = theta(e_obj.theta.is_pm);
M_.Sigma_e = zeros(M_.exo_nbr);
M_.Sigma_e(1:length(e_obj.idu.shocks)+1:end)=theta(e_obj.theta.is_sd).^2;
M_.Sigma_e(e_obj.theta.id(e_obj.theta.is_sd),e_obj.theta.id(e_obj.theta.is_sd)) = M_.Sigma_e;
% update PQ (normal model)
[oo_.dr.ys, M_.params] = feval([M_.fname '.steadystate'],oo_.dr.ys, zeros(M_.exo_nbr,1), M_.params);
oo_.dr = stochastic_solvers(oo_.dr,0,M_,options_,oo_);
% update pw (piece-wise)
%%% PIECE-WISE SOLUTION HOOK
pw = e_obj.pw;
pw.oo_ = oo_;
pw.M_  = M_;

	
%% first include zeros not non-state variables in standard regime
%% and load each submodels
%% (this is a copy paste from build_pw_struct.m)
for i1=1:size(pw.combinations,1)
		
		fname = sprintf('%i',pw.combinations(i1,:))
		if i1 == 1 % if standard model

			%% Since state space (t+1,t,t-1) implies to have 3*n_variables
			%% the steady state must be replicated 3 times for each time horizon
			myss = nan(max(max(pw.M_.lead_lag_incidence)),1);
			for i2=1:size(myss,1)
				[~, col] = find(pw.M_.lead_lag_incidence == i2);
				myss(i2) = pw.oo_.dr.ys(col);
			end


			% GET POLICY FUNCTION
			nvar 	= pw.M_.endo_nbr;
			nstate 	= length(pw.oo_.dr.state_var);
			ghx  	= pw.oo_.dr.ghx(:,:);
			ghu  	= pw.oo_.dr.ghu(:,:);
			pw.QQ = ghu;
			% build PP
			jj = 0;
			for i1 = 1:nvar
				% if state variable, pick in ghx
				if sum( pw.oo_.dr.order_var(i1) ==  pw.oo_.dr.state_var)
					jj = jj + 1;
					PPnew = ghx(:,jj);
				else
					PPnew = zeros(nvar,1);
				end
				if i1 == 1
					pw.PP = PPnew;
				else
					pw.PP = [pw.PP PPnew];
				end
			end
		end
		% alternative regime model
		eval(['[pw.FF' fname ',pw.GG' fname ',pw.HH' fname ',pw.MM' fname ',pw.CC' fname ',pw.Jac' fname '] = read_pq(myss,''' M_.fname '_' fname ''',pw.oo_,pw.M_,pw.options_);'])
		
end
	
% simulate the non-linear model
[y_pw,y_lin] =  simul_pw_mat(pw,ee');

Tsample = size(ee,1);
eval(['y_obs = e_obj.Q*y_pw(2:end,:)'';'])
Nobs = size(dataset_.name,1);
figure
for i1=1:Nobs
	subplot(ceil(sqrt(Nobs)),ceil(sqrt(Nobs)),i1)
	plot(1:Tsample,dataset_.data(:,i1),...
		 1:Tsample,y_obs(i1,:),'or');
    title(M_.endo_names{e_obj.id.obs(i1)})
end


figure;
for i1 = 1:length(e_obj.id.shocks)
    subplot(ceil(sqrt(length(e_obj.id.shocks))),ceil(sqrt(length(e_obj.id.shocks))),i1)
    plot((1:Tsample)',ee(:,i1))
    title(M_.exo_names{e_obj.id.shocks(i1)})
end
