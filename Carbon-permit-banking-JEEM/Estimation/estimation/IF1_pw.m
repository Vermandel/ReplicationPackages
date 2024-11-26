function [llks,ee,res,y_,stop_dist] = IF1_pw(obs,e_obj,oo,M,options_,ee0)
	%DSGE_INVERT_2O Summary of this function goes here
	%   Detailed explanation goes here
	
	% Get sample size
	N = length(e_obj.id.obs);
	T		= size(obs,1);
	% get selection matrix
	Q		= e_obj.Q;
	
	if nargin < 6
		ee0 = zeros(T,N);
	end
	if size(ee0,1) ~= T || size(ee0,2) ~= N
		ee0 = zeros(T,N);
	end
	
	%
	options  = optimoptions('fsolve','Display','off');

	options = optimset('Display','off','Algorithm','Levenberg-Marquardt','TolFun',1e-5,'TolX',1e-5);

	

	k2			= oo.dr.kstate(find(oo.dr.kstate(:,2) <= M.maximum_lag+1),[1 2]);
	k2			= k2(:,1)+(M.maximum_lag+1-k2(:,2))*M.endo_nbr;
	order_var 	= oo.dr.order_var;
	y0			= oo.dr.ys;
	ee 			= zeros(T,N);
	llks 		= zeros(T,1);
	res			= 0;
	y_ 			= zeros(size(y0,1),T+M.maximum_lag);
	y_(:,1) 	= y0;

	% inverse the varcov matrix
	iSigE = inv(M.Sigma_e);
	ighu  = inv(Q*oo.dr.ghu(oo.dr.inv_order_var,:));
	
	if isinf(isinf(ighu))
		error('Error matrix on shocks not invertible.');
	end
	
	exflag = 1;
	persistent v;
	%%% PIECE-WISE SOLUTION HOOK
	pw = e_obj.pw;
	pw.oo_ = oo;
	pw.M_  = M;
	ynew = y0;
	xx=0;
	
	%% first include zeros not non-state variables in standard regime
	%% and load each submodels
	%% (this is a copy paste from build_pw_struct.m)
	for i1=1:size(pw.combinations,1)
		
		fname = sprintf('%i',pw.combinations(i1,:));
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
		else % alternative regime model
			eval(['[pw.FF' fname ',pw.GG' fname ',pw.HH' fname ',pw.MM' fname ',pw.CC' fname ',pw.Jac' fname '] = read_pq(myss,''' M.fname '_' fname ''',pw.oo_,pw.M_,pw.options_);'])
		end
	end

	Rg = pw.Nr;
	
	% Initialize the decision rule structure tree
	vv = [];
	% initialize with the normal regime:
	reg_name=sprintf([repmat('%i',[1 Rg])],pw.combinations(1,:));
	vv.tree=char(reg_name);
	idx     = 1;
	vv.SS(:,:,idx) =   eval(['-inv(pw.FF' reg_name '*pw.PP+pw.GG' reg_name ')']);
	vv.PP(:,:,idx) =   eval(['vv.SS(:,:,idx)*pw.HH' reg_name ]);
	vv.QQ(:,:,idx) =   eval(['vv.SS(:,:,idx)*pw.MM' reg_name ]);
	vv.RR(:,:,idx) =   eval(['vv.SS(:,:,idx)*(pw.CC' reg_name '+pw.FF' reg_name '*zeros(pw.M_.endo_nbr,1))']);
		
		
	% regimes matrix
	rbool   = zeros(Rg,T+1);
	sbool   = '[';
	for i1 = 1:Rg
		if i1 > 1
			sbool   = [sbool ';'];
		end
		sbool   = [sbool eval(['read_bind(pw.c' num2str(i1) '_bind,M)']) ];
	end
	sbool   = [sbool ']'];
	
%	save piie;error('gg')
	% start inversion	
	for i = 2:T+M.maximum_lag
			
		% usual step would be this:
        %y_(order_var,i) = oo.dr.ys(order_var)+dr.ghx*yhat + dr.ghu*epsilon(:,i-1);
		
		% reversion implies that:
		% obs(:,i-1) = Q*oo.dr.ys+Q*(dr.ghx*yhat + dr.ghu*epsilon(:,i-1));
		% (obs(:,i-1)-Q*(oo.dr.ys+dr.ghx*yhat))*(Q*oo_.dr.ghu)^-1 =  *epsilon(:,i-1))
		% get previous state variables x_t = X_t - Xss;
		yhat = y_(order_var(k2),i-1)-y0(order_var(k2));
		yhat2 = oo.dr.ghx*yhat;
		% compute shock using inversion
		ee(i-1,:) = transpose(ighu*((obs(i-1,:)'-Q*(oo.dr.ys+yhat2(oo.dr.inv_order_var)))));
		ghu  	= pw.oo_.dr.ghu(:,:);
		
		% new data
		ynew = oo.dr.ghx*yhat + oo.dr.ghu*(ee(i-1,:)');
		% check if violates some conditions
		rbool(:,i)  = eval(sbool);

		if sum(rbool(:,i))>0 % then some regimes bind
			%save piie; error('STOP');
			%% identify the regime
			[~, idx]=ismember(rbool(:,i)',pw.combinations,'rows');
			%% force discarded observable to be zero
			now_obs = transpose(obs(i-1,:).*pw.obs_per_regime(idx,:));

			% take all variables (including non-states)
			yhat3 = y_(order_var,i-1)-y0(order_var);

			% guess the initial values of the shocks ee(i-1,:)
			% we compute the duration of the constraint using the linear solution
			% then we compute the current decision rules PP(now), RR(now), QQ(now)
			% and invert this equation as iniatial guess: 
			%    	obs = Q*(SS + PP(now)*Y(now-1) + QQ(now)*ee(i-1,now) + RR(now))
			% then : ee(i-1,now) = (Q*QQ(now))^-1 * ( obs - Q*(SS + PP(now)*Y(now-1) + RR(now)))
			[epsilon,vv,ghu,e0,exflag,ynew] = IF_pw_next(vv,pw,obs(i-1,:),sbool,Q,yhat3,ynew);

			if exflag~=1
				%disp('lol')
				
			end
			% guess and try inversion failed
			% use solvers
			%if exflag~=1
				% Solve the observation equation using solver
				%[epsilon,exflag]=csolve(@(x) [obs(i-1,:)'-Q*oo.dr.ys - Q(:,order_var)*pw_next_period22(x,yhat3,pw,vv,M,sbool,Q)],0*(now_S*ee(i-1,:)'),[],1e-6,20);			
			%	[epsilon,exflag] = csolve_grad('test',e0',1e-4,40,yhat3,pw,vv,M,sbool,Q,obs(i-1,:)');
			%	exflag=exflag+1;
				% if csolve did not make it too
			%	if exflag ~= 1
					% try with fsolve (more powerful, but more computational)
			%		[epsilon,fval,exflag] = fsolve(@(x) [obs(i-1,:)'-Q*oo.dr.ys - Q(:,order_var)*pw_next_period22(x,yhat3,pw,vv,M,sbool,Q)],e0',options);
			%	end
			% try with the computed shocks
			[ynew,vv,rbool(:,i),ghu]=pw_next_period22(epsilon,yhat3,pw,vv,M,sbool,Q);
			
			nres= sqrt(sum((obs(i-1,:)'-Q*pw.oo_.dr.ys - Q(:,pw.oo_.dr.order_var)*ynew).^2));
			if i ==T+M.maximum_lag&1==0
			[obs(i-1,:)' Q*pw.oo_.dr.ys + Q(:,pw.oo_.dr.order_var)*ynew]
			end
			if nres >0.000001
				res = res+nres;
			end
			
			if exflag ~= 1 || isempty(ghu) % if inversion failed
				lambda = 1;
				res=10;
			else	
				% if inversion worked:
				% select the current regime
				[~, idx]=ismember(rbool(:,i)',pw.combinations,'rows');
				%% compute the currrent selection matrix
				now_Q = Q;
				% discard some endogenous 
				now_Q(find(pw.shock_per_regime(idx,:)==0),:) = 0;
				%% selection matrix of shocks
				% identity matrix except for disabled shocks
				now_S = eye(M.exo_nbr);
				now_S(1:1+size(now_S,1):end) = pw.shock_per_regime(idx,:);
				% select shock
				ee(i-1,:) = (now_S*epsilon)';
				% standard gradient
				lambda = Q*ghu(oo.dr.inv_order_var,:);
				% remove both obs and shocks
				lambda = lambda(find(pw.obs_per_regime(idx,:)==1),find(pw.shock_per_regime(idx,:)==1));
			end

			idshocks = find(pw.shock_per_regime(idx,:));

			% update likelihood
			llks(i-1) = - 1/2*ee(i-1,idshocks)*inv(M.Sigma_e(idshocks,idshocks))*ee(i-1,idshocks)' - log(abs(det(lambda)));
		else % linear solution is enough
			% determine the jacobian
			lambda = Q*oo.dr.ghu(oo.dr.inv_order_var,:);
			% update likelihood
			llks(i-1) = - 1/2*ee(i-1,:)*iSigE*ee(i-1,:)' -  log(abs(det(lambda)));
		end

		% save the new value of endogenous
		y_(order_var,i) = oo.dr.ys(oo.dr.order_var) + ynew;
		
        % check residuals
        current_gap = sum(abs(obs(i-1,:)'-Q*y_(:,i)));
        if current_gap > 1e-5
            res = res+current_gap;
        end
        
        % [i-1 res]
		% if residual too high, stop inversion here
		if res > 0.1
			% compute residuals
			%res = res+sum(abs(obs(i-1,:)'-Q*y_(:,i)));
			break;
		end
	end

	stop_dist= T+M.maximum_lag-i;
    res = res+stop_dist/T;
	if stop_dist >0
		fprintf('(%i)',i-1)
	end
end

