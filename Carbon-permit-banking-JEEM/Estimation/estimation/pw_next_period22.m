function [ynew,vv,rbool,ghu] = pw_next_period22(ee,yhat,pw,vv,M_,sbool,Q)
%PW_NEXT_PERIOD Summary of this function goes here
%   Detailed explanation goes here

	qMAX = 150; % max number of forward periods to explore
	iterMAX= 150;
	Rg = pw.Nr;
	
	order_var = pw.oo_.dr.order_var;
	
	% get linear case
	ynew 		= pw.PP*yhat+pw.QQ*ee;
	

	rbool  = eval(sbool);
	ghu    = pw.oo_.dr.ghu;

	%% identify the regime
	[~, idx]=ismember(rbool',pw.combinations,'rows');
	%% force discarded observable to be zero
	%now_obs = transpose(obs(i-1,:).*pw.obs_per_regime(idx,:));
	%% compute the currrent selection matrix
	now_Q = Q;
	% discard some endogenous 
	%now_Q(find(pw.shock_per_regime(idx,:)==0),:) = 0;
	%% selection matrix of shocks
	% identity matrix except for disabled shocks
	now_S = eye(pw.M_.exo_nbr);
	%now_S(1:1+size(now_S,1):end) = pw.shock_per_regime(idx,:);

	
	if sum(rbool)>0 % then some regimes bind
       %% step 1 :
		% get linear duration (inital guess)
		b = rbool;
		ib = 0; % counter
		tbool = rbool;
		while sum(b) > 0
				ib = ib + 1;
				ynew=pw.PP*(ynew);
				tbool = [ tbool   eval(sbool) ];
				b =  sum(tbool(:,end));
				if ib > qMAX % then the two regimes are not the same
					ynew=ynew+1000;
					ghu=[];
					q=[];
					% stop here
					return;
				end
		end
		q=ib;

		if ib > qMAX
			ynew=ynew+1000;
				ghu=[];
				q=[];
				% stop here
				return;
		end
		tbool = tbool(:,1:ib); % take away last empty row

		
		iter = 0;
		% step 2: guess and try
		violation = 1;
		while violation == 1

			iter = iter +1;

			if iter > iterMAX
				violation=0;
			end
			
			[vv] = update_dr(vv,q,tbool,pw);
			% from here we have the value of P(t), Q(t) and R(t)
			
			% now let's see if PQ is consistent
			% keep prev history and use state dependent policy rule
			y_test = [yhat zeros(M_.endo_nbr,q)];
			for i1 = 1:q
				% update name (in reverse order here!)
				% we use the longest PQ
				vname = sprintf(['.' repmat('%i',[1 Rg])],flip(tbool(:,i1:end),2));
				idn = strmatch(vname(2:end),vv.tree,'exact');
				% propagation PP(t)*x(-1)
				ynew = vv.PP(:,:,idn)*y_test(:,i1);
				% constant term RR(t)
				ynew = ynew + vv.RR(:,:,idn);
%				vname
				if i1 ==1 % then add the current shock
					%ynew = ynew + eval(['vv.QQ' vname '.m * (now_S*ee) ']);
					ynew = ynew + vv.QQ(:,:,idn)*(now_S*ee);
					% keep the value of the state dependent GHU matrix
					ghu = vv.QQ(:,:,idn);
				end
				%save result
				y_test(:,i1+1) = ynew;
				% check whether the path is consistent with guess
				if i1 > 1 && sum(eval(sbool)-tbool(:,i1))~=0
					%error('err..');
					% if not consistent
					% we are not sure of the outcome in the next periods
					% so delete what's next
					y_test = y_test(:,1:i1);
					q = i1-1;
					tbool = tbool(:,1:(i1-1));
					% stop current loop here
					break; 
				end
			end % for loop ends
			
			% check whether condition is violated at next period
			ynew = pw.PP*y_test(:,end);
			tempbool = eval(sbool);
			
			if sum(tempbool) == 0 || q==qMAX
				% termination condition met
				% leave the loop plz
				violation = 0;
			else
				
				% not finished
				q = q +1;
				% update the booleans
				tbool = [tbool tempbool];
				
				% if still binding now
				% check if we can improve guess
				if sum(tempbool) > 0
					% guess and try again to find new normal regime
					while sum(tempbool) ~= 0
						% check whether condition is violated at next period
						ynew = pw.PP*ynew;
						tempbool = eval(sbool);
						% update the booleans
						if sum(tempbool) > 0
							tbool = [tbool tempbool];
							q=q+1;
						end
						if q > qMAX % then the two regimes are not the same
							ynew=ynew+1000;
							ghu=[];
							q=[];
							% stop here
							return;
						end
					end
				end
			end
			
		end
	if q<qMAX && iter<=iterMAX && size(y_test,2) > 1
		% get last iteration
		ynew=y_test(:,2);
	else
		ynew=yhat*0+10000;
	end
	end
%ynew(pw.oo_.dr.inv_order_var(2))+pw.oo_.dr.ys(2)
	%[round(pw.oo_.dr.ys(17)+ynew(pw.oo_.dr.inv_order_var(17))-0.000000001,8)	rbool]
end

