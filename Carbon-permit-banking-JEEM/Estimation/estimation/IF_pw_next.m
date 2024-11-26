function [epsilon,vv,ghu,e0,exflag,lynew] = IF_pw_next(vv,pw,obs,sbool,Q,yhat3,ynew)
%IF_PW_NEXT Summary of this function goes here
%   Detailed explanation goes here
	
	% provide inversion for a given duration guesss
	iterMAX = 150;
	qMAX    = 150;
	
	ghu = [];
	lynew = ynew+100; % set random value that generate large residuals that we will minimize
	exflag = 0;
	
	y0 = pw.oo_.dr.ys;
	Rg = pw.Nr;
	
	% reset initial guess from linear system
	epsilon = zeros(pw.M_.exo_nbr,1);
	e0      = epsilon;
	
	% guess idx of the system to discard some obs/shocks
	[~, idx]=ismember(eval(sbool)',pw.combinations,'rows');
	% identity matrix except for disabled shocks
	now_S = eye(pw.M_.exo_nbr);
	now_S(1:1+size(now_S,1):end) = pw.shock_per_regime(idx,:);

	%%save next22; error('STOP');
    
	% initialization
	violation = 1;
	iter = 0;
	% get Non-linearity duration
	q = 0; % counter
	tbool = eval(sbool);
	b = eval(sbool);
	while sum(b) > 0
		q = q + 1;
		ynew=pw.PP*(ynew);
		tbool = [ tbool   eval(sbool) ];
		b =  sum(tbool(:,end));
		if q > qMAX
			% way too long duration. Stop here.
			break;
		end
	end
	tbool = tbool(:,1:q); % take away last empty row
	
	if q <= qMAX
		
		while violation == 1

			% count plz
			iter = iter+1;
            
            if iter > iterMAX
                violation = 0;
            end
			% update decision rules with duration q
			[vv,VP,VQ,VR] = update_dr(vv,q,tbool,pw);
			
			% provide inversion given this duration
			% errors of prediction
			VVP = obs'-Q(:,pw.oo_.dr.order_var)*(y0(pw.oo_.dr.order_var)+VP*yhat3+VR);
			VVQ = Q(:,pw.oo_.dr.order_var)*VQ;
			
			% save the initial guess of shocks
			epsilon(find(pw.shock_per_regime(idx,:))) = transpose(inv(VVQ(find(pw.obs_per_regime(idx,:)),find(pw.shock_per_regime(idx,:))))*VVP(find(pw.obs_per_regime(idx,:))));
			
			if iter == 1
				% save first guess that can be re-used
				e0=epsilon';
			end
			
			% now let's see if PQ is consistent
			% keep prev history and use state dependent policy rule
			y_test = [yhat3 zeros(pw.M_.endo_nbr,q)];
			for i1 = 1:q
			
				% update name (in reverse order here!)
				% we use the longest PQ
				vname = sprintf(['.' repmat('%i',[1 Rg])],flip(tbool(:,i1:end),2));
				idn = strmatch(vname(2:end),vv.tree,'exact');
				
				% propagation PP(t)*x(-1)
				ynew = vv.PP(:,:,idn)*y_test(:,i1);
				% constant term RR(t)
				ynew = ynew + vv.RR(:,:,idn);
	%			vname
				if i1 ==1 % then add the current shock
					ynew = ynew + vv.QQ(:,:,idn)*(now_S*epsilon);
					% save best output
					ghu = vv.QQ(:,:,idn);
					lynew = ynew;
				end
				%save result
				y_test(:,i1+1) = ynew;
				% check whether the path is consistent with guess
				if i1 > 1 && sum(eval(sbool)-tbool(:,i1))~=0
%					error('err..');
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

			if (sum(tempbool) == 0 || iter==iterMAX) && size(y_test,2) > 1
				% termination condition met
				% keep new iteration   
				ynew = pw.PP*y_test(:,2);
                % leave the loop plz
				violation = 0;
			else		
				% not finished
				q = q +1;
				% update the booleans
				tbool = [tbool tempbool];
					
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
					if iter == iterMAX || q == qMAX % then the two regimes are not the same
						% stop here
						return; return;
					end
				end
				
			end
				
		end
		
	end
	
	if iter<iterMAX && q < qMAX && sum(obs'-Q*pw.oo_.dr.ys - Q(:,pw.oo_.dr.order_var)*lynew)<0.001
		% then last iteration was ok
		% inversion successful
		exflag = 1;
	else
		exflag = 0;
	end

end

