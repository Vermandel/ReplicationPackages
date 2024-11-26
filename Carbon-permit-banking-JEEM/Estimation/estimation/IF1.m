function [llks,ee,res,y_,stop_dist] = IF1(obs,e_obj,oo,M,options_,ee0)
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
	options2 = optimset('Display','off','Algorithm','Levenberg-Marquardt');

	

	k2			= oo.dr.kstate(find(oo.dr.kstate(:,2) <= M.maximum_lag+1),[1 2]);
	k2			= k2(:,1)+(M.maximum_lag+1-k2(:,2))*M.endo_nbr;
	order_var 	= oo.dr.order_var;
	y0			= oo.dr.ys;
	ee 			= zeros(T,N);
	llks 		= zeros(T,1);
	res			= 0;
	y_ 			= zeros(size(y0,1),T+M.maximum_lag);
	y_(:,1) 	= y0;
%	y_(oo.dr.order_var,:) = y_;
	
	% inverse the varcov matrix
	iSigE = inv(M.Sigma_e);
	ighu  = inv(Q*oo.dr.ghu(oo.dr.inv_order_var,:));
	
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
		
		% from shock, compute new state variables for next period
		y_(order_var,i) = oo.dr.ys(oo.dr.order_var) + oo.dr.ghx*yhat + oo.dr.ghu*(ee(i-1,:)');

		% determine the jacobian
		lambda = Q*oo.dr.ghu(oo.dr.inv_order_var,:);

		% update likelihood
		llks(i-1) = - 1/2*ee(i-1,:)*iSigE*ee(i-1,:)' - log(abs(det(lambda)));
			
		% compute residuals
		res = res+sum(abs(obs(i-1,:)'-Q*y_(:,i)));
			
		% if residual too high, stop inversion here
		if res > 1
			break;
		end
	end
	
	stop_dist= T+M.maximum_lag-i;

end

