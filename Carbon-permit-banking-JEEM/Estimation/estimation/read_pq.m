function [FF,GG,HH,MM,res,Jac] = read_pq(myss,file_name,oo_,M_,options_)
%READ_PQ Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% JACOBIAN PLZ
nvar 	= M_.endo_nbr;
nstate 	= length(oo_.dr.state_var);
% Get now residuals + Jacobian Matrix from Dynare
[res,Jac] = eval([file_name '.dynamic(myss, zeros(1,M_.exo_nbr), M_.params, oo_.dr.ys, 1)']);   
% split the jacobian matrix
Jac_past   = Jac(:,1:(M_.npred+M_.nboth));
Jac_today  = Jac(:,(1+M_.npred+M_.nboth):(M_.npred+M_.nboth+M_.endo_nbr));
Jac_future = Jac(:,(1+M_.npred+M_.nboth+M_.endo_nbr):(end-M_.exo_nbr));
Jac_shocks = Jac(:,1+end-M_.exo_nbr:end);

%% CORRECTION 1: THE EXPECTATION JACOBIAN
% the jacobian matrix is empty for absent expected variable 
% replace these by 0 to have a full ranked matrix
dummy_expected = (M_.lead_lag_incidence(1,:)>0);
Jac_pastfull = zeros(nvar,nvar);
% parsing each variable and check whether it exists
i2 = 0;
for i1 = 1:length(dummy_expected)	
	if dummy_expected(i1) % then try
		i2 = i2 +1;
		Jac_pastfull(:,i1) = Jac_past(:,i2);
	else
		Jac_pastfull(:,i1) = zeros(1,M_.endo_nbr);
	end
end
%% CORRECTION 2: THE EXPECTATION JACOBIAN
% the jacobian matrix is empty for absent expected variable 
% replace these by 0 to have a full ranked matrix
dummy_expected = (M_.lead_lag_incidence(3,:)>0);
Jac_futfull = zeros(nvar,nvar);
% parsing each variable and check whether it exists
i2 = 0;
for i1 = 1:length(dummy_expected)	
	if dummy_expected(i1) % then try
		i2 = i2 +1;
		Jac_futfull(:,i1) = Jac_future(:,i2);
	else
		Jac_futfull(:,i1) = zeros(1,M_.endo_nbr);
	end
end

FF	= Jac_futfull(:,oo_.dr.order_var);
GG	= Jac_today(:,oo_.dr.order_var);
HH	= Jac_pastfull(:,oo_.dr.order_var);
MM	= Jac_shocks(:,:);

end

