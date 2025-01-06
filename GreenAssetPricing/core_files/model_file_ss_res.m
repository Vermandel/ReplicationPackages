function [res] = model_file_ss_res(Y,deltaX,eta1,eta2,beta,sigma,squig,alpha,nY_sts,nH_sts,deltaK,zetaY,kappa,scale_y,iotaX,iotaE,thetaEZ,gamma,y,x,phi,iH,Sg,nu,xfactor,chi,vtax)
%MODEL_FILE_SS_RES Summary of this function goes here
%   Detailed explanation goes here
% load parameters
Y           = real(Y);
gamma       = Y(1);
mu          = Y(2);
zetaY_tilda = Y(3);

% compute auxiliary
btilda     	= beta*(gamma^(thetaEZ-1));
vE  		= eta2*eta1*(mu^(eta2-1))/iotaE;
varrho      = 1*(1-(eta1)*(mu^eta2))-eta1*eta2*(mu^(eta2-1))*(1-mu);
kYyratio   	= btilda*alpha*varrho/(1*(1-btilda*(1-deltaK)));
y          	= ((scale_y*(zetaY_tilda))^(1/(1-alpha)))*((kYyratio)^(alpha/(1-alpha)))*nY_sts; 
x			= (iotaX*((1-mu)*iotaE*y*(1+chi)))/(xfactor*deltaX);       
iH          = gamma-1;
phi         =  varrho*(1-alpha)*y/(1/btilda-(1 + (1-squig)*iH));  
kY     		= kYyratio*y;
kH          = phi*squig*(iH/(varrho*alpha*(y/kY) ));
iT     	   	= (gamma-(1-deltaK))*(kY+kH);
c      	   	= (1-(eta1)*(mu^eta2))*y - iT - Sg*y;
vX      	= ( vtax + nu*c/(xfactor*x) + zetaY*varrho*y +  zetaY*phi*iH  )/(1/(gamma*btilda)-(1-deltaX));

% get residual vector
res = real([ 
	vX*iotaX  - vE;                 
	exp(-zetaY*xfactor*x) -  zetaY_tilda;
	iH  - (zetaY_tilda)*kappa*(kH^squig)*(nH_sts^(1-squig));
]);


end
