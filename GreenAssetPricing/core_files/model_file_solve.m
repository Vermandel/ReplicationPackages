function [gamma,mu,zetaY_tilda] = model_file_solve(deltaX,eta1,eta2,beta,sigma,squig,alpha,nY_sts,nH_sts,deltaK,zetaY,kappa,scale_y,iotaX,iotaE,thetaEZ,gamma,y,x,phi,iH,Sg,nu,xfactor,chi,vtax)

%MODEL_FILE_SOLVE Summary of this function goes here
%   Detailed explanation goes here

btilda     	= beta*(gamma^(thetaEZ-1));
vX      	= ( vtax+nu*0.6*y/(xfactor*x) + zetaY*y +  zetaY*phi*iH  )/(1/(gamma*btilda)-(1-deltaX));
mu     		= 1/2*(iotaX*iotaE/(eta2*eta1)*vX)^(1/(eta2-1));
x      		= (iotaX/(1-(1-deltaX)))*(iotaE*(1-mu)*y*(1+chi))/xfactor;

mus = 0:0.01:1;
res = zeros(size(mus));
for ix=1:length(mus)
       varrho     = 1*(1-(eta1)*(mus(ix)^eta2))-eta2*eta1*(mus(ix)^(eta2-1))*(1-mus(ix));
       kYyratio   = btilda*alpha*varrho/((1-btilda*(1-deltaK)));
       x          = (iotaX/(1-(1-deltaX)))*(iotaE*(1-mus(ix))*y*(1+chi))/xfactor;
       y          = ((scale_y*exp(-zetaY*xfactor*x))^(1/(1-alpha)))*((kYyratio)^(alpha/(1-alpha)))*nY_sts; 
       x          = (iotaX/(1-(1-deltaX)))*(iotaE*(1-mus(ix))*y*(1+chi))/xfactor;
       x0         = [gamma mus(ix) exp(-zetaY*xfactor*x)];
       rez        = model_file_ss_res(x0,deltaX,eta1,eta2,beta,sigma,squig,alpha,nY_sts,nH_sts,deltaK,zetaY,kappa,scale_y,iotaX,iotaE,thetaEZ,gamma,y,x,phi,iH,Sg,nu,xfactor,chi,vtax);
       res(ix)    = sum(abs(rez(1)));
end

idd = find(res==min(res));
mu  = mus(idd);
x   = (iotaX/(1-(1-deltaX)))*(iotaE*(1-mu)*y*(1+chi))/xfactor;
x0  = [gamma mu exp(-zetaY*xfactor*x)];



options = optimoptions('fsolve','Display','off','TolFun',1e-10,'MaxFunEvals',1e5,'Maxiter',1e4);%,'Algorithm','levenberg-marquardt');
xx = fsolve('model_file_ss_res',x0,options ,deltaX,eta1,eta2,beta,sigma,squig,alpha,nY_sts,nH_sts,deltaK,zetaY,kappa,scale_y,iotaX,iotaE,thetaEZ,gamma,y,x,phi,iH,Sg,nu,xfactor,chi,vtax);
gamma       = xx(1);
mu          = xx(2);
zetaY_tilda = xx(3);

if mu>1
    error('mu > 1')
end
if any(imag(xx)) || any(isnan(xx))
    error('Imaginary ss or nans')
end

res        = model_file_ss_res(xx,deltaX,eta1,eta2,beta,sigma,squig,alpha,nY_sts,nH_sts,deltaK,zetaY,kappa,scale_y,iotaX,iotaE,thetaEZ,gamma,y,x,phi,iH,Sg,nu,xfactor,chi,vtax);
%[abs(res) xx]
if any(abs(res) > 0.00001)
    error('residuals')
end