function [y_lf, x_lf, zetaY_tilda_lf, phi_lf, kH_lf, gamma_lf] = EZWH_ss_lf_welfcost(scale_y,zetaY,alpha,xfactor,deltaX,vartheta,kYyratio_lf,kappa,squig,nY_sts,nH_sts,beta,thetaEZ,y_lf, x_lf, zetaY_tilda, phi_lf, kH_lf, gamma_lf,iotaX,iotaE,chi,deltaK)
%MODEL_SS_LF_WELFCOST Summary of this function goes here
%   Detailed explanation goes here

if vartheta ~= 1
    % y_lf x_lf zetaY_tilda phi_lf kH_lf gamma_lf kYyratio_lf
    fun = @(x) [
    %y_lf      	- ((scale_y*vartheta)^(1/(1-alpha)))*((zetaY_tilda)^(1/(1-alpha)))*((kYyratio_lf)^(alpha/(1-alpha)))*nY_sts; 
    x(1)      	- ((scale_y*vartheta)^(1/(1-alpha)))*((x(3))^(1/(1-alpha)))*((x(7))^(alpha/(1-alpha)))*nY_sts; 
    %x_lf     	- ( iotaX * (1+chi)*iotaE  )/deltaX*y_lf/xfactor;
    x(2)     	- ( iotaX * (1+chi)*iotaE  )/deltaX*x(1)/xfactor;
    %zetaY_tilda - exp(-zetaY*(x*xfactor));
    x(3)        - exp(-zetaY*(x(2)*xfactor));
    %phi_lf    	- (btilda_lf*(1-alpha)*y_lf)/(1-btilda_lf*(1+(1-squig)*(gamma_lf-1)));
    x(4)    	- (beta*(x(6)^(thetaEZ-1))*(1-alpha)*x(1))/(1-beta*(x(6)^(thetaEZ-1))*(1+(1-squig)*(x(6)-1)));
    %kH_lf     	- (phi_lf/varrho_lf)*(squig/alpha)*kYyratio_lf*( (gamma_lf-1));
    x(5)     	- (x(4))*(squig/alpha)*x(7)*( (x(6)-1));
    %kappa      	- ( (gamma_lf-1))/((zetaY_tilda)*(kH_lf^squig)*(nH_sts^(1-squig)));
    kappa      	- ( (x(6)-1))/((x(3))*(x(5)^squig)*(nH_sts^(1-squig)));
    %kYyratio_lf = varrho_lf*btilda_lf*alpha/(q*((1-btilda_lf*(1-deltaK))));
    x(7) -          beta*x(6)^(thetaEZ-1)*alpha/((1-beta*(x(6)^(thetaEZ-1)*(1-deltaK))));
    ];
    
    x0 = [y_lf x_lf zetaY_tilda phi_lf kH_lf gamma_lf kYyratio_lf];
    rez = fsolve(fun,x0);
    y_lf     = rez(1);
    x_lf     = rez(2);
    zetaY_tilda_lf = rez(3);
    phi_lf   = rez(4);
    kH_lf    = rez(5);
    gamma_lf = rez(6);
    kYyratio_lf = rez(7);
else
    zetaY_tilda_lf = zetaY_tilda;
end


end

