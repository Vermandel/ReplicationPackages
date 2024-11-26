function [mc,k,c,n,y,pe,mu] = policy_anticipations_ss(e_pe,e_a,phi_1,phi_2,psi,nu,thet,p,A,beta,alpha,h,L,varphi,sigma,delta,chi,mc0,k0,c0,N,Y,pe0,MU,gamma)
%POLICY_ANTICIPATIONS_GV Summary of this function goes here
%   Detailed explanation goes here


	%	mc    k     c    n   y    pe   mu
	%   X(1) X(2) X(3) X(4) X(5) X(6) X(7)
 
	x0 = [mc0    k0     c0    N   Y    pe0   MU];
    options = optimoptions('fsolve','Display','none');
    func = @(X) [ 
				    -e_pe*phi_1*phi_2*X(7)^(phi_2-1)  + X(6)/1000*psi*X(5)^(-nu) ;
					-thet  + psi*(1-X(7))*X(5)^(1-nu);
					-X(1)    + (p-e_pe*phi_1*(X(7))^(phi_2)-X(6)/1000*(1-nu)*psi*(1-X(7))*X(5)^(1-nu)/X(5));
					-X(2)	+ alpha*X(5)*X(1)/(1/(beta*gamma^-sigma)-(1-delta));
					-X(3) 	+ X(5) - phi_1*X(7)^(phi_2)*X(5) - X(6)/1000*thet - (gamma - (1-delta))*X(2);
					-chi*X(4)^(1+varphi) + (X(3)*(1-h/gamma))^-sigma*(1-alpha)*X(5)*X(1);
					-X(5) + e_a*A*(L*X(4))^(1-alpha)*X(2)^alpha; 
 				];
    rez  = fsolve(func,x0,options);
    mc = rez(1);
    k  = rez(2);
    c  = rez(3);
    n  = rez(4);
    y  = rez(5);
    pe = rez(6);   
    mu = rez(7);
end

