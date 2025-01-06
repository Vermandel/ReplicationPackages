function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = call_CARA.static_resid_tt(T, y, x, params);
end
residual = zeros(60, 1);
    residual(1) = (y(32)) - (y(32)*T(1)+T(2)*T(4)^(1-params(17)));
    residual(2) = (y(1)) - (T(5));
    residual(3) = (y(8)) - (T(7)*T(8)-T(6)*y(31)*(1-params(12)));
    residual(4) = (y(31)) - (T(5)*T(9));
residual(5) = y(30);
    residual(6) = (y(14)) - (exp((-params(25))*(y(4)*params(31)-params(3)))*params(10)*y(16)*T(10)*T(11));
    residual(7) = (y(10)) - (y(10)*y(1)*params(39)*(1+y(14)*(1-params(14)))+y(1)*y(11)*(1-params(28))*y(13));
    residual(8) = (y(12)) - (1+y(14));
    residual(9) = (y(11)*params(28)*T(12)) - (params(14)*y(10)*y(14)/y(6));
    residual(10) = (y(18)) - ((1-params(28))*y(13)/params(21));
    residual(11) = (y(15)) - (y(13)*T(13)-y(18)*params(21)-y(7)-y(59)*y(17));
    residual(12) = (y(11)) - (T(13)-y(59)*(1-y(58))*params(4));
    residual(13) = (1) - (y(9)*params(23)*T(15));
    residual(14) = (y(9)) - (y(1)*params(39)*y(9)*T(17)+T(12)*y(1)*y(11)*params(28));
    residual(15) = (y(12)*y(3)) - (y(3)*(1-params(19))+y(3)*(params(23)/(1-params(13))*T(16)+params(24)));
    residual(16) = (y(13)) - (y(16)*exp((-params(25))*(y(4)*params(31)-params(3)))*params(6)*params(33)*T(18)*T(19));
    residual(17) = (y(3)) - (y(6)+y(5));
    residual(18) = (y(4)*params(31)) - (params(31)*y(4)*(1-params(27))+params(5)*(y(17)+params(9)));
    residual(19) = (y(17)) - (y(13)*params(4)*(1-y(58))*params(41));
    residual(20) = (y(13)*T(13)) - (y(2)+y(7)+params(29)*params(42));
    residual(21) = (log(y(16))) - (log(y(16))*params(18)+params(30)*x(1));
    residual(22) = (y(20)) - (params(39)*y(12)*y(1)*(y(15)+y(20)));
    residual(23) = (1/(1+y(19))) - (y(1));
    residual(24) = (y(21)) - ((T(20)-(1+y(19)))*100);
    residual(25) = (y(26)) - (y(19)*100);
    residual(26) = (y(25)) - (y(2)*100/y(13));
    residual(27) = (y(27)) - (y(6)*100/y(5));
    residual(28) = (y(22)) - (100*log(y(12)));
    residual(29) = (y(23)) - (log(y(12)));
    residual(30) = (y(24)) - (100*log(y(12)));
    residual(31) = (y(29)) - (y(2)*params(1)*1000/(y(4)*params(31)));
    residual(32) = (y(28)) - (y(12)*y(1)*T(21));
    residual(33) = (y(59)*params(4)) - (params(43)*params(44)*y(58)^(params(44)-1));
    residual(34) = (params(5)*y(60)) - (y(59));
    residual(35) = (y(60)) - (y(12)*y(1)*T(22));
    residual(36) = (y(33)) - (y(32)*params(32));
residual(37) = y(34);
    residual(38) = (y(35)) - (100*log(y(12)));
residual(39) = y(36);
    residual(40) = (y(37)) - (y(4)*params(31));
    residual(41) = (y(47)) - (T(20)-1);
    residual(42) = (y(54)) - ((T(20)-(1+y(19)))*100);
    residual(43) = (y(42)) - (y(1));
    residual(44) = (y(43)) - (y(1));
    residual(45) = (y(44)) - (T(23));
    residual(46) = (y(45)) - ((y(1)-y(43))^3);
    residual(47) = (y(46)) - (log(y(20)/y(15)));
    residual(48) = (y(48)) - (T(25)/T(24));
    residual(49) = (y(38)) - (y(1)/y(12)*(1+y(38)));
    residual(50) = (y(49)) - ((1+y(38))/y(38));
    residual(51) = (y(50)) - (100*(y(49)-(1+y(19))));
    residual(52) = (y(52)) - (y(28));
    residual(53) = (y(39)) - (log(y(28)/(y(28))));
    residual(54) = (y(40)) - (log(y(2)/(y(2))));
    residual(55) = (y(41)) - (log(y(13)/(y(13))));
    residual(56) = (y(51)) - (100*y(40));
    residual(57) = (y(56)) - (log(T(26)));
    residual(58) = (y(57)) - (log(1+y(19)));
residual(59) = y(55);
    residual(60) = (y(53)) - (log(y(59)));

end
