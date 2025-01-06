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
    T = call_EZ.static_resid_tt(T, y, x, params);
end
residual = zeros(63, 1);
    residual(1) = (y(33)*params(32)) - (T(5)^(1/params(11)));
    residual(2) = (params(7)*y(1)) - ((y(33)*params(32))^(1-params(17)));
    residual(3) = (y(9)) - (T(7)*(T(8)-y(32)*(1-params(12))));
    residual(4) = (y(32)) - (T(12)*(T(8)+y(32)*params(12)));
    residual(5) = (y(2)) - (T(17)*T(20));
    residual(6) = (y(31)*y(13)) - (y(31)*params(12)+y(3)*params(40)*(1-params(12))/T(1));
    residual(7) = (y(15)) - (exp((-params(25))*(y(5)*params(31)-params(3)))*params(10)*y(17)*T(21)*T(22));
    residual(8) = (y(11)) - (y(11)*y(2)*params(39)*(1+y(15)*(1-params(14)))+y(2)*y(12)*(1-params(28))*y(14));
    residual(9) = (y(13)) - (1+y(15));
    residual(10) = (y(12)*params(28)*T(23)) - (params(14)*y(11)*y(15)/y(7));
    residual(11) = (y(19)) - ((1-params(28))*y(14)/params(21));
    residual(12) = (y(16)) - (y(14)*T(24)-y(19)*params(21)-y(8)-y(60)*y(18));
    residual(13) = (y(12)) - (T(24)-y(60)*(1-y(59))*params(4));
    residual(14) = (1) - (y(10)*params(23)*T(26));
    residual(15) = (y(10)) - (y(2)*params(39)*y(10)*T(28)+T(23)*y(2)*y(12)*params(28));
    residual(16) = (y(13)*y(4)) - (y(4)*(1-params(19))+y(4)*(params(23)/(1-params(13))*T(27)+params(24)));
    residual(17) = (y(14)) - (y(17)*exp((-params(25))*(y(5)*params(31)-params(3)))*params(6)*params(33)*T(29)*T(30));
    residual(18) = (y(4)) - (y(7)+y(6));
    residual(19) = (y(5)*params(31)) - (params(31)*y(5)*(1-params(27))+params(5)*(y(18)+params(9)));
    residual(20) = (y(18)) - (y(14)*params(4)*(1-y(59))*params(41));
    residual(21) = (y(14)*T(24)) - (y(3)+y(8)+params(29)*params(42));
    residual(22) = (log(y(17))) - (log(y(17))*params(18)+params(30)*x(1));
    residual(23) = (y(21)) - (params(39)*y(13)*y(2)*(y(16)+y(21)));
    residual(24) = (1/(1+y(20))) - (y(2));
    residual(25) = (y(22)) - ((T(31)-(1+y(20)))*100);
    residual(26) = (y(27)) - (y(20)*100);
    residual(27) = (y(26)) - (y(3)*100/y(14));
    residual(28) = (y(28)) - (y(7)*100/y(6));
    residual(29) = (y(23)) - (100*log(y(13)));
    residual(30) = (y(24)) - (log(y(13)));
    residual(31) = (y(25)) - (100*log(y(13)));
    residual(32) = (y(30)) - (y(3)*params(1)*1000/(y(5)*params(31)));
    residual(33) = (y(29)) - (y(13)*y(2)*T(32));
    residual(34) = (y(60)*params(4)) - (params(43)*params(44)*y(59)^(params(44)-1));
    residual(35) = (params(5)*y(61)) - (y(60));
    residual(36) = (y(61)) - (y(13)*y(2)*T(33));
    residual(37) = (y(34)) - (y(33)*params(32));
residual(38) = y(35);
    residual(39) = (y(36)) - (100*log(y(13)));
residual(40) = y(37);
    residual(41) = (y(38)) - (y(5)*params(31));
    residual(42) = (y(48)) - (T(31)-1);
    residual(43) = (y(55)) - ((T(31)-(1+y(20)))*100);
    residual(44) = (y(43)) - (y(2));
    residual(45) = (y(44)) - (y(2));
    residual(46) = (y(45)) - (T(34));
    residual(47) = (y(46)) - ((y(2)-y(44))^3);
    residual(48) = (y(47)) - (log(y(21)/y(16)));
    residual(49) = (y(49)) - (T(36)/T(35));
    residual(50) = (y(39)) - (y(2)/y(13)*(1+y(39)));
    residual(51) = (y(50)) - ((1+y(39))/y(39));
    residual(52) = (y(51)) - (100*(y(50)-(1+y(20))));
    residual(53) = (y(53)) - (y(29));
    residual(54) = (y(40)) - (log(y(29)/(y(29))));
    residual(55) = (y(41)) - (log(y(3)/(y(3))));
    residual(56) = (y(42)) - (log(y(14)/(y(14))));
    residual(57) = (y(52)) - (100*y(41));
    residual(58) = (y(57)) - (log(T(37)));
    residual(59) = (y(58)) - (log(1+y(20)));
    residual(60) = (y(56)) - (log(y(31)/(y(31))));
    residual(61) = (y(54)) - (log(y(60)));
    residual(62) = (y(62)) - (y(5));
    residual(63) = (y(63)) - (y(31));

end
