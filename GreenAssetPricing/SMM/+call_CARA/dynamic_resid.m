function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = call_CARA.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(60, 1);
    residual(1) = (y(44)) - (T(4)*y(86)+T(1)*T(3)^(1-params(17)));
    residual(2) = (y(13)) - (T(5)*y(20)/y(5));
    residual(3) = (y(20)) - (T(7)*T(8)-T(6)*y(43)*(1-params(12)));
    residual(4) = (y(43)) - (T(5)*T(33));
residual(5) = y(42);
    residual(6) = (y(26)) - (exp((-params(25))*(y(3)*params(31)-params(3)))*params(10)*y(28)*T(11)*T(12));
    residual(7) = (y(22)) - (params(39)*y(73)*y(78)*(1+(1-params(14))*y(81))+y(73)*y(79)*(1-params(28))*y(80));
    residual(8) = (y(24)) - (1+y(26));
    residual(9) = (params(28)*y(23)*y(25)/y(17)) - (params(14)*y(22)*y(26)/y(18));
    residual(10) = (y(30)) - ((1-params(28))*y(25)/params(21));
    residual(11) = (y(27)) - (y(25)*T(13)-y(30)*params(21)-y(19)-y(71)*y(29));
    residual(12) = (y(23)) - (T(13)-y(71)*(1-y(70))*params(4));
    residual(13) = (1) - (y(21)*params(23)*T(15));
    residual(14) = (y(21)) - (params(39)*y(73)*y(77)*(1-params(19)+T(16)*T(18)-params(23)*T(18)+params(24))+y(73)*y(79)*params(28)*T(19));
    residual(15) = (y(24)*y(15)) - (y(2)*(1-params(19))+y(2)*T(20));
    residual(16) = (y(25)) - (y(28)*exp((-params(25))*(y(3)*params(31)-params(3)))*params(6)*params(33)*T(21)*T(22));
    residual(17) = (y(2)) - (y(18)+y(17));
    residual(18) = (params(31)*y(16)) - (params(31)*y(3)*(1-params(27))+params(5)*(y(29)+params(9)));
    residual(19) = (y(29)) - (y(25)*params(4)*(1-y(70))*params(41));
    residual(20) = (y(25)*T(13)) - (y(14)+y(19)+params(29)*params(42));
    residual(21) = (log(y(28))) - (params(18)*log(y(10))+params(30)*x(it_, 1));
    residual(22) = (y(32)) - (params(39)*y(24)*y(73)*(y(82)+y(83)));
    residual(23) = (1/(1+y(31))) - (y(73));
    residual(24) = (y(33)) - ((y(7)*params(39)*(y(27)+y(32))/y(12)-(1+y(31)))*100);
    residual(25) = (y(38)) - (y(31)*100);
    residual(26) = (y(37)) - (y(14)*100/y(25));
    residual(27) = (y(39)) - (y(18)*100/y(17));
    residual(28) = (y(34)) - (100*(log(y(7))+log(y(25))-log(y(8))));
    residual(29) = (y(35)) - (log(y(7))+log(y(14))-log(y(1)));
    residual(30) = (y(36)) - (100*(log(y(7))+log(y(19))-log(y(4))));
    residual(31) = (y(41)) - (y(14)*params(1)*1000/(y(3)*params(31)));
    residual(32) = (y(40)) - (y(24)*y(73)*T(23));
    residual(33) = (y(71)*params(4)) - (params(43)*params(44)*y(70)^(params(44)-1));
    residual(34) = (params(5)*y(72)) - (y(71));
    residual(35) = (y(72)) - (y(24)*y(73)*T(24));
    residual(36) = (y(45)) - (y(44)*params(32));
    residual(37) = (y(46)) - (100*(log(y(21))-log(y(6))));
    residual(38) = (y(47)) - (100*(log(y(7))+log(y(27))-log(y(9))));
    residual(39) = (y(48)) - (100*(log(y(29))-log(y(11))));
    residual(40) = (y(49)) - (params(31)*y(16));
    residual(41) = (y(59)) - (y(7)*params(39)*(y(27)+y(32))/y(12)-1);
    residual(42) = (y(66)) - (100*(T(25)-(1+y(31))));
    residual(43) = (y(54)) - (y(73));
    residual(44) = (y(55)) - (y(73));
    residual(45) = (y(56)) - (T(26));
    residual(46) = (y(57)) - ((y(73)-y(55))^3);
    residual(47) = (y(58)) - (log(T(27)));
    residual(48) = (y(60)) - ((y(40)/y(14)-T(28))/T(28));
    residual(49) = (y(50)) - (y(73)/y(24)*(1+y(87)));
    residual(50) = (y(61)) - ((1+y(87))/y(50));
    residual(51) = (y(62)) - (100*(y(61)-(1+y(31))));
    residual(52) = (y(64)) - (y(84));
    residual(53) = (y(51)) - (log(T(29)));
    residual(54) = (y(52)) - (log(T(30)));
    residual(55) = (y(53)) - (log(T(31)));
    residual(56) = (y(63)) - (100*y(52));
    residual(57) = (y(68)) - (log(T(32)));
    residual(58) = (y(69)) - (log(1+y(31)));
residual(59) = y(67);
    residual(60) = (y(65)) - (log(y(71)));

end
