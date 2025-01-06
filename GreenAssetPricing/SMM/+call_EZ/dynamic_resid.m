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
    T = call_EZ.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(63, 1);
    residual(1) = (y(49)*params(32)) - (T(7)^T(8));
    residual(2) = (params(7)*y(17)) - ((params(32)*y(93))^(1-params(17)));
    residual(3) = (y(25)) - (T(11)*T(12));
    residual(4) = (y(48)) - (T(17)*T(20));
    residual(5) = (y(18)) - (T(52)*T(56));
    residual(6) = (y(29)*y(47)) - (y(13)*params(12)+y(19)*params(40)*(1-params(12))/T(2));
    residual(7) = (y(31)) - (exp((-params(25))*(y(4)*params(31)-params(3)))*params(10)*y(33)*T(27)*T(28));
    residual(8) = (y(27)) - (params(39)*y(80)*y(85)*(1+(1-params(14))*y(88))+y(80)*y(86)*(1-params(28))*y(87));
    residual(9) = (y(29)) - (1+y(31));
    residual(10) = (params(28)*y(28)*y(30)/y(22)) - (params(14)*y(27)*y(31)/y(23));
    residual(11) = (y(35)) - ((1-params(28))*y(30)/params(21));
    residual(12) = (y(32)) - (y(30)*T(29)-y(35)*params(21)-y(24)-y(76)*y(34));
    residual(13) = (y(28)) - (T(29)-y(76)*(1-y(75))*params(4));
    residual(14) = (1) - (y(26)*params(23)*T(31));
    residual(15) = (y(26)) - (params(39)*y(80)*y(84)*(1-params(19)+T(32)*T(34)-params(23)*T(34)+params(24))+y(80)*y(86)*params(28)*T(35));
    residual(16) = (y(29)*y(20)) - (y(3)*(1-params(19))+y(3)*T(36));
    residual(17) = (y(30)) - (y(33)*exp((-params(25))*(y(4)*params(31)-params(3)))*params(6)*params(33)*T(37)*T(38));
    residual(18) = (y(3)) - (y(23)+y(22));
    residual(19) = (params(31)*y(21)) - (params(31)*y(4)*(1-params(27))+params(5)*(y(34)+params(9)));
    residual(20) = (y(34)) - (y(30)*params(4)*(1-y(75))*params(41));
    residual(21) = (y(30)*T(29)) - (y(19)+y(24)+params(29)*params(42));
    residual(22) = (log(y(33))) - (params(18)*log(y(10))+params(30)*x(it_, 1));
    residual(23) = (y(37)) - (params(39)*y(29)*y(80)*(y(89)+y(90)));
    residual(24) = (1/(1+y(36))) - (y(80));
    residual(25) = (y(38)) - ((y(7)*params(39)*(y(32)+y(37))/y(12)-(1+y(36)))*100);
    residual(26) = (y(43)) - (y(36)*100);
    residual(27) = (y(42)) - (y(19)*100/y(30));
    residual(28) = (y(44)) - (y(23)*100/y(22));
    residual(29) = (y(39)) - (100*(log(y(7))+log(y(30))-log(y(8))));
    residual(30) = (y(40)) - (log(y(7))+log(y(19))-log(y(2)));
    residual(31) = (y(41)) - (100*(log(y(7))+log(y(24))-log(y(5))));
    residual(32) = (y(46)) - (y(19)*params(1)*1000/(y(4)*params(31)));
    residual(33) = (y(45)) - (y(29)*y(80)*T(39));
    residual(34) = (y(76)*params(4)) - (params(43)*params(44)*y(75)^(params(44)-1));
    residual(35) = (params(5)*y(77)) - (y(76));
    residual(36) = (y(77)) - (y(29)*y(80)*T(40));
    residual(37) = (y(50)) - (y(49)*params(32));
    residual(38) = (y(51)) - (100*(log(y(26))-log(y(6))));
    residual(39) = (y(52)) - (100*(log(y(7))+log(y(32))-log(y(9))));
    residual(40) = (y(53)) - (100*(log(y(34))-log(y(11))));
    residual(41) = (y(54)) - (params(31)*y(21));
    residual(42) = (y(64)) - (y(7)*params(39)*(y(32)+y(37))/y(12)-1);
    residual(43) = (y(71)) - (100*(T(41)-(1+y(36))));
    residual(44) = (y(59)) - (y(80));
    residual(45) = (y(60)) - (y(80));
    residual(46) = (y(61)) - (T(42));
    residual(47) = (y(62)) - ((y(80)-y(60))^3);
    residual(48) = (y(63)) - (log(T(43)));
    residual(49) = (y(65)) - ((y(45)/y(19)-T(44))/T(44));
    residual(50) = (y(55)) - (y(80)/y(29)*(1+y(94)));
    residual(51) = (y(66)) - ((1+y(94))/y(55));
    residual(52) = (y(67)) - (100*(y(66)-(1+y(36))));
    residual(53) = (y(69)) - (y(91));
    residual(54) = (y(56)) - (log(T(45)));
    residual(55) = (y(57)) - (log(T(46)));
    residual(56) = (y(58)) - (log(T(47)));
    residual(57) = (y(68)) - (100*y(57));
    residual(58) = (y(73)) - (log(T(48)));
    residual(59) = (y(74)) - (log(1+y(36)));
    residual(60) = (y(72)) - (log(T(49)));
    residual(61) = (y(70)) - (log(y(76)));
    residual(62) = (y(78)) - (y(4));
    residual(63) = (y(79)) - (y(13));

end
