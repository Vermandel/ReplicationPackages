function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 130);

T = call_CARA.dynamic_g1_tt(T, y, x, params, steady_state, it_);

T(77) = getPowerDeriv(T(3),1-params(17),2);
T(78) = params(31)*params(31)*getPowerDeriv(y(3)*params(31),params(1),2);
T(79) = (T(2)*T(2)*(-(params(40)*y(14)*T(78)))-(-(params(40)*y(14)*T(47)))*(T(2)*T(47)+T(2)*T(47)))/(T(2)*T(2)*T(2)*T(2));
T(80) = params(15)*getPowerDeriv(y(24),1-params(17),2);
T(81) = params(15)*getPowerDeriv(y(7),(-params(17)),2);
T(82) = getPowerDeriv(y(14),(-params(17)),2);
T(83) = (T(2)*T(2)*(-(params(40)*T(78)))-(-(params(40)*T(47)))*(T(2)*T(47)+T(2)*T(47)))/(T(2)*T(2)*T(2)*T(2));
T(84) = getPowerDeriv(T(6),1-params(17),2);
T(85) = T(49)*T(84);
T(86) = T(50)*T(83)+T(49)*T(85);
T(87) = getPowerDeriv(T(10),(-params(17)),2);
T(88) = (-(params(40)*T(52)))/(T(9)*T(9));
T(89) = params(31)*params(31)*getPowerDeriv(params(31)*y(16),params(1),2);
T(90) = (T(9)*T(9)*(-(params(40)*y(74)*T(89)))-(-(params(40)*y(74)*T(52)))*(T(9)*T(52)+T(9)*T(52)))/(T(9)*T(9)*T(9)*T(9));
T(91) = getPowerDeriv(y(18),params(14),2);
T(92) = params(43)*getPowerDeriv(y(70),params(44),2);
T(93) = (-((-y(19))*(y(2)+y(2))))/(y(2)*y(2)*y(2)*y(2));
T(94) = getPowerDeriv(T(14),(-params(13)),2);
T(95) = (-1)/(y(2)*y(2));
T(96) = (-((-y(76))*(y(15)+y(15))))/(y(15)*y(15)*y(15)*y(15));
T(97) = getPowerDeriv(T(17),1-params(13),2);
T(98) = T(45)*T(96)+T(44)*T(44)*T(97);
T(99) = (-1)/(y(15)*y(15));
T(100) = T(45)*T(99)+T(44)*T(60)*T(97);
T(101) = (-((-y(80))*(y(75)+y(75))))/(y(75)*y(75)*y(75)*y(75));
T(102) = T(60)*T(60)*T(97);
T(103) = getPowerDeriv(T(14),1-params(13),2);
T(104) = getPowerDeriv(y(17),params(28),2);
T(105) = (1+y(31))*(1+y(31))*(1+y(31))*(1+y(31));
T(106) = (-((-(y(7)*params(39)*(y(27)+y(32))))*(y(12)+y(12))))/(y(12)*y(12)*y(12)*y(12));
T(107) = (-(100*(-1)/(y(7)*y(7))));
T(108) = (-(params(31)*params(1)))/(params(31)*y(16)*params(31)*y(16));
T(109) = params(5)*1000*T(108);
T(110) = (-((-(params(31)*params(1)*y(74)))*(params(31)*params(31)*y(16)+params(31)*params(31)*y(16))))/(params(31)*y(16)*params(31)*y(16)*params(31)*y(16)*params(31)*y(16));
T(111) = (-(params(39)*(y(82)+y(83))))/(y(32)*y(32));
T(112) = (-(y(24)*params(39)))/(y(32)*y(32));
T(113) = (-((-((y(82)+y(83))*y(24)*params(39)))*(y(32)+y(32))))/(y(32)*y(32)*y(32)*y(32));
T(114) = (-((-y(32))*(y(27)+y(27))))/(y(27)*y(27)*y(27)*y(27));
T(115) = T(27)*T(114)-T(68)*T(68);
T(116) = (-(T(74)*T(74)));
T(117) = (-(T(36)*T(36)));
T(118) = (-(T(65)*T(65)));
T(119) = params(39)/y(32)/(1+y(31));
T(120) = (-((T(32)*T(119)-T(64)*T(69))/(T(32)*T(32))));
T(121) = (-(params(39)*(y(82)+y(83))/y(32)))/((1+y(31))*(1+y(31)));
T(122) = T(111)/(1+y(31));
T(123) = (-((-(T(69)*T(69)))/(T(32)*T(32))));
T(124) = (-(y(24)*params(39)/y(32)))/((1+y(31))*(1+y(31)));
T(125) = (-((T(32)*T(124)-T(69)*T(71))/(T(32)*T(32))));
T(126) = T(112)/(1+y(31));
T(127) = (-((T(32)*T(126)-T(69)*T(73))/(T(32)*T(32))));
T(128) = (-((-T(25))*(1+y(31)+1+y(31))))/T(105);
T(129) = (-T(72))/((1+y(31))*(1+y(31)));
T(130) = T(113)/(1+y(31));

end
