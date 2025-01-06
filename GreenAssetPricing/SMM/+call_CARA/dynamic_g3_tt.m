function T = dynamic_g3_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g3_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 158);

T = call_CARA.dynamic_g2_tt(T, y, x, params, steady_state, it_);

T(131) = getPowerDeriv(T(3),1-params(17),3);
T(132) = params(31)*params(31)*params(31)*getPowerDeriv(y(3)*params(31),params(1),3);
T(133) = params(15)*getPowerDeriv(y(7),(-params(17)),3);
T(134) = (T(2)*T(2)*T(2)*T(2)*((T(2)*T(47)+T(2)*T(47))*(-(params(40)*T(78)))+T(2)*T(2)*(-(params(40)*T(132)))-((T(2)*T(47)+T(2)*T(47))*(-(params(40)*T(78)))+(-(params(40)*T(47)))*(T(47)*T(47)+T(2)*T(78)+T(47)*T(47)+T(2)*T(78))))-(T(2)*T(2)*(-(params(40)*T(78)))-(-(params(40)*T(47)))*(T(2)*T(47)+T(2)*T(47)))*(T(2)*T(2)*(T(2)*T(47)+T(2)*T(47))+T(2)*T(2)*(T(2)*T(47)+T(2)*T(47))))/(T(2)*T(2)*T(2)*T(2)*T(2)*T(2)*T(2)*T(2));
T(135) = getPowerDeriv(T(10),(-params(17)),3);
T(136) = params(43)*getPowerDeriv(y(70),params(44),3);
T(137) = (y(2)*y(2)*y(2)*y(2)*(-(2*(-y(19))))-(-((-y(19))*(y(2)+y(2))))*(y(2)*y(2)*(y(2)+y(2))+y(2)*y(2)*(y(2)+y(2))))/(y(2)*y(2)*y(2)*y(2)*y(2)*y(2)*y(2)*y(2));
T(138) = getPowerDeriv(T(14),(-params(13)),3);
T(139) = getPowerDeriv(T(17),1-params(13),3);
T(140) = T(96)*T(44)*T(97)+T(45)*(y(15)*y(15)*y(15)*y(15)*(-(2*(-y(76))))-(-((-y(76))*(y(15)+y(15))))*(y(15)*y(15)*(y(15)+y(15))+y(15)*y(15)*(y(15)+y(15))))/(y(15)*y(15)*y(15)*y(15)*y(15)*y(15)*y(15)*y(15))+T(96)*T(44)*T(97)+T(44)*(T(96)*T(97)+T(44)*T(44)*T(139));
T(141) = T(96)*T(60)*T(97)+T(45)*(y(15)+y(15))/(y(15)*y(15)*y(15)*y(15))+T(44)*T(97)*T(99)+T(44)*(T(97)*T(99)+T(44)*T(60)*T(139));
T(142) = getPowerDeriv(T(14),1-params(13),3);
T(143) = (-(100*(-((-(y(7)*params(39)))*(y(12)+y(12))))/(y(12)*y(12)*y(12)*y(12))));
T(144) = (y(12)*y(12)*y(12)*y(12)*(-(2*(-(y(7)*params(39)*(y(27)+y(32))))))-(-((-(y(7)*params(39)*(y(27)+y(32))))*(y(12)+y(12))))*(y(12)*y(12)*(y(12)+y(12))+y(12)*y(12)*(y(12)+y(12))))/(y(12)*y(12)*y(12)*y(12)*y(12)*y(12)*y(12)*y(12));
T(145) = (-(100*(y(7)+y(7))/(y(7)*y(7)*y(7)*y(7))));
T(146) = (-((-(params(31)*params(1)))*(params(31)*params(31)*y(16)+params(31)*params(31)*y(16))))/(params(31)*y(16)*params(31)*y(16)*params(31)*y(16)*params(31)*y(16));
T(147) = (params(31)*y(16)*params(31)*y(16)*params(31)*y(16)*params(31)*y(16)*(-((-(params(31)*params(1)*y(74)))*(params(31)*params(31)+params(31)*params(31))))-(-((-(params(31)*params(1)*y(74)))*(params(31)*params(31)*y(16)+params(31)*params(31)*y(16))))*(params(31)*y(16)*params(31)*y(16)*(params(31)*params(31)*y(16)+params(31)*params(31)*y(16))+params(31)*y(16)*params(31)*y(16)*(params(31)*params(31)*y(16)+params(31)*params(31)*y(16))))/(params(31)*y(16)*params(31)*y(16)*params(31)*y(16)*params(31)*y(16)*params(31)*y(16)*params(31)*y(16)*params(31)*y(16)*params(31)*y(16));
T(148) = (-(100*(-((-(y(24)*params(39)))*(y(32)+y(32))))/(y(32)*y(32)*y(32)*y(32))));
T(149) = (y(32)*y(32)*y(32)*y(32)*(-(2*(-((y(82)+y(83))*y(24)*params(39)))))-(-((-((y(82)+y(83))*y(24)*params(39)))*(y(32)+y(32))))*(y(32)*y(32)*(y(32)+y(32))+y(32)*y(32)*(y(32)+y(32))))/(y(32)*y(32)*y(32)*y(32)*y(32)*y(32)*y(32)*y(32));
T(150) = T(68)*T(114);
T(151) = T(32)*T(32)*T(32)*T(32);
T(152) = (-((T(32)*T(32)*(-(T(64)*T(119)+T(64)*T(119)))-(-(T(64)*T(64)))*(T(32)*T(69)+T(32)*T(69)))/T(151)));
T(153) = (-((-((T(32)*T(119)-T(64)*T(69))*(T(32)*T(69)+T(32)*T(69))))/T(151)));
T(154) = T(32)*(-(params(39)/y(32)))/((1+y(31))*(1+y(31)));
T(155) = T(32)*(-params(39))/(y(32)*y(32))/(1+y(31));
T(156) = (-((-((-(T(69)*T(69)))*(T(32)*T(69)+T(32)*T(69))))/T(151)));
T(157) = T(32)*(-((1+y(31)+1+y(31))*(-(y(24)*params(39)/y(32)))))/T(105);
T(158) = T(32)*(-((-(y(24)*params(39)))*(y(32)+y(32))))/(y(32)*y(32)*y(32)*y(32))/(1+y(31));

end
