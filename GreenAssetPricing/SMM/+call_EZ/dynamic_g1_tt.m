function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 141);

T = call_EZ.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(57) = params(7)*getPowerDeriv(params(7)*y(1),T(22),1);
T(58) = (-(y(49)*params(32)*T(57)));
T(59) = T(58)/(T(23)*T(23));
T(60) = getPowerDeriv(T(24),1-params(17)-params(11),1);
T(61) = T(51)*T(21)*T(59)*T(60);
T(62) = params(7)*getPowerDeriv(params(7)*y(17),T(5),1);
T(63) = T(4)*T(62);
T(64) = getPowerDeriv(T(7),T(8),1);
T(65) = params(7)*getPowerDeriv(params(7)*y(17),T(14),1);
T(66) = T(16)*T(13)*T(65);
T(67) = params(40)/T(53);
T(68) = getPowerDeriv(T(54),params(11)-1,1);
T(69) = T(1)*T(67)*T(68);
T(70) = T(55)*T(55);
T(71) = (-(T(12)*T(69)))/T(70);
T(72) = params(40)/T(2);
T(73) = getPowerDeriv(T(3),params(11),1);
T(74) = T(1)*T(72)*T(73);
T(75) = getPowerDeriv(T(3),params(11)-1,1);
T(76) = T(1)*T(72)*T(75)/T(55);
T(77) = 1/(steady_state(3));
T(78) = params(40)/T(18);
T(79) = getPowerDeriv(T(19),params(11)-1,1);
T(80) = params(1)/(params(31)*y(21));
T(81) = params(5)*1000*T(80);
T(82) = (-y(24))/(y(3)*y(3));
T(83) = getPowerDeriv(T(30),(-params(13)),1);
T(84) = getPowerDeriv(T(30),1-params(13),1);
T(85) = (-y(83))/(y(20)*y(20));
T(86) = getPowerDeriv(T(33),1-params(13),1);
T(87) = T(85)*T(86);
T(88) = params(31)*getPowerDeriv(y(4)*params(31),params(1),1);
T(89) = (-(params(40)*y(19)*T(88)))/(T(2)*T(2));
T(90) = T(9)*(-T(88))/(T(2)*T(2));
T(91) = params(31)/(params(31)*y(15));
T(92) = getPowerDeriv(T(50),(-params(1)),1);
T(93) = T(91)*T(92);
T(94) = T(26)*T(93);
T(95) = T(1)*T(75)*T(89)/T(55);
T(96) = params(31)*getPowerDeriv(params(31)*y(21),params(1),1);
T(97) = (-(params(40)*y(81)*T(96)))/(T(18)*T(18));
T(98) = T(1)*T(79)*T(97);
T(99) = (-(params(31)*params(1)*y(81)))/(params(31)*y(21)*params(31)*y(21));
T(100) = params(5)*1000*T(99);
T(101) = getPowerDeriv(y(22),params(28),1);
T(102) = (-y(87))/(y(82)*y(82));
T(103) = getPowerDeriv(y(23),params(14),1);
T(104) = 1/y(3);
T(105) = 1/y(20);
T(106) = T(32)*T(86)*T(105)-params(23)*T(86)*T(105);
T(107) = params(15)*getPowerDeriv(y(7),params(11)-1,1);
T(108) = T(51)*T(25)*T(107);
T(109) = params(15)*getPowerDeriv(y(29),params(11),1);
T(110) = params(15)*getPowerDeriv(y(29),params(11)-1,1);
T(111) = params(39)*(y(89)+y(90))/y(37)/(1+y(36));
T(112) = 1/(steady_state(14));
T(113) = 1/y(82);
T(114) = 1/y(32);
T(115) = (-y(37))/(y(32)*y(32));
T(116) = y(29)*params(39)/y(37)/(1+y(36));
T(117) = (-(T(116)/T(48)));
T(118) = (-T(41))/((1+y(36))*(1+y(36)));
T(119) = (-((y(89)+y(90))*y(29)*params(39)))/(y(37)*y(37));
T(120) = T(119)/(1+y(36));
T(121) = 1/(steady_state(29));
T(122) = T(1)*(-T(75))/T(55);
T(123) = 1/(steady_state(31));
T(124) = (-(1-params(12)));
T(125) = (-(T(12)*T(124)))/T(70);
T(126) = T(124)/T(55);
T(127) = params(32)*getPowerDeriv(y(49)*params(32),1-params(11),1);
T(128) = T(10)*T(127);
T(129) = params(32)/T(23);
T(130) = T(51)*T(21)*T(60)*T(129);
T(131) = params(32)*getPowerDeriv(params(32)*y(93),1-params(17)-params(11),1);
T(132) = params(43)*getPowerDeriv(y(75),params(44),1);
T(133) = (-T(132));
T(134) = (-(params(31)*y(4)*params(31)))/(params(31)*y(15)*params(31)*y(15));
T(135) = T(92)*T(134);
T(136) = T(26)*T(135);
T(137) = params(31)*getPowerDeriv(params(31)*y(15),params(1),1);
T(138) = (-(params(40)*y(2)*T(137)))/(T(53)*T(53));
T(139) = T(1)*T(68)*T(138);
T(140) = (-(T(12)*T(139)))/T(70);
T(141) = (-(T(12)*T(1)*(-T(68))))/T(70);

end
