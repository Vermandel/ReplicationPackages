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

assert(length(T) >= 298);

T = call_EZ.dynamic_g2_tt(T, y, x, params, steady_state, it_);

T(262) = getPowerDeriv(T(7),T(8),3);
T(263) = T(144)*T(1)*T(145)+T(1)*(-T(73))*T(1)*(-T(73))*T(262);
T(264) = getPowerDeriv(T(3),params(11),3);
T(265) = (T(2)*T(2)*(-(params(40)*T(148)))-(-(params(40)*T(88)))*(T(2)*T(88)+T(2)*T(88)))/(T(2)*T(2)*T(2)*T(2));
T(266) = params(31)*params(31)*params(31)*getPowerDeriv(y(4)*params(31),params(1),3);
T(267) = (T(2)*T(2)*T(2)*T(2)*((-(params(40)*y(19)*T(148)))*(T(2)*T(88)+T(2)*T(88))+T(2)*T(2)*(-(params(40)*y(19)*T(266)))-((-(params(40)*y(19)*T(148)))*(T(2)*T(88)+T(2)*T(88))+(-(params(40)*y(19)*T(88)))*(T(88)*T(88)+T(2)*T(148)+T(88)*T(88)+T(2)*T(148))))-(T(2)*T(2)*(-(params(40)*y(19)*T(148)))-(-(params(40)*y(19)*T(88)))*(T(2)*T(88)+T(2)*T(88)))*(T(2)*T(2)*(T(2)*T(88)+T(2)*T(88))+T(2)*T(2)*(T(2)*T(88)+T(2)*T(88))))/(T(2)*T(2)*T(2)*T(2)*T(2)*T(2)*T(2)*T(2));
T(268) = getPowerDeriv(T(3),params(11)-1,3);
T(269) = getPowerDeriv(T(19),params(11)-1,3);
T(270) = getPowerDeriv(T(24),1-params(17)-params(11),3);
T(271) = getPowerDeriv(T(54),params(11)-1,3);
T(272) = T(139)*T(139)+T(55)*T(1)*(T(138)*T(138)*T(173)+T(68)*T(217))+T(139)*T(139)+T(55)*T(1)*(T(138)*T(138)*T(173)+T(68)*T(217));
T(273) = T(139)*T(1)*(-T(68))+T(55)*T(1)*T(138)*(-T(173))+T(139)*T(1)*(-T(68))+T(55)*T(1)*T(138)*(-T(173));
T(274) = T(1)*(-T(68))*T(1)*(-T(68))+T(55)*T(1)*T(173)+T(1)*(-T(68))*T(1)*(-T(68))+T(55)*T(1)*T(173);
T(275) = getPowerDeriv(T(50),(-params(1)),3);
T(276) = params(43)*getPowerDeriv(y(75),params(44),3);
T(277) = (y(3)*y(3)*y(3)*y(3)*(-(2*(-y(24))))-(-((-y(24))*(y(3)+y(3))))*(y(3)*y(3)*(y(3)+y(3))+y(3)*y(3)*(y(3)+y(3))))/(y(3)*y(3)*y(3)*y(3)*y(3)*y(3)*y(3)*y(3));
T(278) = getPowerDeriv(T(30),(-params(13)),3);
T(279) = getPowerDeriv(T(33),1-params(13),3);
T(280) = T(226)*T(85)*T(227)+T(86)*(y(20)*y(20)*y(20)*y(20)*(-(2*(-y(83))))-(-((-y(83))*(y(20)+y(20))))*(y(20)*y(20)*(y(20)+y(20))+y(20)*y(20)*(y(20)+y(20))))/(y(20)*y(20)*y(20)*y(20)*y(20)*y(20)*y(20)*y(20))+T(226)*T(85)*T(227)+T(85)*(T(226)*T(227)+T(85)*T(85)*T(279));
T(281) = T(226)*T(105)*T(227)+T(86)*(y(20)+y(20))/(y(20)*y(20)*y(20)*y(20))+T(85)*T(227)*T(229)+T(85)*(T(227)*T(229)+T(85)*T(105)*T(279));
T(282) = getPowerDeriv(T(30),1-params(13),3);
T(283) = (-(100*(-((-(y(7)*params(39)))*(y(12)+y(12))))/(y(12)*y(12)*y(12)*y(12))));
T(284) = (y(12)*y(12)*y(12)*y(12)*(-(2*(-(y(7)*params(39)*(y(32)+y(37))))))-(-((-(y(7)*params(39)*(y(32)+y(37))))*(y(12)+y(12))))*(y(12)*y(12)*(y(12)+y(12))+y(12)*y(12)*(y(12)+y(12))))/(y(12)*y(12)*y(12)*y(12)*y(12)*y(12)*y(12)*y(12));
T(285) = (-(100*(y(7)+y(7))/(y(7)*y(7)*y(7)*y(7))));
T(286) = (-((-(params(31)*params(1)))*(params(31)*params(31)*y(21)+params(31)*params(31)*y(21))))/(params(31)*y(21)*params(31)*y(21)*params(31)*y(21)*params(31)*y(21));
T(287) = (params(31)*y(21)*params(31)*y(21)*params(31)*y(21)*params(31)*y(21)*(-((-(params(31)*params(1)*y(81)))*(params(31)*params(31)+params(31)*params(31))))-(-((-(params(31)*params(1)*y(81)))*(params(31)*params(31)*y(21)+params(31)*params(31)*y(21))))*(params(31)*y(21)*params(31)*y(21)*(params(31)*params(31)*y(21)+params(31)*params(31)*y(21))+params(31)*y(21)*params(31)*y(21)*(params(31)*params(31)*y(21)+params(31)*params(31)*y(21))))/(params(31)*y(21)*params(31)*y(21)*params(31)*y(21)*params(31)*y(21)*params(31)*y(21)*params(31)*y(21)*params(31)*y(21)*params(31)*y(21));
T(288) = (-(100*(-((-(y(29)*params(39)))*(y(37)+y(37))))/(y(37)*y(37)*y(37)*y(37))));
T(289) = (y(37)*y(37)*y(37)*y(37)*(-(2*(-((y(89)+y(90))*y(29)*params(39)))))-(-((-((y(89)+y(90))*y(29)*params(39)))*(y(37)+y(37))))*(y(37)*y(37)*(y(37)+y(37))+y(37)*y(37)*(y(37)+y(37))))/(y(37)*y(37)*y(37)*y(37)*y(37)*y(37)*y(37)*y(37));
T(290) = T(115)*T(244);
T(291) = T(48)*T(48)*T(48)*T(48);
T(292) = (-((T(48)*T(48)*(-(T(111)*T(249)+T(111)*T(249)))-(-(T(111)*T(111)))*(T(48)*T(116)+T(48)*T(116)))/T(291)));
T(293) = (-((-((T(48)*T(249)-T(111)*T(116))*(T(48)*T(116)+T(48)*T(116))))/T(291)));
T(294) = T(48)*(-(params(39)/y(37)))/((1+y(36))*(1+y(36)));
T(295) = T(48)*(-params(39))/(y(37)*y(37))/(1+y(36));
T(296) = (-((-((-(T(116)*T(116)))*(T(48)*T(116)+T(48)*T(116))))/T(291)));
T(297) = T(48)*(-((1+y(36)+1+y(36))*(-(y(29)*params(39)/y(37)))))/T(235);
T(298) = T(48)*(-((-(y(29)*params(39)))*(y(37)+y(37))))/(y(37)*y(37)*y(37)*y(37))/(1+y(36));

end
