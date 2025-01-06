function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 56);

T(1) = 1-params(15);
T(2) = (y(4)*params(31))^params(1);
T(3) = params(40)*y(19)/T(2)-y(13);
T(4) = params(15)*y(29)^params(11);
T(5) = params(11)/(1-params(17));
T(6) = (params(7)*y(17))^T(5);
T(7) = T(1)*T(3)^params(11)+T(4)*T(6);
T(8) = 1/params(11);
T(9) = (y(49)*params(32))^(1-params(11));
T(10) = 1/T(2);
T(11) = T(9)*T(10);
T(12) = T(1)*T(3)^(params(11)-1)-y(48)*(1-params(12));
T(13) = params(15)*y(29)^(params(11)-1);
T(14) = (params(11)-(1-params(17)))/(1-params(17));
T(15) = (params(7)*y(17))^T(14);
T(16) = (params(32)*y(93))^(1-params(17)-params(11));
T(17) = T(13)*T(15)*T(16);
T(18) = (params(31)*y(21))^params(1);
T(19) = params(40)*y(81)/T(18)-y(47);
T(20) = T(1)*T(19)^(params(11)-1)+params(12)*y(92);
T(21) = params(15)*y(7)^(params(11)-1);
T(22) = 1/(1-params(17));
T(23) = (params(7)*y(1))^T(22);
T(24) = y(49)*params(32)/T(23);
T(25) = T(24)^(1-params(17)-params(11));
T(26) = T(21)*T(25);
T(27) = y(23)^params(14);
T(28) = params(22)^(1-params(14));
T(29) = 1-params(43)*y(75)^params(44);
T(30) = y(24)/y(3);
T(31) = T(30)^(-params(13));
T(32) = params(23)/(1-params(13));
T(33) = y(83)/y(20);
T(34) = T(33)^(1-params(13));
T(35) = y(87)/y(82);
T(36) = params(24)+T(32)*T(30)^(1-params(13));
T(37) = y(22)^params(28);
T(38) = params(21)^(1-params(28));
T(39) = params(5)*1000*(params(34)+params(1)*y(81)/(params(31)*y(21))+y(87)*params(25)*y(86)+y(88)*params(25)*y(85))+(1-params(27))*y(91);
T(40) = params(34)+params(1)*y(81)/(params(31)*y(21))+y(87)*params(25)*y(86)+y(88)*params(25)*y(85)+(1-params(27))*y(95);
T(41) = (y(89)+y(90))*y(29)*params(39)/y(37);
T(42) = (y(80)-y(60))^2;
T(43) = y(37)/y(32);
T(44) = (steady_state(29))/(steady_state(3));
T(45) = y(45)/(steady_state(29));
T(46) = y(19)/(steady_state(3));
T(47) = y(30)/(steady_state(14));
T(48) = T(41)/(1+y(36));
T(49) = y(47)/(steady_state(31));
T(50) = y(4)*params(31)/(params(31)*y(15));
T(51) = T(50)^(-params(1));
T(52) = T(26)*T(51);
T(53) = (params(31)*y(15))^params(1);
T(54) = params(40)*y(2)/T(53)-y(16);
T(55) = T(1)*T(54)^(params(11)-1)-(1-params(12))*y(14);
T(56) = T(12)/T(55);

end
