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

assert(length(T) >= 76);

T = call_CARA.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(34) = getPowerDeriv(T(3),1-params(17),1);
T(35) = getPowerDeriv(y(14),(-params(17)),1);
T(36) = 1/(steady_state(2));
T(37) = params(40)/T(9);
T(38) = getPowerDeriv(T(10),(-params(17)),1);
T(39) = params(1)/(params(31)*y(16));
T(40) = params(5)*1000*T(39);
T(41) = (-y(19))/(y(2)*y(2));
T(42) = getPowerDeriv(T(14),(-params(13)),1);
T(43) = getPowerDeriv(T(14),1-params(13),1);
T(44) = (-y(76))/(y(15)*y(15));
T(45) = getPowerDeriv(T(17),1-params(13),1);
T(46) = T(44)*T(45);
T(47) = params(31)*getPowerDeriv(y(3)*params(31),params(1),1);
T(48) = (-(params(40)*y(14)*T(47)))/(T(2)*T(2));
T(49) = (-(params(40)*T(47)))/(T(2)*T(2));
T(50) = getPowerDeriv(T(6),1-params(17),1);
T(51) = T(49)*T(50);
T(52) = params(31)*getPowerDeriv(params(31)*y(16),params(1),1);
T(53) = (-(params(40)*y(74)*T(52)))/(T(9)*T(9));
T(54) = (-(params(31)*params(1)*y(74)))/(params(31)*y(16)*params(31)*y(16));
T(55) = params(5)*1000*T(54);
T(56) = getPowerDeriv(y(17),params(28),1);
T(57) = (-y(80))/(y(75)*y(75));
T(58) = getPowerDeriv(y(18),params(14),1);
T(59) = 1/y(2);
T(60) = 1/y(15);
T(61) = T(16)*T(45)*T(60)-params(23)*T(45)*T(60);
T(62) = params(15)*getPowerDeriv(y(7),(-params(17)),1);
T(63) = params(15)*getPowerDeriv(y(24),1-params(17),1);
T(64) = params(39)*(y(82)+y(83))/y(32)/(1+y(31));
T(65) = 1/(steady_state(13));
T(66) = 1/y(75);
T(67) = 1/y(27);
T(68) = (-y(32))/(y(27)*y(27));
T(69) = y(24)*params(39)/y(32)/(1+y(31));
T(70) = (-(T(69)/T(32)));
T(71) = (-T(25))/((1+y(31))*(1+y(31)));
T(72) = (-((y(82)+y(83))*y(24)*params(39)))/(y(32)*y(32));
T(73) = T(72)/(1+y(31));
T(74) = 1/(steady_state(28));
T(75) = params(43)*getPowerDeriv(y(70),params(44),1);
T(76) = (-T(75));

end
