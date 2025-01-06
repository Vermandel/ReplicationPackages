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

assert(length(T) >= 33);

T(1) = (1-params(17))^(-1);
T(2) = (y(3)*params(31))^params(1);
T(3) = params(40)*y(14)/T(2);
T(4) = params(15)*y(24)^(1-params(17));
T(5) = params(15)*y(7)^(-params(17));
T(6) = params(40)/T(2);
T(7) = T(6)^(1-params(17));
T(8) = y(14)^(-params(17));
T(9) = (params(31)*y(16))^params(1);
T(10) = params(40)*y(74)/T(9);
T(11) = y(18)^params(14);
T(12) = params(22)^(1-params(14));
T(13) = 1-params(43)*y(70)^params(44);
T(14) = y(19)/y(2);
T(15) = T(14)^(-params(13));
T(16) = params(23)/(1-params(13));
T(17) = y(76)/y(15);
T(18) = T(17)^(1-params(13));
T(19) = y(80)/y(75);
T(20) = params(24)+T(16)*T(14)^(1-params(13));
T(21) = y(17)^params(28);
T(22) = params(21)^(1-params(28));
T(23) = params(5)*1000*(params(34)+params(1)*y(74)/(params(31)*y(16))+y(80)*params(25)*y(79)+y(81)*params(25)*y(78))+(1-params(27))*y(84);
T(24) = params(34)+params(1)*y(74)/(params(31)*y(16))+y(80)*params(25)*y(79)+y(81)*params(25)*y(78)+(1-params(27))*y(88);
T(25) = (y(82)+y(83))*y(24)*params(39)/y(32);
T(26) = (y(73)-y(55))^2;
T(27) = y(32)/y(27);
T(28) = (steady_state(28))/(steady_state(2));
T(29) = y(40)/(steady_state(28));
T(30) = y(14)/(steady_state(2));
T(31) = y(25)/(steady_state(13));
T(32) = T(25)/(1+y(31));
T(33) = params(12)*y(85)+T(10)^(-params(17));

end
