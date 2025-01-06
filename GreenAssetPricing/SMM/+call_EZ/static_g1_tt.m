function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 53);

T = call_EZ.static_resid_tt(T, y, x, params);

T(38) = getPowerDeriv(T(5),1/params(11),1);
T(39) = getPowerDeriv(y(33)*params(32)/T(13),1-params(17)-params(11),1);
T(40) = getPowerDeriv(T(2),params(11),1);
T(41) = getPowerDeriv(T(2),params(11)-1,1);
T(42) = (1-params(15))*params(40)/T(1)*T(41);
T(43) = getPowerDeriv(params(40)*y(3)/T(18)-y(63),params(11)-1,1);
T(44) = getPowerDeriv(T(25),(-params(13)),1);
T(45) = getPowerDeriv(T(25),1-params(13),1);
T(46) = (-y(8))/(y(4)*y(4))*T(45);
T(47) = params(31)*getPowerDeriv(y(5)*params(31),params(1),1);
T(48) = (1-params(15))*T(41)*(-(params(40)*y(3)*T(47)))/(T(1)*T(1));
T(49) = getPowerDeriv(y(5)*params(31)/(params(31)*y(62)),(-params(1)),1);
T(50) = params(23)/(1-params(13))*T(45)*1/y(4);
T(51) = params(15)*getPowerDeriv(y(13),params(11)-1,1);
T(52) = (y(21)*y(13)*params(39)-(y(16)+y(21))*y(13)*params(39))/(y(21)*y(21));
T(53) = (-(params(43)*getPowerDeriv(y(59),params(44),1)));

end
