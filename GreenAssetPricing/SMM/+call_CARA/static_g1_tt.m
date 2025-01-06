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

assert(length(T) >= 36);

T = call_CARA.static_resid_tt(T, y, x, params);

T(27) = getPowerDeriv(T(4),1-params(17),1);
T(28) = getPowerDeriv(T(4),(-params(17)),1);
T(29) = getPowerDeriv(T(14),(-params(13)),1);
T(30) = getPowerDeriv(T(14),1-params(13),1);
T(31) = (-y(7))/(y(3)*y(3))*T(30);
T(32) = params(31)*getPowerDeriv(y(4)*params(31),params(1),1);
T(33) = params(23)/(1-params(13))*T(30)*1/y(3);
T(34) = params(15)*getPowerDeriv(y(12),(-params(17)),1);
T(35) = (y(20)*y(12)*params(39)-(y(15)+y(20))*y(12)*params(39))/(y(20)*y(20));
T(36) = (-(params(43)*getPowerDeriv(y(58),params(44),1)));

end
