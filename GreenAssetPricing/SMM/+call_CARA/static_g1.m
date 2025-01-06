function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = call_CARA.static_g1_tt(T, y, x, params);
end
g1 = zeros(60, 60);
g1(1,2)=(-(T(2)*T(6)*T(27)));
g1(1,4)=(-(T(2)*T(27)*(-(params(40)*y(2)*T(32)))/(T(3)*T(3))));
g1(1,12)=(-(y(32)*params(15)*getPowerDeriv(y(12),1-params(17),1)));
g1(1,32)=1-T(1);
g1(2,1)=1;
g1(2,12)=(-T(34));
g1(3,2)=(-(T(7)*getPowerDeriv(y(2),(-params(17)),1)));
g1(3,4)=(-(T(8)*(-(params(40)*T(32)))/(T(3)*T(3))*getPowerDeriv(T(6),1-params(17),1)-(1-params(12))*y(31)*(-(params(40)*T(32)))/(T(3)*T(3))));
g1(3,8)=1;
g1(3,31)=T(6)*(1-params(12));
g1(4,2)=(-(T(5)*T(6)*T(28)));
g1(4,4)=(-(T(5)*T(28)*(-(params(40)*y(2)*T(32)))/(T(3)*T(3))));
g1(4,12)=(-(T(9)*T(34)));
g1(4,31)=1-T(5)*params(12);
g1(5,30)=1;
g1(6,4)=(-(T(11)*T(10)*y(16)*params(10)*exp((-params(25))*(y(4)*params(31)-params(3)))*params(31)*(-params(25))));
g1(6,6)=(-(T(11)*exp((-params(25))*(y(4)*params(31)-params(3)))*params(10)*y(16)*getPowerDeriv(y(6),params(14),1)));
g1(6,14)=1;
g1(6,16)=(-(T(11)*exp((-params(25))*(y(4)*params(31)-params(3)))*params(10)*T(10)));
g1(7,1)=(-((1+y(14)*(1-params(14)))*y(10)*params(39)+y(13)*y(11)*(1-params(28))));
g1(7,10)=1-y(1)*params(39)*(1+y(14)*(1-params(14)));
g1(7,11)=(-(y(13)*y(1)*(1-params(28))));
g1(7,13)=(-(y(1)*y(11)*(1-params(28))));
g1(7,14)=(-((1-params(14))*y(10)*y(1)*params(39)));
g1(8,12)=1;
g1(8,14)=(-1);
g1(9,5)=y(11)*params(28)*(-y(13))/(y(5)*y(5));
g1(9,6)=(-(params(14)*y(10)*(-y(14))/(y(6)*y(6))));
g1(9,10)=(-(params(14)*y(14)/y(6)));
g1(9,11)=params(28)*T(12);
g1(9,13)=y(11)*params(28)*1/y(5);
g1(9,14)=(-(params(14)*y(10)*1/y(6)));
g1(10,13)=(-((1-params(28))/params(21)));
g1(10,18)=1;
g1(11,7)=1;
g1(11,13)=(-T(13));
g1(11,15)=1;
g1(11,17)=y(59);
g1(11,18)=params(21);
g1(11,58)=(-(y(13)*T(36)));
g1(11,59)=y(17);
g1(12,11)=1;
g1(12,58)=(-(T(36)-params(4)*(-y(59))));
g1(12,59)=(1-y(58))*params(4);
g1(13,3)=(-(y(9)*params(23)*(-y(7))/(y(3)*y(3))*T(29)));
g1(13,7)=(-(y(9)*params(23)*T(29)*1/y(3)));
g1(13,9)=(-(params(23)*T(15)));
g1(14,1)=(-(y(11)*params(28)*T(12)+T(17)*params(39)*y(9)));
g1(14,3)=(-(y(1)*params(39)*y(9)*(params(23)/(1-params(13))*T(31)-params(23)*T(31))));
g1(14,5)=(-(y(1)*y(11)*params(28)*(-y(13))/(y(5)*y(5))));
g1(14,7)=(-(y(1)*params(39)*y(9)*(T(33)-params(23)*T(30)*1/y(3))));
g1(14,9)=1-y(1)*params(39)*T(17);
g1(14,11)=(-(T(12)*y(1)*params(28)));
g1(14,13)=(-(y(1)*y(11)*params(28)*1/y(5)));
g1(15,3)=y(12)-(1-params(19)+params(23)/(1-params(13))*T(16)+params(24)+y(3)*params(23)/(1-params(13))*T(31));
g1(15,7)=(-(y(3)*T(33)));
g1(15,12)=y(3);
g1(16,4)=(-(T(19)*T(18)*y(16)*params(6)*params(33)*exp((-params(25))*(y(4)*params(31)-params(3)))*params(31)*(-params(25))));
g1(16,5)=(-(T(19)*y(16)*exp((-params(25))*(y(4)*params(31)-params(3)))*params(6)*params(33)*getPowerDeriv(y(5),params(28),1)));
g1(16,13)=1;
g1(16,16)=(-(T(19)*exp((-params(25))*(y(4)*params(31)-params(3)))*params(6)*params(33)*T(18)));
g1(17,3)=1;
g1(17,5)=(-1);
g1(17,6)=(-1);
g1(18,4)=params(31)-params(31)*(1-params(27));
g1(18,17)=(-params(5));
g1(19,13)=(-(params(4)*(1-y(58))*params(41)));
g1(19,17)=1;
g1(19,58)=(-(y(13)*params(4)*(-params(41))));
g1(20,2)=(-1);
g1(20,7)=(-1);
g1(20,13)=T(13);
g1(20,58)=y(13)*T(36);
g1(21,16)=1/y(16)-params(18)*1/y(16);
g1(22,1)=(-((y(15)+y(20))*y(12)*params(39)));
g1(22,12)=(-(y(1)*params(39)*(y(15)+y(20))));
g1(22,15)=(-(params(39)*y(12)*y(1)));
g1(22,20)=1-params(39)*y(12)*y(1);
g1(23,1)=(-1);
g1(23,19)=(-1)/((1+y(19))*(1+y(19)));
g1(24,12)=(-(100*params(39)*(y(15)+y(20))/y(20)));
g1(24,15)=(-(100*y(12)*params(39)/y(20)));
g1(24,19)=100;
g1(24,20)=(-(100*T(35)));
g1(24,21)=1;
g1(25,19)=(-100);
g1(25,26)=1;
g1(26,2)=(-(100/y(13)));
g1(26,13)=(-((-(y(2)*100))/(y(13)*y(13))));
g1(26,25)=1;
g1(27,5)=(-((-(y(6)*100))/(y(5)*y(5))));
g1(27,6)=(-(100/y(5)));
g1(27,27)=1;
g1(28,12)=(-(100*1/y(12)));
g1(28,22)=1;
g1(29,12)=(-(1/y(12)));
g1(29,23)=1;
g1(30,12)=(-(100*1/y(12)));
g1(30,24)=1;
g1(31,2)=(-(params(1)*1000/(y(4)*params(31))));
g1(31,4)=(-((-(params(31)*y(2)*params(1)*1000))/(y(4)*params(31)*y(4)*params(31))));
g1(31,29)=1;
g1(32,1)=(-(y(12)*T(21)));
g1(32,2)=(-(y(12)*y(1)*params(5)*1000*params(1)/(y(4)*params(31))));
g1(32,4)=(-(y(12)*y(1)*params(5)*1000*(-(params(31)*y(2)*params(1)))/(y(4)*params(31)*y(4)*params(31))));
g1(32,10)=(-(y(12)*y(1)*params(5)*1000*y(14)*params(25)));
g1(32,11)=(-(y(12)*y(1)*params(5)*1000*params(25)*y(13)));
g1(32,12)=(-(y(1)*T(21)));
g1(32,13)=(-(y(12)*y(1)*params(5)*1000*params(25)*y(11)));
g1(32,14)=(-(y(12)*y(1)*params(5)*1000*params(25)*y(10)));
g1(32,28)=1-(1-params(27))*y(12)*y(1);
g1(33,58)=(-(params(43)*params(44)*getPowerDeriv(y(58),params(44)-1,1)));
g1(33,59)=params(4);
g1(34,59)=(-1);
g1(34,60)=params(5);
g1(35,1)=(-(y(12)*T(22)));
g1(35,2)=(-(y(12)*y(1)*params(1)/(y(4)*params(31))));
g1(35,4)=(-(y(12)*y(1)*(-(params(31)*y(2)*params(1)))/(y(4)*params(31)*y(4)*params(31))));
g1(35,10)=(-(y(12)*y(1)*y(14)*params(25)));
g1(35,11)=(-(y(12)*y(1)*params(25)*y(13)));
g1(35,12)=(-(y(1)*T(22)));
g1(35,13)=(-(y(12)*y(1)*params(25)*y(11)));
g1(35,14)=(-(y(12)*y(1)*params(25)*y(10)));
g1(35,60)=1-(1-params(27))*y(12)*y(1);
g1(36,32)=(-params(32));
g1(36,33)=1;
g1(37,34)=1;
g1(38,12)=(-(100*1/y(12)));
g1(38,35)=1;
g1(39,36)=1;
g1(40,4)=(-params(31));
g1(40,37)=1;
g1(41,12)=(-(params(39)*(y(15)+y(20))/y(20)));
g1(41,15)=(-(y(12)*params(39)/y(20)));
g1(41,20)=(-T(35));
g1(41,47)=1;
g1(42,12)=(-(100*params(39)*(y(15)+y(20))/y(20)));
g1(42,15)=(-(100*y(12)*params(39)/y(20)));
g1(42,19)=100;
g1(42,20)=(-(100*T(35)));
g1(42,54)=1;
g1(43,1)=(-1);
g1(43,42)=1;
g1(44,1)=(-1);
g1(44,43)=1;
g1(45,1)=(-(2*(y(1)-y(43))));
g1(45,43)=2*(y(1)-y(43));
g1(45,44)=1;
g1(46,1)=(-(3*T(23)));
g1(46,43)=3*T(23);
g1(46,45)=1;
g1(47,15)=(-((-y(20))/(y(15)*y(15))/(y(20)/y(15))));
g1(47,20)=(-(1/y(15)/(y(20)/y(15))));
g1(47,46)=1;
g1(48,2)=(-((T(24)*((-y(28))/(y(2)*y(2))-(-(y(28)))/((y(2))*(y(2))))-T(25)*(-(y(28)))/((y(2))*(y(2))))/(T(24)*T(24))));
g1(48,28)=(-((T(24)*(1/y(2)-1/(y(2)))-T(25)*1/(y(2)))/(T(24)*T(24))));
g1(48,48)=1;
g1(49,1)=(-((1+y(38))*1/y(12)));
g1(49,12)=(-((1+y(38))*(-y(1))/(y(12)*y(12))));
g1(49,38)=1-y(1)/y(12);
g1(50,38)=(-((y(38)-(1+y(38)))/(y(38)*y(38))));
g1(50,49)=1;
g1(51,19)=100;
g1(51,49)=(-100);
g1(51,50)=1;
g1(52,28)=(-1);
g1(52,52)=1;
g1(53,28)=(-(((y(28))-y(28))/((y(28))*(y(28)))/(y(28)/(y(28)))));
g1(53,39)=1;
g1(54,2)=(-(((y(2))-y(2))/((y(2))*(y(2)))/(y(2)/(y(2)))));
g1(54,40)=1;
g1(55,13)=(-(((y(13))-y(13))/((y(13))*(y(13)))/(y(13)/(y(13)))));
g1(55,41)=1;
g1(56,40)=(-100);
g1(56,51)=1;
g1(57,12)=(-(params(39)*(y(15)+y(20))/y(20)/(1+y(19))/T(26)));
g1(57,15)=(-(y(12)*params(39)/y(20)/(1+y(19))/T(26)));
g1(57,19)=(-((-T(20))/((1+y(19))*(1+y(19)))/T(26)));
g1(57,20)=(-(T(35)/(1+y(19))/T(26)));
g1(57,56)=1;
g1(58,19)=(-(1/(1+y(19))));
g1(58,57)=1;
g1(59,55)=1;
g1(60,53)=1;
g1(60,59)=(-(1/y(59)));

end
