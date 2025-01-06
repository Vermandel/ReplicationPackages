function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 26
    T = [T; NaN(26 - size(T, 1), 1)];
end
T(1) = params(15)*y(12)^(1-params(17));
T(2) = (1-params(17))^(-1);
T(3) = (y(4)*params(31))^params(1);
T(4) = params(40)*y(2)/T(3);
T(5) = params(15)*y(12)^(-params(17));
T(6) = params(40)/T(3);
T(7) = T(6)^(1-params(17));
T(8) = y(2)^(-params(17));
T(9) = y(31)*params(12)+T(4)^(-params(17));
T(10) = y(6)^params(14);
T(11) = params(22)^(1-params(14));
T(12) = y(13)/y(5);
T(13) = 1-params(43)*y(58)^params(44);
T(14) = y(7)/y(3);
T(15) = T(14)^(-params(13));
T(16) = T(14)^(1-params(13));
T(17) = 1-params(19)+params(23)/(1-params(13))*T(16)-params(23)*T(16)+params(24);
T(18) = y(5)^params(28);
T(19) = params(21)^(1-params(28));
T(20) = (y(15)+y(20))*y(12)*params(39)/y(20);
T(21) = params(5)*1000*(params(34)+y(2)*params(1)/(y(4)*params(31))+y(13)*params(25)*y(11)+y(14)*params(25)*y(10))+(1-params(27))*y(28);
T(22) = params(34)+y(2)*params(1)/(y(4)*params(31))+y(13)*params(25)*y(11)+y(14)*params(25)*y(10)+(1-params(27))*y(60);
T(23) = (y(1)-y(43))^2;
T(24) = (y(28))/(y(2));
T(25) = y(28)/y(2)-T(24);
T(26) = T(20)/(1+y(19));
end
