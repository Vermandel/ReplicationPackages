function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 37
    T = [T; NaN(37 - size(T, 1), 1)];
end
T(1) = (y(5)*params(31))^params(1);
T(2) = params(40)*y(3)/T(1)-y(31);
T(3) = params(15)*y(13)^params(11);
T(4) = (params(7)*y(1))^(params(11)/(1-params(17)));
T(5) = (1-params(15))*T(2)^params(11)+T(3)*T(4);
T(6) = (y(33)*params(32))^(1-params(11));
T(7) = T(6)*1/T(1);
T(8) = (1-params(15))*T(2)^(params(11)-1);
T(9) = params(15)*y(13)^(params(11)-1);
T(10) = (params(7)*y(1))^((params(11)-(1-params(17)))/(1-params(17)));
T(11) = (y(33)*params(32))^(1-params(17)-params(11));
T(12) = T(9)*T(10)*T(11);
T(13) = (params(7)*y(1))^(1/(1-params(17)));
T(14) = (y(33)*params(32)/T(13))^(1-params(17)-params(11));
T(15) = T(9)*T(14);
T(16) = (y(5)*params(31)/(params(31)*y(62)))^(-params(1));
T(17) = T(15)*T(16);
T(18) = (params(31)*y(62))^params(1);
T(19) = (1-params(15))*(params(40)*y(3)/T(18)-y(63))^(params(11)-1)-y(32)*(1-params(12));
T(20) = (T(8)-y(32)*(1-params(12)))/T(19);
T(21) = y(7)^params(14);
T(22) = params(22)^(1-params(14));
T(23) = y(14)/y(6);
T(24) = 1-params(43)*y(59)^params(44);
T(25) = y(8)/y(4);
T(26) = T(25)^(-params(13));
T(27) = T(25)^(1-params(13));
T(28) = 1-params(19)+params(23)/(1-params(13))*T(27)-params(23)*T(27)+params(24);
T(29) = y(6)^params(28);
T(30) = params(21)^(1-params(28));
T(31) = (y(16)+y(21))*y(13)*params(39)/y(21);
T(32) = params(5)*1000*(params(34)+y(3)*params(1)/(y(5)*params(31))+y(14)*params(25)*y(12)+y(15)*params(25)*y(11))+(1-params(27))*y(29);
T(33) = params(34)+y(3)*params(1)/(y(5)*params(31))+y(14)*params(25)*y(12)+y(15)*params(25)*y(11)+(1-params(27))*y(61);
T(34) = (y(2)-y(44))^2;
T(35) = (y(29))/(y(3));
T(36) = y(29)/y(3)-T(35);
T(37) = T(31)/(1+y(20));
end
