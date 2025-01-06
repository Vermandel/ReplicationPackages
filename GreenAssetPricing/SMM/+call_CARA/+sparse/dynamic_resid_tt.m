function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 33
    T = [T; NaN(33 - size(T, 1), 1)];
end
T(1) = (1-params(17))^(-1);
T(2) = (y(4)*params(31))^params(1);
T(3) = params(40)*y(62)/T(2);
T(4) = params(15)*y(72)^(1-params(17));
T(5) = params(15)*y(12)^(-params(17));
T(6) = params(40)/T(2);
T(7) = T(6)^(1-params(17));
T(8) = y(62)^(-params(17));
T(9) = (params(31)*y(64))^params(1);
T(10) = params(40)*y(122)/T(9);
T(11) = y(66)^params(14);
T(12) = params(22)^(1-params(14));
T(13) = 1-params(43)*y(118)^params(44);
T(14) = y(67)/y(3);
T(15) = T(14)^(-params(13));
T(16) = params(23)/(1-params(13));
T(17) = y(127)/y(63);
T(18) = T(17)^(1-params(13));
T(19) = y(133)/y(125);
T(20) = params(24)+T(16)*T(14)^(1-params(13));
T(21) = y(65)^params(28);
T(22) = params(21)^(1-params(28));
T(23) = params(5)*1000*(params(34)+params(1)*y(122)/(params(31)*y(64))+y(133)*params(25)*y(131)+y(134)*params(25)*y(130))+(1-params(27))*y(148);
T(24) = params(34)+params(1)*y(122)/(params(31)*y(64))+y(133)*params(25)*y(131)+y(134)*params(25)*y(130)+(1-params(27))*y(180);
T(25) = (y(135)+y(140))*y(72)*params(39)/y(80);
T(26) = (y(121)-y(103))^2;
T(27) = y(80)/y(75);
T(28) = (steady_state(28))/(steady_state(2));
T(29) = y(88)/(steady_state(28));
T(30) = y(62)/(steady_state(2));
T(31) = y(73)/(steady_state(13));
T(32) = T(25)/(1+y(79));
T(33) = params(12)*y(151)+T(10)^(-params(17));
end
