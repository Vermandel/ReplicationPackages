function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = call_CARA.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 76
    T = [T; NaN(76 - size(T, 1), 1)];
end
T(34) = getPowerDeriv(T(3),1-params(17),1);
T(35) = getPowerDeriv(y(62),(-params(17)),1);
T(36) = 1/(steady_state(2));
T(37) = params(40)/T(9);
T(38) = getPowerDeriv(T(10),(-params(17)),1);
T(39) = params(1)/(params(31)*y(64));
T(40) = params(5)*1000*T(39);
T(41) = (-y(67))/(y(3)*y(3));
T(42) = getPowerDeriv(T(14),(-params(13)),1);
T(43) = getPowerDeriv(T(14),1-params(13),1);
T(44) = (-y(127))/(y(63)*y(63));
T(45) = getPowerDeriv(T(17),1-params(13),1);
T(46) = T(44)*T(45);
T(47) = params(31)*getPowerDeriv(y(4)*params(31),params(1),1);
T(48) = (-(params(40)*y(62)*T(47)))/(T(2)*T(2));
T(49) = (-(params(40)*T(47)))/(T(2)*T(2));
T(50) = getPowerDeriv(T(6),1-params(17),1);
T(51) = T(49)*T(50);
T(52) = params(31)*getPowerDeriv(params(31)*y(64),params(1),1);
T(53) = (-(params(40)*y(122)*T(52)))/(T(9)*T(9));
T(54) = (-(params(31)*params(1)*y(122)))/(params(31)*y(64)*params(31)*y(64));
T(55) = params(5)*1000*T(54);
T(56) = getPowerDeriv(y(65),params(28),1);
T(57) = (-y(133))/(y(125)*y(125));
T(58) = getPowerDeriv(y(66),params(14),1);
T(59) = 1/y(3);
T(60) = 1/y(63);
T(61) = T(16)*T(45)*T(60)-params(23)*T(45)*T(60);
T(62) = params(15)*getPowerDeriv(y(12),(-params(17)),1);
T(63) = params(15)*getPowerDeriv(y(72),1-params(17),1);
T(64) = params(39)*(y(135)+y(140))/y(80)/(1+y(79));
T(65) = 1/(steady_state(13));
T(66) = 1/y(125);
T(67) = 1/y(75);
T(68) = (-y(80))/(y(75)*y(75));
T(69) = y(72)*params(39)/y(80)/(1+y(79));
T(70) = (-(T(69)/T(32)));
T(71) = (-T(25))/((1+y(79))*(1+y(79)));
T(72) = (-((y(135)+y(140))*y(72)*params(39)))/(y(80)*y(80));
T(73) = T(72)/(1+y(79));
T(74) = 1/(steady_state(28));
T(75) = params(43)*getPowerDeriv(y(118),params(44),1);
T(76) = (-T(75));
end
