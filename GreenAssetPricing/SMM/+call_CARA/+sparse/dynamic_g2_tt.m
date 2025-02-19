function [T_order, T] = dynamic_g2_tt(y, x, params, steady_state, T_order, T)
if T_order >= 2
    return
end
[T_order, T] = call_CARA.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
T_order = 2;
if size(T, 1) < 130
    T = [T; NaN(130 - size(T, 1), 1)];
end
T(77) = getPowerDeriv(T(3),1-params(17),2);
T(78) = params(31)*params(31)*getPowerDeriv(y(4)*params(31),params(1),2);
T(79) = (T(2)*T(2)*(-(params(40)*y(62)*T(78)))-(-(params(40)*y(62)*T(47)))*(T(2)*T(47)+T(2)*T(47)))/(T(2)*T(2)*T(2)*T(2));
T(80) = params(15)*getPowerDeriv(y(72),1-params(17),2);
T(81) = params(15)*getPowerDeriv(y(12),(-params(17)),2);
T(82) = getPowerDeriv(y(62),(-params(17)),2);
T(83) = (T(2)*T(2)*(-(params(40)*T(78)))-(-(params(40)*T(47)))*(T(2)*T(47)+T(2)*T(47)))/(T(2)*T(2)*T(2)*T(2));
T(84) = getPowerDeriv(T(6),1-params(17),2);
T(85) = T(49)*T(84);
T(86) = T(50)*T(83)+T(49)*T(85);
T(87) = getPowerDeriv(T(10),(-params(17)),2);
T(88) = (-(params(40)*T(52)))/(T(9)*T(9));
T(89) = params(31)*params(31)*getPowerDeriv(params(31)*y(64),params(1),2);
T(90) = (T(9)*T(9)*(-(params(40)*y(122)*T(89)))-(-(params(40)*y(122)*T(52)))*(T(9)*T(52)+T(9)*T(52)))/(T(9)*T(9)*T(9)*T(9));
T(91) = getPowerDeriv(y(66),params(14),2);
T(92) = params(43)*getPowerDeriv(y(118),params(44),2);
T(93) = (-((-y(67))*(y(3)+y(3))))/(y(3)*y(3)*y(3)*y(3));
T(94) = getPowerDeriv(T(14),(-params(13)),2);
T(95) = (-1)/(y(3)*y(3));
T(96) = (-((-y(127))*(y(63)+y(63))))/(y(63)*y(63)*y(63)*y(63));
T(97) = getPowerDeriv(T(17),1-params(13),2);
T(98) = T(45)*T(96)+T(44)*T(44)*T(97);
T(99) = (-1)/(y(63)*y(63));
T(100) = T(45)*T(99)+T(44)*T(60)*T(97);
T(101) = (-((-y(133))*(y(125)+y(125))))/(y(125)*y(125)*y(125)*y(125));
T(102) = T(60)*T(60)*T(97);
T(103) = getPowerDeriv(T(14),1-params(13),2);
T(104) = getPowerDeriv(y(65),params(28),2);
T(105) = (1+y(79))*(1+y(79))*(1+y(79))*(1+y(79));
T(106) = (-((-(y(12)*params(39)*(y(75)+y(80))))*(y(20)+y(20))))/(y(20)*y(20)*y(20)*y(20));
T(107) = (-(100*(-1)/(y(12)*y(12))));
T(108) = (-(params(31)*params(1)))/(params(31)*y(64)*params(31)*y(64));
T(109) = params(5)*1000*T(108);
T(110) = (-((-(params(31)*params(1)*y(122)))*(params(31)*params(31)*y(64)+params(31)*params(31)*y(64))))/(params(31)*y(64)*params(31)*y(64)*params(31)*y(64)*params(31)*y(64));
T(111) = (-(params(39)*(y(135)+y(140))))/(y(80)*y(80));
T(112) = (-(y(72)*params(39)))/(y(80)*y(80));
T(113) = (-((-((y(135)+y(140))*y(72)*params(39)))*(y(80)+y(80))))/(y(80)*y(80)*y(80)*y(80));
T(114) = (-((-y(80))*(y(75)+y(75))))/(y(75)*y(75)*y(75)*y(75));
T(115) = T(27)*T(114)-T(68)*T(68);
T(116) = (-(T(74)*T(74)));
T(117) = (-(T(36)*T(36)));
T(118) = (-(T(65)*T(65)));
T(119) = params(39)/y(80)/(1+y(79));
T(120) = (-((T(32)*T(119)-T(64)*T(69))/(T(32)*T(32))));
T(121) = (-(params(39)*(y(135)+y(140))/y(80)))/((1+y(79))*(1+y(79)));
T(122) = T(111)/(1+y(79));
T(123) = (-((-(T(69)*T(69)))/(T(32)*T(32))));
T(124) = (-(y(72)*params(39)/y(80)))/((1+y(79))*(1+y(79)));
T(125) = (-((T(32)*T(124)-T(69)*T(71))/(T(32)*T(32))));
T(126) = T(112)/(1+y(79));
T(127) = (-((T(32)*T(126)-T(69)*T(73))/(T(32)*T(32))));
T(128) = (-((-T(25))*(1+y(79)+1+y(79))))/T(105);
T(129) = (-T(72))/((1+y(79))*(1+y(79)));
T(130) = T(113)/(1+y(79));
end
