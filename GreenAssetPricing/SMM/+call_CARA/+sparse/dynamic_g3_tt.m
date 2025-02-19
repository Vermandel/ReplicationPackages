function [T_order, T] = dynamic_g3_tt(y, x, params, steady_state, T_order, T)
if T_order >= 3
    return
end
[T_order, T] = call_CARA.sparse.dynamic_g2_tt(y, x, params, steady_state, T_order, T);
T_order = 3;
if size(T, 1) < 158
    T = [T; NaN(158 - size(T, 1), 1)];
end
T(131) = getPowerDeriv(T(3),1-params(17),3);
T(132) = params(31)*params(31)*params(31)*getPowerDeriv(y(4)*params(31),params(1),3);
T(133) = params(15)*getPowerDeriv(y(12),(-params(17)),3);
T(134) = (T(2)*T(2)*T(2)*T(2)*((T(2)*T(47)+T(2)*T(47))*(-(params(40)*T(78)))+T(2)*T(2)*(-(params(40)*T(132)))-((T(2)*T(47)+T(2)*T(47))*(-(params(40)*T(78)))+(-(params(40)*T(47)))*(T(47)*T(47)+T(2)*T(78)+T(47)*T(47)+T(2)*T(78))))-(T(2)*T(2)*(-(params(40)*T(78)))-(-(params(40)*T(47)))*(T(2)*T(47)+T(2)*T(47)))*(T(2)*T(2)*(T(2)*T(47)+T(2)*T(47))+T(2)*T(2)*(T(2)*T(47)+T(2)*T(47))))/(T(2)*T(2)*T(2)*T(2)*T(2)*T(2)*T(2)*T(2));
T(135) = getPowerDeriv(T(10),(-params(17)),3);
T(136) = params(43)*getPowerDeriv(y(118),params(44),3);
T(137) = (y(3)*y(3)*y(3)*y(3)*(-(2*(-y(67))))-(-((-y(67))*(y(3)+y(3))))*(y(3)*y(3)*(y(3)+y(3))+y(3)*y(3)*(y(3)+y(3))))/(y(3)*y(3)*y(3)*y(3)*y(3)*y(3)*y(3)*y(3));
T(138) = getPowerDeriv(T(14),(-params(13)),3);
T(139) = getPowerDeriv(T(17),1-params(13),3);
T(140) = T(96)*T(44)*T(97)+T(45)*(y(63)*y(63)*y(63)*y(63)*(-(2*(-y(127))))-(-((-y(127))*(y(63)+y(63))))*(y(63)*y(63)*(y(63)+y(63))+y(63)*y(63)*(y(63)+y(63))))/(y(63)*y(63)*y(63)*y(63)*y(63)*y(63)*y(63)*y(63))+T(96)*T(44)*T(97)+T(44)*(T(96)*T(97)+T(44)*T(44)*T(139));
T(141) = T(96)*T(60)*T(97)+T(45)*(y(63)+y(63))/(y(63)*y(63)*y(63)*y(63))+T(44)*T(97)*T(99)+T(44)*(T(97)*T(99)+T(44)*T(60)*T(139));
T(142) = getPowerDeriv(T(14),1-params(13),3);
T(143) = (-(100*(-((-(y(12)*params(39)))*(y(20)+y(20))))/(y(20)*y(20)*y(20)*y(20))));
T(144) = (y(20)*y(20)*y(20)*y(20)*(-(2*(-(y(12)*params(39)*(y(75)+y(80))))))-(-((-(y(12)*params(39)*(y(75)+y(80))))*(y(20)+y(20))))*(y(20)*y(20)*(y(20)+y(20))+y(20)*y(20)*(y(20)+y(20))))/(y(20)*y(20)*y(20)*y(20)*y(20)*y(20)*y(20)*y(20));
T(145) = (-(100*(y(12)+y(12))/(y(12)*y(12)*y(12)*y(12))));
T(146) = (-((-(params(31)*params(1)))*(params(31)*params(31)*y(64)+params(31)*params(31)*y(64))))/(params(31)*y(64)*params(31)*y(64)*params(31)*y(64)*params(31)*y(64));
T(147) = (params(31)*y(64)*params(31)*y(64)*params(31)*y(64)*params(31)*y(64)*(-((-(params(31)*params(1)*y(122)))*(params(31)*params(31)+params(31)*params(31))))-(-((-(params(31)*params(1)*y(122)))*(params(31)*params(31)*y(64)+params(31)*params(31)*y(64))))*(params(31)*y(64)*params(31)*y(64)*(params(31)*params(31)*y(64)+params(31)*params(31)*y(64))+params(31)*y(64)*params(31)*y(64)*(params(31)*params(31)*y(64)+params(31)*params(31)*y(64))))/(params(31)*y(64)*params(31)*y(64)*params(31)*y(64)*params(31)*y(64)*params(31)*y(64)*params(31)*y(64)*params(31)*y(64)*params(31)*y(64));
T(148) = (-(100*(-((-(y(72)*params(39)))*(y(80)+y(80))))/(y(80)*y(80)*y(80)*y(80))));
T(149) = (y(80)*y(80)*y(80)*y(80)*(-(2*(-((y(135)+y(140))*y(72)*params(39)))))-(-((-((y(135)+y(140))*y(72)*params(39)))*(y(80)+y(80))))*(y(80)*y(80)*(y(80)+y(80))+y(80)*y(80)*(y(80)+y(80))))/(y(80)*y(80)*y(80)*y(80)*y(80)*y(80)*y(80)*y(80));
T(150) = T(68)*T(114);
T(151) = T(32)*T(32)*T(32)*T(32);
T(152) = (-((T(32)*T(32)*(-(T(64)*T(119)+T(64)*T(119)))-(-(T(64)*T(64)))*(T(32)*T(69)+T(32)*T(69)))/T(151)));
T(153) = (-((-((T(32)*T(119)-T(64)*T(69))*(T(32)*T(69)+T(32)*T(69))))/T(151)));
T(154) = T(32)*(-(params(39)/y(80)))/((1+y(79))*(1+y(79)));
T(155) = T(32)*(-params(39))/(y(80)*y(80))/(1+y(79));
T(156) = (-((-((-(T(69)*T(69)))*(T(32)*T(69)+T(32)*T(69))))/T(151)));
T(157) = T(32)*(-((1+y(79)+1+y(79))*(-(y(72)*params(39)/y(80)))))/T(105);
T(158) = T(32)*(-((-(y(72)*params(39)))*(y(80)+y(80))))/(y(80)*y(80)*y(80)*y(80))/(1+y(79));
end
