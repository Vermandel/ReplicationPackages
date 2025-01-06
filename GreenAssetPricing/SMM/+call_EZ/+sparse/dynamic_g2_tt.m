function [T_order, T] = dynamic_g2_tt(y, x, params, steady_state, T_order, T)
if T_order >= 2
    return
end
[T_order, T] = call_EZ.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
T_order = 2;
if size(T, 1) < 261
    T = [T; NaN(261 - size(T, 1), 1)];
end
T(142) = params(7)*params(7)*getPowerDeriv(params(7)*y(64),T(5),2);
T(143) = T(4)*T(142);
T(144) = getPowerDeriv(T(7),T(8),2);
T(145) = getPowerDeriv(T(3),params(11),2);
T(146) = (-(params(40)*T(88)))/(T(2)*T(2));
T(147) = T(1)*(T(73)*T(146)+T(72)*T(89)*T(145));
T(148) = params(31)*params(31)*getPowerDeriv(y(5)*params(31),params(1),2);
T(149) = (T(2)*T(2)*(-(params(40)*y(66)*T(148)))-(-(params(40)*y(66)*T(88)))*(T(2)*T(88)+T(2)*T(88)))/(T(2)*T(2)*T(2)*T(2));
T(150) = params(15)*getPowerDeriv(y(76),params(11),2);
T(151) = getPowerDeriv(T(3),params(11)-1,2);
T(152) = (T(2)*T(2)*(-T(148))-(-T(88))*(T(2)*T(88)+T(2)*T(88)))/(T(2)*T(2)*T(2)*T(2));
T(153) = T(9)*T(152);
T(154) = (-T(88))/(T(2)*T(2))*T(127);
T(155) = params(32)*params(32)*getPowerDeriv(y(96)*params(32),1-params(11),2);
T(156) = T(10)*T(155);
T(157) = params(7)*params(7)*getPowerDeriv(params(7)*y(64),T(14),2);
T(158) = T(16)*T(13)*T(157);
T(159) = getPowerDeriv(T(19),params(11)-1,2);
T(160) = (-(params(40)*T(96)))/(T(18)*T(18));
T(161) = params(31)*params(31)*getPowerDeriv(params(31)*y(68),params(1),2);
T(162) = (T(18)*T(18)*(-(params(40)*y(129)*T(161)))-(-(params(40)*y(129)*T(96)))*(T(18)*T(96)+T(18)*T(96)))/(T(18)*T(18)*T(18)*T(18));
T(163) = T(1)*(T(97)*T(97)*T(159)+T(79)*T(162));
T(164) = params(15)*getPowerDeriv(y(76),params(11)-1,2);
T(165) = params(32)*params(32)*getPowerDeriv(params(32)*y(159),1-params(17)-params(11),2);
T(166) = params(7)*params(7)*getPowerDeriv(params(7)*y(1),T(22),2);
T(167) = (-(y(96)*params(32)*T(166)));
T(168) = (T(23)*T(23)*T(167)-T(58)*(T(23)*T(57)+T(23)*T(57)))/(T(23)*T(23)*T(23)*T(23));
T(169) = getPowerDeriv(T(24),1-params(17)-params(11),2);
T(170) = T(51)*T(21)*(T(60)*T(168)+T(59)*T(59)*T(169));
T(171) = (-(params(32)*T(57)))/(T(23)*T(23));
T(172) = T(51)*T(21)*(T(60)*T(171)+T(59)*T(129)*T(169));
T(173) = getPowerDeriv(T(54),params(11)-1,2);
T(174) = T(70)*(-(T(12)*T(1)*T(67)*T(67)*T(173)))-(-(T(12)*T(69)))*(T(55)*T(69)+T(55)*T(69));
T(175) = T(70)*T(70);
T(176) = T(174)/T(175);
T(177) = (-(T(69)*T(1)*T(72)*T(75)))/T(70);
T(178) = (-(T(69)*T(1)*T(75)*T(89)))/T(70);
T(179) = (-(T(69)*T(1)*(-T(75))))/T(70);
T(180) = T(55)*T(124)+T(55)*T(124);
T(181) = (-((-(T(12)*T(69)))*T(180)))/T(175);
T(182) = (-(T(69)*T(124)))/T(70);
T(183) = (-(params(40)*T(137)))/(T(53)*T(53));
T(184) = T(55)*T(139)+T(55)*T(139);
T(185) = T(70)*(-(T(12)*T(1)*(T(68)*T(183)+T(67)*T(138)*T(173))))-(-(T(12)*T(69)))*T(184);
T(186) = T(185)/T(175);
T(187) = T(55)*T(1)*(-T(68))+T(55)*T(1)*(-T(68));
T(188) = (T(70)*(-(T(12)*T(1)*T(67)*(-T(173))))-(-(T(12)*T(69)))*T(187))/T(175);
T(189) = T(1)*T(72)*T(72)*T(151)/T(55);
T(190) = T(1)*(T(75)*T(146)+T(72)*T(89)*T(151))/T(55);
T(191) = T(1)*T(72)*(-T(151))/T(55);
T(192) = (-(T(1)*T(72)*T(75)*T(124)))/T(70);
T(193) = (-(T(1)*T(72)*T(75)*T(139)))/T(70);
T(194) = (-(T(1)*T(72)*T(75)*T(1)*(-T(68))))/T(70);
T(195) = getPowerDeriv(T(50),(-params(1)),2);
T(196) = T(1)*(T(89)*T(89)*T(151)+T(75)*T(149))/T(55);
T(197) = T(1)*T(89)*(-T(151))/T(55);
T(198) = (-(T(1)*T(75)*T(89)*T(124)))/T(70);
T(199) = (-(params(31)*params(31)))/(params(31)*y(62)*params(31)*y(62));
T(200) = T(26)*(T(92)*T(199)+T(91)*T(134)*T(195));
T(201) = (-(T(1)*T(75)*T(89)*T(139)))/T(70);
T(202) = (-(T(1)*T(75)*T(89)*T(1)*(-T(68))))/T(70);
T(203) = params(15)*getPowerDeriv(y(13),params(11)-1,2);
T(204) = T(1)*T(151)/T(55);
T(205) = (-(T(1)*(-T(75))*T(124)))/T(70);
T(206) = (-(T(1)*(-T(75))*T(139)))/T(70);
T(207) = (-(T(1)*(-T(75))*T(1)*(-T(68))))/T(70);
T(208) = (-((-(T(12)*T(124)))*T(180)))/T(175);
T(209) = (-(T(124)*T(124)))/T(70);
T(210) = (-((-(T(12)*T(124)))*T(184)))/T(175);
T(211) = (-((-(T(12)*T(124)))*T(187)))/T(175);
T(212) = (-(T(124)*T(139)))/T(70);
T(213) = (-(T(124)*T(1)*(-T(68))))/T(70);
T(214) = (-((-(params(31)*y(5)*params(31)))*(params(31)*params(31)*y(62)+params(31)*params(31)*y(62))))/(params(31)*y(62)*params(31)*y(62)*params(31)*y(62)*params(31)*y(62));
T(215) = T(134)*T(134)*T(195)+T(92)*T(214);
T(216) = params(31)*params(31)*getPowerDeriv(params(31)*y(62),params(1),2);
T(217) = (T(53)*T(53)*(-(params(40)*y(3)*T(216)))-(-(params(40)*y(3)*T(137)))*(T(53)*T(137)+T(53)*T(137)))/(T(53)*T(53)*T(53)*T(53));
T(218) = (T(70)*(-(T(12)*T(1)*(T(138)*T(138)*T(173)+T(68)*T(217))))-(-(T(12)*T(139)))*T(184))/T(175);
T(219) = (T(70)*(-(T(12)*T(1)*T(138)*(-T(173))))-(-(T(12)*T(139)))*T(187))/T(175);
T(220) = (T(70)*(-(T(12)*T(1)*T(173)))-(-(T(12)*T(1)*(-T(68))))*T(187))/T(175);
T(221) = getPowerDeriv(y(70),params(14),2);
T(222) = params(43)*getPowerDeriv(y(122),params(44),2);
T(223) = (-((-y(71))*(y(4)+y(4))))/(y(4)*y(4)*y(4)*y(4));
T(224) = getPowerDeriv(T(30),(-params(13)),2);
T(225) = (-1)/(y(4)*y(4));
T(226) = (-((-y(134))*(y(67)+y(67))))/(y(67)*y(67)*y(67)*y(67));
T(227) = getPowerDeriv(T(33),1-params(13),2);
T(228) = T(86)*T(226)+T(85)*T(85)*T(227);
T(229) = (-1)/(y(67)*y(67));
T(230) = T(86)*T(229)+T(85)*T(105)*T(227);
T(231) = (-((-y(140))*(y(132)+y(132))))/(y(132)*y(132)*y(132)*y(132));
T(232) = T(105)*T(105)*T(227);
T(233) = getPowerDeriv(T(30),1-params(13),2);
T(234) = getPowerDeriv(y(69),params(28),2);
T(235) = (1+y(83))*(1+y(83))*(1+y(83))*(1+y(83));
T(236) = (-((-(y(13)*params(39)*(y(79)+y(84))))*(y(21)+y(21))))/(y(21)*y(21)*y(21)*y(21));
T(237) = (-(100*(-1)/(y(13)*y(13))));
T(238) = (-(params(31)*params(1)))/(params(31)*y(68)*params(31)*y(68));
T(239) = params(5)*1000*T(238);
T(240) = (-((-(params(31)*params(1)*y(129)))*(params(31)*params(31)*y(68)+params(31)*params(31)*y(68))))/(params(31)*y(68)*params(31)*y(68)*params(31)*y(68)*params(31)*y(68));
T(241) = (-(params(39)*(y(142)+y(147))))/(y(84)*y(84));
T(242) = (-(y(76)*params(39)))/(y(84)*y(84));
T(243) = (-((-((y(142)+y(147))*y(76)*params(39)))*(y(84)+y(84))))/(y(84)*y(84)*y(84)*y(84));
T(244) = (-((-y(84))*(y(79)+y(79))))/(y(79)*y(79)*y(79)*y(79));
T(245) = T(43)*T(244)-T(115)*T(115);
T(246) = (-(T(121)*T(121)));
T(247) = (-(T(77)*T(77)));
T(248) = (-(T(112)*T(112)));
T(249) = params(39)/y(84)/(1+y(83));
T(250) = (-((T(48)*T(249)-T(111)*T(116))/(T(48)*T(48))));
T(251) = (-(params(39)*(y(142)+y(147))/y(84)))/((1+y(83))*(1+y(83)));
T(252) = T(241)/(1+y(83));
T(253) = (-((-(T(116)*T(116)))/(T(48)*T(48))));
T(254) = (-(y(76)*params(39)/y(84)))/((1+y(83))*(1+y(83)));
T(255) = (-((T(48)*T(254)-T(116)*T(118))/(T(48)*T(48))));
T(256) = T(242)/(1+y(83));
T(257) = (-((T(48)*T(256)-T(116)*T(120))/(T(48)*T(48))));
T(258) = (-((-T(41))*(1+y(83)+1+y(83))))/T(235);
T(259) = (-T(119))/((1+y(83))*(1+y(83)));
T(260) = T(243)/(1+y(83));
T(261) = (-(T(123)*T(123)));
end
