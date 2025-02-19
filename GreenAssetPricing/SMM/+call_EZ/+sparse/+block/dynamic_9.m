function [y, T] = dynamic_9(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(114)=100*(y(113)-(1+y(83)));
  y(112)=(y(92)/y(66)-(steady_state(29))/(steady_state(3)))/((steady_state(29))/(steady_state(3)));
  y(110)=log(y(84)/y(79));
  y(109)=(y(128)-y(107))^3;
  y(103)=log(y(92)/(steady_state(29)));
  y(104)=log(y(66)/(steady_state(3)));
  y(105)=log(y(77)/(steady_state(14)));
  y(115)=100*y(104);
  y(120)=log((y(142)+y(147))*y(76)*params(39)/y(84)/(1+y(83)));
  y(121)=log(1+y(83));
  y(119)=log(y(94)/(steady_state(31)));
  y(117)=log(y(123));
  y(108)=(y(128)-y(107))^2;
  y(106)=y(128);
  y(118)=100*((y(142)+y(147))*y(76)*params(39)/y(84)-(1+y(83)));
  y(111)=y(13)*params(39)*(y(79)+y(84))/y(21)-1;
  y(101)=params(31)*y(68);
  y(100)=100*(log(y(81))-log(y(18)));
  T(45)=log(y(13));
  y(99)=100*(T(45)+log(y(79))-log(y(16)));
  y(98)=100*(log(y(73))-log(y(10)));
  y(97)=y(96)*params(32);
  y(93)=y(66)*params(1)*1000/(y(5)*params(31));
  y(88)=100*(T(45)+log(y(71))-log(y(8)));
  y(87)=T(45)+log(y(66))-log(y(3));
  y(86)=100*(T(45)+log(y(77))-log(y(14)));
  y(91)=y(70)*100/y(69);
  y(89)=y(66)*100/y(77);
  y(90)=y(83)*100;
  y(85)=(y(13)*params(39)*(y(79)+y(84))/y(21)-(1+y(83)))*100;
  y(72)=(y(96)*params(32))^(1-params(11))*1/T(4)*T(14);
end
