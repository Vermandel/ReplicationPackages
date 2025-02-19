function [y, T, residual, g1] = dynamic_4(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(19, 1);
  residual(1)=(y(119)*params(4))-(params(43)*params(44)*y(118)^(params(44)-1));
  T(1)=exp((-params(25))*(y(4)*params(31)-params(3)));
  T(2)=y(66)^params(14);
  T(3)=params(22)^(1-params(14));
  residual(2)=(y(74))-(T(1)*params(10)*y(76)*T(2)*T(3));
  residual(3)=(params(31)*y(64))-(params(31)*y(4)*(1-params(27))+params(5)*(y(77)+params(9)));
  residual(4)=(y(77))-(y(73)*params(4)*(1-y(118))*params(41));
  residual(5)=(y(72))-(1+y(74));
  T(4)=params(34)+params(1)*y(122)/(params(31)*y(64))+y(133)*params(25)*y(131)+y(134)*params(25)*y(130)+(1-params(27))*y(180);
  residual(6)=(y(120))-(y(72)*y(121)*T(4));
  T(5)=(y(127)/y(63))^(1-params(13));
  residual(7)=(y(69))-(params(39)*y(121)*y(129)*(1-params(19)+params(23)/(1-params(13))*T(5)-params(23)*T(5)+params(24))+y(121)*y(131)*params(28)*y(133)/y(125));
  T(6)=(y(4)*params(31))^params(1);
  T(7)=params(40)/T(6);
  T(8)=T(7)^(1-params(17));
  T(9)=y(62)^(-params(17));
  residual(8)=(y(68))-(T(8)*T(9)-T(7)*y(91)*(1-params(12)));
  T(10)=params(15)*y(12)^(-params(17));
  T(11)=(params(31)*y(64))^params(1);
  T(12)=params(40)*y(122)/T(11);
  T(13)=params(12)*y(151)+T(12)^(-params(17));
  residual(9)=(y(91))-(T(10)*T(13));
  residual(10)=(y(70))-(params(39)*y(121)*y(130)*(1+(1-params(14))*y(134))+y(121)*y(131)*(1-params(28))*y(133));
  residual(11)=(params(28)*y(71)*y(73)/y(65))-(params(14)*y(70)*y(74)/y(66));
  residual(12)=(params(5)*y(120))-(y(119));
  T(14)=1-params(43)*y(118)^params(44);
  residual(13)=(y(71))-(T(14)-y(119)*(1-y(118))*params(4));
  T(15)=(y(67)/y(3))^(-params(13));
  residual(14)=(1)-(y(69)*params(23)*T(15));
  T(16)=params(24)+params(23)/(1-params(13))*(y(67)/y(3))^(1-params(13));
  residual(15)=(y(72)*y(63))-(y(3)*(1-params(19))+y(3)*T(16));
  T(17)=y(65)^params(28);
  T(18)=params(21)^(1-params(28));
  residual(16)=(y(73))-(y(76)*T(1)*params(6)*params(33)*T(17)*T(18));
  residual(17)=(y(3))-(y(66)+y(65));
  residual(18)=(y(73)*T(14))-(y(62)+y(67)+params(29)*params(42));
  residual(19)=(y(61))-(T(10)*y(68)/y(8));
  T(19)=getPowerDeriv(T(12),(-params(17)),1);
  T(20)=getPowerDeriv(y(67)/y(3),(-params(13)),1);
  T(21)=getPowerDeriv(y(67)/y(3),1-params(13),1);
  T(22)=getPowerDeriv(y(127)/y(63),1-params(13),1);
  T(23)=(-y(127))/(y(63)*y(63))*T(22);
  T(24)=(-(params(40)*params(31)*getPowerDeriv(y(4)*params(31),params(1),1)))/(T(6)*T(6));
  T(25)=params(15)*getPowerDeriv(y(12),(-params(17)),1);
  T(26)=(-(params(43)*getPowerDeriv(y(118),params(44),1)));
if nargout > 3
    g1_v = NaN(78, 1);
g1_v(1)=(-(T(13)*T(25)));
g1_v(2)=(-(y(68)*T(25)/y(8)));
g1_v(3)=(-(T(3)*T(2)*y(76)*params(10)*T(1)*params(31)*(-params(25))));
g1_v(4)=(-(params(31)*(1-params(27))));
g1_v(5)=(-(T(9)*T(24)*getPowerDeriv(T(7),1-params(17),1)-(1-params(12))*y(91)*T(24)));
g1_v(6)=(-(T(18)*T(17)*y(76)*params(6)*params(33)*T(1)*params(31)*(-params(25))));
g1_v(7)=(-(y(69)*params(23)*(-y(67))/(y(3)*y(3))*T(20)));
g1_v(8)=(-(1-params(19)+T(16)+y(3)*params(23)/(1-params(13))*(-y(67))/(y(3)*y(3))*T(21)));
g1_v(9)=1;
g1_v(10)=(-((-(T(10)*y(68)))/(y(8)*y(8))));
g1_v(11)=params(4);
g1_v(12)=(-1);
g1_v(13)=(1-y(118))*params(4);
g1_v(14)=(-(T(3)*T(1)*params(10)*y(76)*getPowerDeriv(y(66),params(14),1)));
g1_v(15)=(-(params(14)*y(70)*(-y(74))/(y(66)*y(66))));
g1_v(16)=(-1);
g1_v(17)=(-params(5));
g1_v(18)=1;
g1_v(19)=(-(params(43)*params(44)*getPowerDeriv(y(118),params(44)-1,1)));
g1_v(20)=(-(y(73)*params(4)*(-params(41))));
g1_v(21)=(-(T(26)-params(4)*(-y(119))));
g1_v(22)=y(73)*T(26);
g1_v(23)=1;
g1_v(24)=(-(y(121)*T(4)));
g1_v(25)=y(63);
g1_v(26)=params(31);
g1_v(27)=(-(y(72)*y(121)*(-(params(31)*params(1)*y(122)))/(params(31)*y(64)*params(31)*y(64))));
g1_v(28)=(-(T(10)*T(19)*(-(params(40)*y(122)*params(31)*getPowerDeriv(params(31)*y(64),params(1),1)))/(T(11)*T(11))));
g1_v(29)=(-(params(39)*y(121)*y(129)*(params(23)/(1-params(13))*T(23)-params(23)*T(23))));
g1_v(30)=y(72);
g1_v(31)=1;
g1_v(32)=(-(T(10)/y(8)));
g1_v(33)=T(7)*(1-params(12));
g1_v(34)=1;
g1_v(35)=1;
g1_v(36)=(-(params(14)*y(74)/y(66)));
g1_v(37)=1;
g1_v(38)=(-1);
g1_v(39)=(-(params(14)*y(70)*1/y(66)));
g1_v(40)=1;
g1_v(41)=params(5);
g1_v(42)=params(28)*y(73)/y(65);
g1_v(43)=1;
g1_v(44)=1;
g1_v(45)=(-(params(23)*T(15)));
g1_v(46)=(-(y(69)*params(23)*T(20)*1/y(3)));
g1_v(47)=(-(y(3)*params(23)/(1-params(13))*T(21)*1/y(3)));
g1_v(48)=(-1);
g1_v(49)=(-(params(4)*(1-y(118))*params(41)));
g1_v(50)=params(28)*y(71)*1/y(65);
g1_v(51)=1;
g1_v(52)=T(14);
g1_v(53)=params(28)*y(71)*(-y(73))/(y(65)*y(65));
g1_v(54)=(-(T(18)*y(76)*T(1)*params(6)*params(33)*getPowerDeriv(y(65),params(28),1)));
g1_v(55)=(-1);
g1_v(56)=(-(T(8)*getPowerDeriv(y(62),(-params(17)),1)));
g1_v(57)=(-1);
g1_v(58)=1;
g1_v(59)=(-(T(10)*params(12)));
g1_v(60)=(-(y(72)*y(121)*params(25)*y(134)));
g1_v(61)=(-(params(39)*y(121)*(1+(1-params(14))*y(134))));
g1_v(62)=(-(y(72)*y(121)*params(25)*y(130)));
g1_v(63)=(-((1-params(14))*params(39)*y(121)*y(130)));
g1_v(64)=(-((1-params(27))*y(72)*y(121)));
g1_v(65)=(-(y(72)*y(121)*params(25)*y(133)));
g1_v(66)=(-(y(133)/y(125)*y(121)*params(28)));
g1_v(67)=(-(y(133)*y(121)*(1-params(28))));
g1_v(68)=(-(params(39)*y(121)*(1-params(19)+params(23)/(1-params(13))*T(5)-params(23)*T(5)+params(24))));
g1_v(69)=(-(params(39)*y(121)*y(129)*(params(23)/(1-params(13))*T(22)*1/y(63)-params(23)*T(22)*1/y(63))));
g1_v(70)=(-(y(72)*y(121)*params(25)*y(131)));
g1_v(71)=(-(y(121)*y(131)*params(28)*1/y(125)));
g1_v(72)=(-(y(121)*y(131)*(1-params(28))));
g1_v(73)=(-(y(121)*y(131)*params(28)*(-y(133))/(y(125)*y(125))));
g1_v(74)=(-(y(72)*y(121)*params(1)/(params(31)*y(64))));
g1_v(75)=(-(T(10)*params(40)/T(11)*T(19)));
g1_v(76)=(-(y(72)*T(4)));
g1_v(77)=(-((1-params(19)+params(23)/(1-params(13))*T(5)-params(23)*T(5)+params(24))*params(39)*y(129)+y(133)/y(125)*y(131)*params(28)));
g1_v(78)=(-((1+(1-params(14))*y(134))*params(39)*y(130)+y(133)*y(131)*(1-params(28))));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 19, 57);
end
end
