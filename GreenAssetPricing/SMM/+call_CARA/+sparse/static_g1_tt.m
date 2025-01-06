function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = call_CARA.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 36
    T = [T; NaN(36 - size(T, 1), 1)];
end
T(27) = getPowerDeriv(T(4),1-params(17),1);
T(28) = getPowerDeriv(T(4),(-params(17)),1);
T(29) = getPowerDeriv(T(14),(-params(13)),1);
T(30) = getPowerDeriv(T(14),1-params(13),1);
T(31) = (-y(7))/(y(3)*y(3))*T(30);
T(32) = params(31)*getPowerDeriv(y(4)*params(31),params(1),1);
T(33) = params(23)/(1-params(13))*T(30)*1/y(3);
T(34) = params(15)*getPowerDeriv(y(12),(-params(17)),1);
T(35) = (y(20)*y(12)*params(39)-(y(15)+y(20))*y(12)*params(39))/(y(20)*y(20));
T(36) = (-(params(43)*getPowerDeriv(y(58),params(44),1)));
end
