function [y, T] = dynamic_11(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  T(28)=params(15)*y(72)^(1-params(17));
  y(92)=T(28)*y(152)+(1-params(17))^(-1)*(params(40)*y(62)/T(6))^(1-params(17));
end
