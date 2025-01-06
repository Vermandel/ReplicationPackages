function [y, T, residual, g1] = static_12(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(1, 1);
  T(16)=params(15)*y(12)^(1-params(17));
  T(17)=(y(4)*params(31))^params(1);
  T(18)=params(40)*y(2)/T(17);
  residual(1)=(y(32))-(y(32)*T(16)+(1-params(17))^(-1)*T(18)^(1-params(17)));
if nargout > 3
    g1_v = NaN(1, 1);
g1_v(1)=1-T(16);
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 1, 1);
end
end
