function [y, T, residual, g1] = static_6(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(1, 1);
  residual(1)=(y(29))-(y(13)*y(2)*(params(5)*1000*(params(34)+y(3)*params(1)/(y(5)*params(31))+y(14)*params(25)*y(12)+y(15)*params(25)*y(11))+(1-params(27))*y(29)));
if nargout > 3
    g1_v = NaN(1, 1);
g1_v(1)=1-(1-params(27))*y(13)*y(2);
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 1, 1);
end
end
