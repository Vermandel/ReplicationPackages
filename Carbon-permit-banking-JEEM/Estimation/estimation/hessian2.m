function hessian_mat = hessian(func,x, gstep, varargin) % --*-- Unitary tests --*--

% Computes second order partial derivatives
%
% INPUTS
%    func        [string]   name of the function
%    x           [double]   vector, the Hessian of "func" is evaluated at x.
%    gstep       [double]   scalar, size of epsilon.
%    varargin    [void]     list of additional arguments for "func".
%
% OUTPUTS
%    hessian_mat [double]   Hessian matrix
%
% ALGORITHM
%    Uses Abramowitz and Stegun (1965) formulas 25.3.23
% \[
%     \frac{\partial^2 f_{0,0}}{\partial {x^2}} = \frac{1}{h^2}\left( f_{1,0} - 2f_{0,0} + f_{ - 1,0} \right)
% \]
% and 25.3.27 p. 884
%
% \[
%     \frac{\partial ^2f_{0,0}}{\partial x\partial y} = \frac{-1}{2h^2}\left(f_{1,0} + f_{-1,0} + f_{0,1} + f_{0,-1} - 2f_{0,0} - f_{1,1} - f_{-1,-1} \right)
% \]
%
% SPECIAL REQUIREMENTS
%    none
%

% Copyright (C) 2001-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if ~isa(func, 'function_handle')
    func = str2func(func);
end

n   = size(x,1);
h1  = max(abs(x), sqrt(gstep(1))*ones(n, 1))*eps^(1/6)*gstep(2);
h_1 = h1;
xh1 = x+h1;
h1  = xh1-x;
xh1 = x-h_1;
h_1 = x-xh1;
xh1 = x;
f0  = feval(func, x, varargin{:});
f1  = zeros(size(f0, 1), n);
f_1 = f1;

%for i=1:n
    %do step up
%    xh1(i)   = x(i)+h1(i);
%    f1(:,i)  = feval(func, xh1, varargin{:});
    %do step down
%    xh1(i)   = x(i)-h_1(i);
%    f_1(:,i) = feval(func, xh1, varargin{:});
%    %reset parameter
%    xh1(i)   = x(i);
%end

bayestopt_=varargin{6};

parfor i=1:n
    lnprior = priordens(x,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4,1);
    my_xh1 = xh1;
    %do step up
    my_xh1(i)   = x(i)+h1(i);
    f1(:,i)  = feval(func, my_xh1, varargin{:});
    %do step downxh1
    my_xh1(i)   = x(i)-h_1(i);
    f_1(:,i) = feval(func, my_xh1, varargin{:});
end

xh_1 = xh1;
temp = f1+f_1-f0*ones(1, n); %term f_(1,0)+f_(-1,0)-f_(0,0) used later

hessian_mat = zeros(size(f0,1), n*n);

for i=1:n
    if i > 1
        %fill symmetric part of Hessian based on previously computed results
        k = [i:n:n*(i-1)];
        hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1) = hessian_mat(:,k);
    end
    hessian_mat(:,(i-1)*n+i) = (f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i)); %formula 25.3.23
    
    myhess = nan(n,1);
    parfor j=i+1:n
        my_xh1 = xh1;
        my_xh_1 = xh_1;
        %step in up direction
        my_xh1(i) = x(i)+h1(i);
        my_xh1(j) = x(j)+h_1(j);
        %step in down direction
        my_xh_1(i) = x(i)-h1(i);
        my_xh_1(j) = x(j)-h_1(j);
        
        %hessian_mat(:,(i-1)*n+j) =-(-feval(func, my_xh1, varargin{:})-feval(func, my_xh_1, varargin{:})+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j)); %formula 25.3.27                                                                                                                                %reset grid points
        myhess(j) =-(-feval(func, my_xh1, varargin{:})-feval(func, my_xh_1, varargin{:})+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j)); %formula 25.3.27                             
    end
    j=i+1:n;
    hessian_mat(:,(i-1)*n+j)=myhess(j);
end
