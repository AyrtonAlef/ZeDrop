function p = polyfitcon(x,y,n,y0)
%POLYFITCON Fits a polynomial function to a dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function extracted from Dropen_V01.m code 
% (https://board.unimib.it/datasets/wzchzbm58p/3). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input:
% x and y: coordinates of the data to be fitted
% n: polynomial order 
% y0: constant of the polynomial equation
%   Output:
% p: coefficients of the polynomial function

%Load data
x = x(:);
y = y(:);
V(:,n) = zeros(length(x),1,class(x));
%Formation of the polynomial function
for j = 1:n
   V(:,j) = x.^(n+1-j);
end
%Fitting polynomial function (least square method)
p_unknown = V\(y - y0);
p = [p_unknown; y0];

end

