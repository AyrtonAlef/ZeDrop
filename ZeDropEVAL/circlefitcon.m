function [xc,yc,r] = circlefitcon(x,y)
%CIRCLEFIT Fit a circle to a dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function extracted from Dropen_V01.m code 
% (https://board.unimib.it/datasets/wzchzbm58p/3). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
% x and y: coordinates of the data to be fitted
%   OUTPUT:
% xc and yc: coordinates of the circle center
% r: circle radius

%Load data
x = x(:);
y = y(:);
V(:,1) = x;
V(:,2) = y;
V(:,3) = 1;
q = -(x.^2 + y.^2);

%Estimation of circle coefficients
circle_coeff = V\q;
a = circle_coeff(1);
b = circle_coeff(2);
c = circle_coeff(3);
%Determination of circle parameters
xc = -a/2;
yc = -b/2;
r = sqrt(a^2/4 + b^2/4 - c);
end

