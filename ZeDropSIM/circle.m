function [x,y] = circle(xc,yc,r,theta,N)
%CIRCLE Return circle coordinedes given radius, center and apparentCA
%   INPUT
% xc and yc - coordinates of the center
% r - radius
% N - number of points

% theta = linspace(0,2*pi,N);
phi_conv = 90-theta;
phi_min = deg2rad(phi_conv);
phi_max = deg2rad(180-phi_conv);
phi = linspace(phi_min,phi_max,N);
x = r*cos(phi)+xc;
y = r*sin(phi)+yc;
end

