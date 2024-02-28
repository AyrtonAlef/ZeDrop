function [dy]=laplace_r0_c(s,y,r0,c)
%LAPLACE defines the ordinary differential equations to be solved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Derek Fultz and Russ Stacy
%7−24−07
% z=drop height
% x=distance from axis to drop interface
% phi=contact angle
% r0=radius of curvature at apex
% s=arc length
%
% NON−DIMENSIONALIZE Z,X,S,B USING C^(1/2)
% R0=r0*c^(1/2)
% X=x*c^(1/2)
% Z=z*c^(1/2)
% S=s*c^(1/2)
%
% No need to define X,Z,S for equations
%
% KEY OF VARIABLE DEFINITIONS FOR SYSTEM
% Z=y(1); Z'=dy(1); Z' is with respect to S
% X=y(2); X'=dy(2); X' is with respect to S
% phi=y(3); phi'=dy(3); phi' is with respect to S
%
% Z'=sin(phi)
% X'=cos(phi)
% phi'=2/R0+Z−(sin(phi)/X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy = zeros(3,1); % a column vector
R0 = r0*c^.5; % non−dimensionalized curvature at apex
dy(1)=sin(y(3)); % first geometrical relationships
dy(2)=cos(y(3)); % second geometrical relationships
dy(3)=(2/R0)+y(1)-(sin(y(3))/y(2)); % Young-Laplace equation for axisymmetric drops
end

