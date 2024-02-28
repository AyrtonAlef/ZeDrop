function [dy]=laplace_pendant_b_c(s,y,b,c)
%LAPLACE defines the ordinary differential equations to be solved for pendant drops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified code from Konduru 2010:
%(i) modified the input variable r0 (curvature radius at the apex) to b
%(ii) include the simultaneous calculus of drop volume (V) and drop superficial area (A)
%(iii) generation of a pendant drop

%Equations presented in del Río et al 1997
%System of ordinary differential equation as a function of arc length s
%dx/ds = cos(phi)
%dz/ds = sin(phi)
%dphi/ds = 2*b+c*z-(sin(phi)/x)
%dv/ds = pi*(x^2)*sin(phi)
%da/ds = 2*pi*x

%Derek Fultz and Russ Stacy
%7−24−07
% z=drop height
% x=distance from axis to drop interface
% phi=contact angle
% b=curvature at apex
% s=arc length
% v=volume of the drop
% a=surface area of the drop

% NON−DIMENSIONALIZE Z,X,S,B,V,A USING c
% B=b/c^(1/2)
% X=x*c^(1/2)
% Z=z*c^(1/2)
% S=s*c^(1/2)
% V=v*c^(3/2)
% A=a*c

% No need to define X,Z,S for equations
%
% KEY OF VARIABLE DEFINITIONS FOR SYSTEM
% Z=y(1); Z'=dy(1); Z' is with respect to S
% X=y(2); X'=dy(2); X' is with respect to S
% phi=y(3); phi'=dy(3); phi' is with respect to S
% V=y(4); V'=dy(4); V' is with respect to S
% A=y(5); A'=dy(5); A' is with respect to S

% Z'=sin(phi)
% X'=cos(phi)
% phi'=(2*B)+Z−(sin(phi)/X)
% V'=pi*(X^2)*sin(phi)*(c^(3/4))
% A'=2*pi*X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy = zeros(3,1); % a column vector
B = b/c^.5; % non−dimensionalized curvature at apex
dy(1) = sin(y(3)); % first geometrical relationships
dy(2) = cos(y(3)); % second geometrical relationships
dy(3) = (2*B)-y(1)-(sin(y(3))/y(2)); % Young-Laplace equation for axisymmetric drops
%dy(4)=pi*(y(2)^2)*sin(y(3))*(c^(3/4)); %Calculation of volume
dy(4) = pi*(y(2)^2)*sin(y(3)); %Calculation of volume
dy(5) = 2*pi*y(2); %Calculation of area
end

