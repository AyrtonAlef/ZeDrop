function [dropVol] = dropVolumeSph(baseL,baseR,CA,scale)
% DROPVOLUMESPH Calculates the volume considering a spherical cap
%   INPUT
% contactDiameter - contact diameter of the drop in pixels 
% CA - contact angle in degrees
% scale - scal in pixels/mm

x0L=baseL(1);
y0L=baseL(2);
x0R=baseR(1);
y0R=baseR(2);
d_pxs = sqrt((x0R-x0L)^2+(y0R-y0L)^2); %Estimation of of contact diameter in pixels (d_pxs)
d = d_pxs/scale; %Calculation of contact diameter (d) in mm
a = d/2; % calculation of drop radius (a) in mm
if CA <= 90
    beta = CA;
    h = a*(1-cos(deg2rad(beta)))/sin(deg2rad(beta)); %Calculation of drop height (h) in mm
elseif CA > 90
    beta = 180 - CA;
    h = a*(1+cos(deg2rad(beta)))/sin(deg2rad(beta)); %Calculation of drop height (h) in mm
end
dropVol = (1/6)*pi*h*(3*a^2+h^2); %Estimation of drop volume in mm³ or uL
end

