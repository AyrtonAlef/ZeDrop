function [rotm,a,b,h] = spherical_cap_dim(V,theta)
%SPHERICAL_CAP_DIM Return spherical cap dimensions for given values of drop volume and apparent CA
%   INPUT
% V - drop volume in mm³ or uL
% theta - apparent CA in degrees

%   OUTPUT
% rotm - sphere radius in mm
% a - wetted radius in mm
% b - distance from the center of the sphere to the center of the spherical cap in mm
% h - drop height in mm

Vcalc = 0;
tol = 1e-3; %Precision to find volume
diffV = 1;
r = 0; %Initial sphere radius
inc = 1;
while abs(diffV) > tol 
    if theta <= 90 %Se CA <= 90°
        beta = theta;
        a = r*sin(deg2rad(beta)); %Drop wetted radius
        b = r*cos(deg2rad(beta));
        h = r*(1-cos(deg2rad(beta))); %Drop height
    else %Se CA > 90°
        beta = 180 - theta;
        a = r*sin(deg2rad(beta)); %Drop wetted radius
        b = r*cos(deg2rad(beta));
        h = r*(1+cos(deg2rad(beta))); %Drop height
    end
    Vcalc = (1/3)*pi*(h^2)*(3*r-h); %Drop volume
    diffV = V-Vcalc;
    if diffV > 0
        rant = r;
        r = r + inc;
    else
        r = rant;
        inc = inc/2;
        r = r + inc;
    end
    rotm = rant;
    %fprintf('r = %.5f mm, h = %.4f mm and V = %.4f uL\n',rotm,h,Vcalc);
end 
end

