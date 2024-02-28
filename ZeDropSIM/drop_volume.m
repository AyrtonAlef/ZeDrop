function [drop_vol] = drop_volume(z_ly_cm,x_ly_cm)
% This function finds the volume using washers scheme
%   INPUT:
% z_ly_l - z coordinates of the drop profile in cm
% x_ly_l - x coordinates of the drop profile

V = zeros(length(z_ly_cm)-1,1); % Size of Volume matrix
for n = 2:length(z_ly_cm)
    V(n-1) = pi*(((x_ly_cm(n)+x_ly_cm(n-1))/2)^2)*(z_ly_cm(n)-z_ly_cm(n-1));
end
drop_vol = sum(V)*1000; %conversion to mm³ or uL
end

