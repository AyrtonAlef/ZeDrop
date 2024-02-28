function [dropVol] = dropVolumeWasher(edgeL,edgeR,baseL,baseR,scale)
% DROPVOLUMEWASHER Calculates the volume using washers scheme
%   INPUT
% z_ly_l - z coordinates of the drop profile in cm
% x_ly_l - x coordinates of the drop profile

left=[edgeL.x,edgeL.y];
right=[edgeR.x,edgeR.y];
%{
x0L=baseL(1);
y0L=baseL(2);
x0R=baseR(1);
y0R=baseR(2);
%}
indexL = find_index(left,baseR,baseL);
indexR = find_index(right,baseR,baseL);
leftNorm_mm = abs(left(1:indexL,:)-left(1,:))/scale; %Changing reference and trasforming to mm
rightNorm_mm = abs(right(1:indexR,:)-right(1,:))/scale; %Changing reference and trasforming to mm

%Estimating drop volume using left drop profile
vLeft = zeros(length(leftNorm_mm)-1,1);
for n = 2:length(leftNorm_mm)
    vLeft(n-1) = pi*(((leftNorm_mm(n,1)+leftNorm_mm(n-1,1))/2)^2)*(leftNorm_mm(n,2)-leftNorm_mm(n-1,2));
end
dropVolLeft = sum(vLeft); %Estimating drop volume using left edge
%Estimating drop volume using right drop profile
vRight = zeros(length(rightNorm_mm)-1,1);
for n = 2:length(rightNorm_mm)
    vRight(n-1) = pi*(((rightNorm_mm(n,1)+rightNorm_mm(n-1,1))/2)^2)*(rightNorm_mm(n,2)-rightNorm_mm(n-1,2));
end
dropVolRight = sum(vRight); %Estimating drop volume using left edge
dropVol = (dropVolLeft+dropVolRight)/2;
%{
V = zeros(length(z_ly_cm)-1,1); % Size of Volume matrix
for n = 2:length(z_ly_cm)
    V(n-1) = pi*(((x_ly_cm(n)+x_ly_cm(n-1))/2)^2)*(z_ly_cm(n)-z_ly_cm(n-1));
end
drop_vol = sum(V)*1000; %conversion to mm³ or uL
%}
end

