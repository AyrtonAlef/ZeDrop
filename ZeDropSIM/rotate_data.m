function [x_rot,y_rot] = rotate_data(x,y,x_center,y_center,phi)
%ROTATE_DATA Rotate data counterclockwise in the plane by a certain angle about a specific point in the
%same plane
%   INPUT
% x and y - data to be rotated
% xc and yc - center of rotation
% phi - rotation angle in degrees

%Definition of the rotation angle
phi_rad = deg2rad(phi);
R = [cos(phi_rad) -sin(phi_rad) 0; sin(phi_rad) cos(phi_rad) 0; 0 0 1];
%Define affine transformation for translation
a = [1 0 x_center;0 1 y_center; 0 0 1];
c = [1 0 -x_center;0 1 -y_center; 0 0 1];
M = a*R*c;
for i=1:length(x)
    rot(:,i) = M*[x(i) y(i) 1]';
end
%Pick out the vectors of rotated x- and y-data
x_rot = rot(1,:);
y_rot = rot(2,:);
%{
%Make a plot
plot(x,y,'k-',x_rot, y_rot, 'r-', x_center, y_center, 'bo');
axis equal
%}
end

