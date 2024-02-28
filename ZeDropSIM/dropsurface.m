function [X,Y,Z,gof] = dropsurface(x,z,n)
%CREATE_MESH Create a revolution surface from drop profile

%   INPUT
% x - x coordinates of the drop profile
% z - z coordinate of the drop profile
% n - number of points around the cilinder circunference and along z-axis
%   OUTPUT
% X,Y,Z - coordinates of the 3D drop profile
% gof - goodness-of-fit of cubic interpolation

%new_z_ly_sim_cm = -z_ly_sim_cm;
%[f,gof]=fit(z,x,'cubicinterp'); %Fit a piecewise cubic interpolation to data
[f,gof]=fit(z,x,'smoothingspline');
%gof - godness of fit
%plot(f,z,x) %Plot fit and data points
z_equal = linspace(0,z(end),n); %Generate an equally spaced z array
r = f(z_equal); %Generate an equally spaced x data according to z array
[X,Y,Z] = cylinder(r,n); %Generate the 3D matrix
Z = Z*z(end); %Adjust the height of the 3D drop profile
%mesh(X,Y,Z) %Ploting 3D mesh
end

