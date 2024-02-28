function [xn,zn,g] = grayscale_edge(x_data,z_data,n_normal,g1,g2,W)
%GRAYSCALE_EDGE Calculate the grayscale profile along the edge
%   INPUT:
% x_data - x coordinates of the edge data
% z_data - z coordinates of the edge data
% n_normal - number of points normal to the edge
% g1 and g2 - corresponds to the plateau grey values at the two sides of the edge
% W - measure of the dge width

xn = zeros(1,length(x_data)*n_normal);
zn = zeros(1,length(z_data)*n_normal);
g = zeros(1,length(x_data)*n_normal);
scale_n = linspace(-floor(n_normal/2),floor(n_normal/2),n_normal);
y0 = 0;
for i=1:(length(z_data))
    if i == 1 %For the inital point in the edge
        dx = (x_data(i+1)-x_data(i)); %step 1: dv/dt
        dy =(z_data(i+1)-z_data(i));    
    elseif i == length(z_data) %For the last point in the edge
        dx = (x_data(i)-x_data(i-1)); %step 1: dv/dt
        dy =(z_data(i)-z_data(i-1));
    else     
        dx = (x_data(i+1)-x_data(i-1)); %step 1: dv/dt
        dy =(z_data(i+1)-z_data(i-1));
    end
    %Calculate normal to edge points
    dvdt = [dx;dy];
    dvdtRot = [-dy ;dx]; %step 2: rotate 90°
    dvdtRotNorm = dvdtRot/norm(dvdtRot); %step 3:scale it to magnitude 1
    k = 0;
    %Calculate normal point
    for j = (((i-1)*n_normal)+1):(i*n_normal) 
        k = k + 1;
        xn(j) = x_data(i) + (scale_n(k)*(W/2)*dvdtRotNorm(1));
        zn(j) = z_data(i) + (scale_n(k)*(W/2)*dvdtRotNorm(2));
        g(j) = ((g1-g2)/(1+exp((scale_n(k)-y0)/W)))+g2;
    end 
end
end

