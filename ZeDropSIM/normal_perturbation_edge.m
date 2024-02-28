function [xn_pert,zn_pert,xn_pert_main,zn_pert_main] = normal_perturbation_edge(xn_data,zn_data,pert,n_normal)
%NORMAL_PERTURBATION_EDGE Induce a normal perturbation (addition) in the edge
%   INPUT:
% xn_data - x coordinates of the edge data
% zn_data - z coordinates of the edge data
% pert - magnitude of the induced perturbation
% n_normal - number of points normal to the edge
%   OUTPUT:
% xn_pert - x coordinates of all perturbed data
% zn_pert - z coordinates of all perturbed data


xn_pert = zeros(1,length(xn_data));
zn_pert = zeros(1,length(zn_data));
xn_pert_main = zeros(1,(length(xn_data)/n_normal));
zn_pert_main = zeros(1,(length(zn_data)/n_normal));
%{
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
    
    %Calculate normal point
    rand_pert = -pert+2*pert*rand(1); %Generate uniformly distributed random numbers
    x_pert(i) = x_data(i) + rand_pert*dvdtRotNorm(1);
    z_pert(i) = z_data(i) + rand_pert*dvdtRotNorm(2);
end
%}
k = 0;
for i=1:(length(xn_data)/n_normal)
    if i == 1 %For the inital point in the edge
        %(i-1)*n_normal + ceil(n_normal/2)
        dx = xn_data(i*n_normal+ceil(n_normal/2))- xn_data((i-1)*n_normal+ceil(n_normal/2)); %step 1: dv/dt
        dy = zn_data(i*n_normal+ceil(n_normal/2)) - zn_data((i-1)*n_normal+ceil(n_normal/2));    
    elseif i == (length(xn_data)/n_normal) %For the last point in the edge
        dx = xn_data((i-1)*n_normal+ceil(n_normal/2)) - xn_data((i-2)*n_normal+ceil(n_normal/2)); %step 1: dv/dt
        dy = zn_data((i-1)*n_normal+ceil(n_normal/2)) - zn_data((i-2)*n_normal+ceil(n_normal/2));
    else     
        dx = xn_data(i*n_normal+ceil(n_normal/2)) - xn_data((i-2)*n_normal+ceil(n_normal/2)); %step 1: dv/dt
        dy = zn_data(i*n_normal+ceil(n_normal/2)) - zn_data((i-2)*n_normal+ceil(n_normal/2));
    end
    %Calculate normal to edge points
    dvdt = [dx;dy];
    dvdtRot = [-dy ;dx]; %step 2: rotate 90°
    dvdtRotNorm = dvdtRot/norm(dvdtRot); %step 3:scale it to magnitude 1
    rand_pert = -pert+2*pert*rand(1); %Generate uniformly distributed random numbers
    
    %Random normal perturbation in the main drop profile 
    k = k+1;
    xn_pert_main(k) = xn_data((i-1)*n_normal+ceil(n_normal/2));
    zn_pert_main(k) = zn_data((i-1)*n_normal+ceil(n_normal/2));
    
    %Random normal perturbation at each point of the edge 
    for j = ((i-1)*n_normal+1):(i*n_normal)
        xn_pert(j) = xn_data(j) + rand_pert*dvdtRotNorm(1);
        zn_pert(j) = zn_data(j) + rand_pert*dvdtRotNorm(2);
    end
end
end

