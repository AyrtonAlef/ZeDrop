function [p,pos0L,pos0R,pos0L_cut,pos0R_cut,index_left,index_right,index_cut_left,index_cut_right] = plane_mod(bound_l,bound_r,percent)
%PLANE_MOD Finds the z-plane for the image.
% The z-plane is defined as the surface in which the drop sets.
% left_plane_ind: the indices that fall on the left side of the apex and that are within +- 1 pixel of the first element on the left.
% right_plane_ind: the indices that fall on the right side of the apex and that are within +- 1 pixel of the first element on the right.
% left_plane: values at the indices for the left.
% right_plane: values at the indices for the right.
% left_avg and right_avg: averages of the left_plane and right_plane vectors, respectively.
%   INPUT:
% bound_l: contain the coordinates of the left drop profile + solidsubstrate
% bound_l(:,1): x coordinates of the left profile (drop profile + surface substrate)
% bound_l(:,2) z coordinates of the lçeft profile (drop profile + surface substrate)
% bound_r: contain the coordinates of the right drop profile + solidsubstrate
% bound_r(:,1): x coordinates of the right profile (drop profile + surface substrate)
% bound_r(:,2) z coordinates of the right profile (drop profile + surface substrate)
% percent: percent of z-plane height where the contact points are going to be defined
%   OUTPUT:
% p: average of left_avg and right_avg. Used as the z plane being the mean of both sides of the drop.
% index_left: indice of the left contact point.
% index_right: indice of the right contact point.
% index_cut_left: indice of the left contact point using cutoff value.
% index_cut_right: indice of the right contact point using cutoff value.
% pos0L: coordinates of the left contact point.
% pos0R: coordinates of the right contact point. 
% pos0L_cut: coordinates of the left contact point using cutoff value.
% pos0R_cut: coordinates of the right contact point using cutoff value. 


left = bound_l(end,2); %Pick z coordinate of the first point in the left boundary
right = bound_r(end,2); %Pick up z coordinate of the first point in the right boundary
plane_limit=20; % compare left and right sides of z-plane. If not within limit end program and print "Image is not level."
% This is to ensure that the left and right z-planes are within a specific range of each other to guard against being unlevel.
if abs(left-right)>plane_limit
    fprintf('Image is not level! ');
    fprintf('left = %f ',left);
    fprintf('right = %f \n',right);
end
left_plane_ind = find(bound_l(:,2)>=left-1 & bound_l(:,2)<=left+1); %Identifying indices of the points in the boundary that are within +1 pixel of the extremely left surface substrate point
right_plane_ind = find(bound_r(:,2)>=right-1 & bound_r(:,2)<=right+1); %Identifying indices of the points in the boundary that are within +1 pixel of the extremely left surface substrate point
left_plane = bound_l(left_plane_ind,2);
right_plane = bound_r(right_plane_ind,2);
left_avg = mean(left_plane);
right_avg = mean(right_plane);
lr = [left_avg right_avg];
p = round(mean(lr)); %z plane in pixels

%Identification of contact points
index_left = min(left_plane_ind);
index_right = min(right_plane_ind);

x0L = bound_l(index_left,1);
z0L = bound_l(index_left,2);
pos0L = [x0L, z0L];
x0R = bound_r(index_right,1);
z0R = bound_r(index_right,2);
pos0R = [x0R,z0R];

%Identification of contact points (using cutoff)
cutoff = round(percent*p);
ind_left = find(bound_l(:,2)<cutoff);
ind_right = find(bound_r(:,2)<cutoff);
index_cut_left = ind_left(end);
index_cut_right = ind_right(end);

x0L_cut = bound_l(index_cut_left,1);
z0L_cut = bound_l(index_cut_left,2);
pos0L_cut = [x0L_cut, z0L_cut];
x0R_cut = bound_r(index_cut_right,1);
z0R_cut = bound_r(index_cut_right,2);
pos0R_cut = [x0R_cut, z0R_cut];

%{
x_ind=round(mean(find(bound(:,2)==z_apex))); % x indice of apex
left=z(1); % plane_l = surface on left of drop
right=z(end); % plane_r = surface on right of drop
plane_limit=20; % compare left and right sides of z-plane. If not within limit end program and print "Image is not level."
% This is to ensure that the left and right z-planes are within a specific range of each other to guard against being unlevel.
if abs(left-right)>plane_limit
    fprintf('Image is not level! nn');
    fprintf('left=%fnn',left);
    fprintf('right=%fnn',right);
end
plane_ind_left=find(z>=left-1 & z<=left+1);
plane_ind_right=find(z>=right-1 & z<=right+1);
left_plane_ind=plane_ind_left(find(plane_ind_left<x_ind));
right_plane_ind=plane_ind_right(find(plane_ind_right>x_ind));
left_plane=z(left_plane_ind);
right_plane=z(right_plane_ind);
left_avg=mean(left_plane);
right_avg=mean(right_plane);
lr=[left_avg right_avg];
p=round(mean(lr)); % z plane in pixels
p_cm=p/s; % z plane in cm
p_nd=p_cm*const^.5; % z plane non-dimensionalized
% To find wetted radius, find the maximum of the left_plane_ind(x_value)and the minimum of the right_plane_ind(x_value) and evalute the x and z data at that index.
% The evaluation takes place in height interface function.
index_left=max(left_plane_ind);
index_right=min(right_plane_ind);
%}
end

