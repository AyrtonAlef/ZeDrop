function [CI,CA_circ, CA_pol, D, CP, circle_l, circle_r, drop_left_fit, drop_right_fit] = PolynomialFitting (ID,poly_ord,x_CPL,y_CPL,x_CPR,y_CPR,TiltAngle_S)
%POLYNOMIALFITTING Calculation of contact angle using polynomial fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function extracted from Dropen_V01.m code 
% (https://board.unimib.it/datasets/wzchzbm58p/3) and modified to analyze
% left and right drop profile edges separately. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
% ID: image
% poly_ord: polynomial order (the authors recommend: (i) for CA < 60° -> poly_ord = 2 and (ii) for CA > 60° -> poly_ord = 3)
% Apex_x and Apex_y: z and y coordinates of the drop apex
% x_CPL,y_CPL,x_CPR and y_CPR: coordinates of the left and right contact points, respectively.
% TiltAngle_S: baseline tilt angle
%   OUTPUT:
% CI: cropped rotated image with no black frame
% CA_circ: contact angle calculated by circle fitting
% CA_pol: contact angle calculated by polynomial fitting
% D: contact diameter in pixels
% circle_l: circle fitting parameters for the left drop profile 
% circle_r: circle fitting parameters for the right drop profile
% drop_left_fit: estimated coordinates from polynomial fitting for left drop edge
% drop_right_fit: estimated coordinates from polynomial fitting for right drop edge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Load vairables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_fit = 100; %number of fitting points
x_base = round((x_CPL+ x_CPR) /2);
y_base = round((y_CPL + y_CPR) /2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Rotate image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Irot_uncrop = imrotate(ID,TiltAngle_S,'bilinear');
%-- The rotated image is cropped to avoid black frame
crop_dim2 = floor((size(Irot_uncrop,2) - size(ID,2))/2);
crop_dim1 = floor((size(Irot_uncrop,1) - size(ID,1))/2);
%-- Cropped rotated image with no black frame
I_rot = Irot_uncrop(floor(2*crop_dim1+1):end-2*crop_dim1-1,2*crop_dim2+1:end-2*crop_dim2-1);
CI=I_rot;
y_baseline=x_base-1;

[CA_circ, CA_pol, D, CP, circle_l, circle_r, drop_left_fit, drop_right_fit] = ...
    Image2CA(CI, y_baseline, n_fit, poly_ord);
end

