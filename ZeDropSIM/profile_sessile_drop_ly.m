function [x_ly_sim_cm,z_ly_sim_cm,theta_sim,wetted_radius_cm,r_eq_cm,drop_height_cm,A_cm2,V_uL] = profile_sessile_drop_ly(b,c,theta,S_span,initialCo)
%PROFILE_SESSILE_DROP_LY Calculation of sessile drop profile by the numerical integration
%of the Young-Laplace equation given curvature at the apex (b), liquid
%capillary constant (c) and contact angle (theta)

%   INPUT:
% b - curvature at the apex in cm-1
% c - capillary constant in cm-2
% theta - apparent CA in degrees
% S-span - step variable for ode45 solver
% initialCo - initial conditions for integration
%   OUTPUT:
% x_ly_sim_cm - x coordinates of pendant drop profile trimmed to the holder radius in cm
% z_ly-sim_cm - z coordinates of pendant drop profile trimmed to the holder radius in cm
% wetted_radius_cm - raio molhado em cm
% r_eq_cm - raio equatorial (máximo) da gota em cm
% drop_height_cm - altura da gota em cm
% A_cm2 - drop surface area in cm²
% V_uL - drop volume in uL calculated concomitantly with the Laplace curve generation

[S,Y] = ode45(@laplace_b_c,S_span,initialCo,[],b,c); % Second input might need to be increased if curve does not loop
Z = Y(:,1); % the first column of Y is the Z data. Needs to be defined to find cutoff point below.
% This loop stores the value, i, where the Laplace-Young data starts to 'loop' out of control.
m = 1;
while Z(m)<Z(m+1) % Finds where contact angle is 180 degrees (cutoff point)
    m = m+1;
end
% Output Variables
x_ly = Y(1:m,2);               % Drop width dimensionless
z_ly = Y(1:m,1);               % Drop height dimensionless
S_ly = S(1:m,1);               % Arc length dimensionless
Phi_ly = Y(1:m,3)*(180/pi);    % Contact angle in degrees
V_ly = Y(1:m,4);               % Drop volume dimensionless
A_ly = Y(1:m,5);               % Drop surface area dimensionless
   
%Trim the simulated data to the desired range corresponding to theta
indx = find((abs(Phi_ly-theta)) == min(abs(Phi_ly-theta)),1,'last'); %Find the index of thetha desired
x_ly_sim = x_ly(1:indx);   
z_ly_sim = z_ly(1:indx);
S_ly_sim = S_ly(1:indx);
Phi_ly_sim = Phi_ly(1:indx);
V_ly_sim = V_ly(1:indx);
A_ly_sim = A_ly(1:indx);

%Conversion trimmed data to dimensions
x_ly_sim_cm = x_ly_sim/(c^.5);   % Drop width in cm
z_ly_sim_cm = z_ly_sim/(c^.5);   % Drop height in cm
S_ly_sim_cm = S_ly_sim/(c^.5);   % Arc length in cm
V_ly_sim_uL = (V_ly_sim/(c^(3/2)))*1000;    % Drop volume in mm³ or uL
A_ly_sim_cm2 = A_ly_sim/c; 

%Drop properties
theta_sim = Phi_ly_sim(end);
wetted_radius_cm = x_ly_sim_cm(end);
%wetted_radius_pxs(i,j) = x_ly_sim_pxs(end);
r_eq_cm = max(x_ly_sim_cm);
%r_eq_pxs(i,j) = max(x_ly_sim_pxs);
drop_height_cm = z_ly_sim_cm(end);
%drop_height_pxs(i,j) = z_ly_sim_pxs(end);
A_cm2 = A_ly_sim_cm2(end);
V_uL = V_ly_sim_uL(end); %Volume calculated concomitantly with the Laplace curve generation
end

