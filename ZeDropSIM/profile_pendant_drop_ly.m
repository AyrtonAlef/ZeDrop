function [x_ly_cm,z_ly_cm,x_ly_sim_cm,z_ly_sim_cm,indx_last,indx_second,indx_deq,indx_ds,d_eq_cm,h_eq_cm,h_cm,d_s_cm,A_cm2,V_uL] = profile_pendant_drop_ly(b,c,r_h_cm,S_span,initialCo,tol)
%PENDANT_DROP_LY Calculation of pendant drop profile by the numerical integration
%of the Young-Laplace equation
%   
%   INPUT:
% b - curvature at the apex in cm-1
% c - capillary constant in cm-2
% r_h_cm - holder radius in cm
% S-span - step variable for ode45 solver
% initialCo - initial conditions for integration
% tol - tolerance to find the intersection points of the drop with the holder radius
%   OUTPUT:
% x_ly_cm - x coordinates of pendant drop profile in cm
% z_ly_cm - z coordinates of pendant drop profile in cm
% x_ly_sim_cm - x coordinates of pendant drop profile trimmed to the holder radius in cm
% z_ly-sim_cm - z coordinates of pendant drop profile trimmed to the holder radius in cm
% indx_last - index corresponding to the last(third) interception of the Laplace curve with the holder radius
% indx_second - index corresponding to second interception of the Laplace curve with the holder radius
% indx_deq - index corresponding to the maximum diameter
% indx_ds - index corresponding to the drop height desired (h = d_eq) 
% d_eq_cm - drop maximum diameter in cm
% h_eq_cm - distance from the drop apex to the horizontal position of drop maximum diameter in cm
% h_cm - drop height (distance from the drop apex to the needle tip) in cm
% d_s_cm - drop diameter from a d_eq distance from the drop apex in cm
% A_cm2 - drop surface area in cm²
% V_uL - drop volume in mm³ or uL

[S,Y] = ode45(@laplace_pendant_b_c,S_span,initialCo,[],b,c); % Second input might need to be increased if curve does not loop
Phi = Y(:,3)*(180/pi); % the third column of Y is the Phi data. Needs to be defined to find cutoff point below.
m = 1;
while Phi(m)>=0 %Finds where contact angle becomes negative
    m = m + 1;
end
% Output Variables
x_ly = Y(1:m,2);               % Drop width dimensionless
z_ly = Y(1:m,1);               % Drop height dimensionless
S_ly = S(1:m,1);               % Arc length dimensionless
Phi_ly = Y(1:m,3)*(180/pi);    % Contact angle in degrees
V_ly = Y(1:m,4);               % Drop volume dimensionless
A_ly = Y(1:m,5);               % Drop surface area dimensionless

%Conversion data to dimensions
x_ly_cm = x_ly/(c^.5);   % Drop width in cm
z_ly_cm = z_ly/(c^.5);   % Drop height in cm
S_ly_cm = S_ly/(c^.5);   % Arc length in cm
V_ly_uL = (V_ly/(c^(3/2)))*1000;    % Drop volume in mm³ or uL
A_ly_cm2 = A_ly/c;       % Drop surface area in cm²

%Find the position corresponding to the last(third) interception of the Laplace curve with the holder radius
indx_last = length(x_ly_cm);
while abs(x_ly_cm(indx_last) - r_h_cm) > tol
    indx_last = indx_last-1;
end
        
%Find the position corresponding to the maximum diameter
indx_deq = 1;
exit = 0;
while x_ly_cm(indx_deq) < x_ly_cm(indx_deq+1)
    indx_deq = indx_deq+1;
    if indx_deq+1 == length(x_ly_cm(1:indx_last)) %Stop while if the diameter of the drop reachs holder outer diameter
        exit = 1; 
        break
    end
end
        
%Find the position corresponding to second interception of the Laplace curve with the holder radius
indx_second = indx_deq;
while abs(x_ly_cm(indx_second) - r_h_cm) > tol && exit==0
    indx_second = indx_second + 1;
end

%Trim the simulated data to the second interception of the Laplace curve with the holder radius
x_ly_sim_cm = x_ly_cm(1:indx_second);   
z_ly_sim_cm = z_ly_cm(1:indx_second);
%S_ly_sim_cm = S_ly_cm(1:indx_second);
%Phi_ly_sim = Phi_ly(1:indx_second);
V_ly_sim_uL = V_ly_uL(1:indx_second);
A_ly_sim_cm2 = A_ly_cm2(1:indx_second);

%Drop properties
d_eq_cm = 2*x_ly_sim_cm(indx_deq); %drop maximum diameter in cm 
h_eq_cm = z_ly_sim_cm(indx_deq); %distance from the drop apex to the horizontal position of drop maximum diameter in cm
h_cm = z_ly_sim_cm(end); %drop height (distance from the drop apex to the needle tip) in cm
%indx_ds = 0;
if d_eq_cm > h_cm %drop height greater than maximum diameter
    d_s_cm = 0; %It is not possible to calculate d_s_cm
else
    indx_ds = find((abs(z_ly_sim_cm-d_eq_cm)) == min(abs(z_ly_sim_cm-d_eq_cm)),1,'last'); %Find the index of the drop height
    d_s_cm = 2*x_ly_sim_cm(indx_ds); %drop diameter at a d_eq distance from the drop apex in cm
end
A_cm2 = A_ly_sim_cm2(end);
V_uL = V_ly_sim_uL(end);
end

