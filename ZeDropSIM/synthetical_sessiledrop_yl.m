function [I,theta_sim,r_w_cm,r_w_pxs,r_eq_cm,r_eq_pxs,h_cm,h_pxs,A_cm2,V_uL] = synthetical_sessiledrop_yl(b,c,theta,S_span,initialCo,resolution,hsubstrato,scale_cm,needle_diameter_mm,needle_in)
%Construction of a syhtetical sessile drop binary image
%   INPUTS:
% b - curvature at apex in cm-1
% c - cappilary constant in cm-2
% S-span - integration interval
% initialCo - initial condition for integration
% resolution - image resolution [res_h res_v]
% hsubstrato - height of solute substrate in pxs
% scale_cm - image scale in cm
% needle_diameter_cm - outer needle diameter in mm
% needle_in - generation of synthetical sessile drop image with needle_in (=1) or needle_out (=0)

color = 'black'; %Color of the inserted elements
opacity = 1; %Opacity of the inserted elements

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
A_ly_sim_cm2 = A_ly_sim/c;       % Drop surface area in cm²
        
%Conversion to pixels
x_ly_sim_pxs = x_ly_sim_cm*scale_cm;   % Drop width in pixels
z_ly_sim_pxs = z_ly_sim_cm*scale_cm;   % Drop height in pixels
S_ly_sim_pxs = S_ly_sim_cm*scale_cm;   % Arc length in pixels

%Drop properties
theta_sim = Phi_ly_sim(end);
r_w_cm = x_ly_sim_cm(end);
r_w_pxs = x_ly_sim_pxs(end);
r_eq_cm = max(x_ly_sim_cm);
r_eq_pxs = max(x_ly_sim_pxs);
h_cm = z_ly_sim_cm(end);
h_pxs = z_ly_sim_pxs(end);
A_cm2 = A_ly_sim_cm2(end);
V_uL = V_ly_sim_uL(end); %Volume calculated concomitantly with the Laplace curve generation
        
%Construction of synthetical image
res_v = resolution(2);
res_h = resolution(1);
xc = round(res_h/2); %x coordinate of the drop center in pixels
M = ones(res_v,res_h)*255; %Array with image resolution size
I = uint8(M); %Image creation
I = insertShape(I,'FilledRectangle',[1 (res_v-hsubstrato) res_h hsubstrato],'Color',color,'Opacity',opacity); %Insert rectangle that corresponds to the substrate

x_plot_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)
z_plot = z_ly_sim_pxs - (hsubstrato + h_pxs - res_v);
z_plot_lap_full = [wrev(z_plot);z_plot]; %Creation of drop profile (z coordinate)
ly_polygon = [x_plot_lap_full z_plot_lap_full];
I = insertShape(I,'FilledPolygon',ly_polygon,'Color',color,'Opacity',opacity); %Insert polygon that corresponds to the drop

if needle_in %Construction of needle
    needle_diameter_pxs = (needle_diameter_mm/10)*scale_cm; %Needle outer diameter in pixels
    I = insertShape(I,'FilledRectangle',[(xc-needle_diameter_pxs/2) 1 needle_diameter_pxs (z_plot(end)-h_pxs/2)],'Color',color,'Opacity',opacity);%Insert rectangle that corresponds to the needle
end  
end

