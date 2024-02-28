function [I,d_eq_cm,h_eq_cm,d_s_cm,h_cm,A_cm2,V_uL] = synthetical_pendantdrop_yl(b,c,r_h_cm,S_span,initialCo,resolution,scale_cm,needle)
%Construction of a syhtetical sessile drop binary image
%   INPUTS:
% b - curvature at apex in cm-1
% c - cappilary constant in cm-2
% r_h_cm - holder outer radius in cm
% S-span - integration interval
% initialCo - initial condition for integration
% resolution - image resolution [res_h res_v]
% scale_cm - image scale in cm
% needle_diameter_cm - outer needle diameter in mm
% needle - generation of synthetical pendant drop image with needle (=1) or without (=0)

color = 'black'; %Cor dos elementos a serem inseridos (preto)
opacity = 1; %Opacidade dos elementos a serem inseridos na imagem
tol = 10^-5; %tolerance to find the intersection points of the drop with the holder radius

[S,Y] = ode45(@laplace_pendant_b_c,S_span,initialCo,[],b,c); % Second input might need to be increased if curve does not loop
Phi = Y(:,3)*(180/pi); % the third column of Y is the Phi data. Needs to be defined to find cutoff point below.
m = 1;
while Phi(m)>=0 %Finds where contact angle becomes negative
    m = m+1;
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
while x_ly_cm(indx_deq)<x_ly_cm(indx_deq+1)
    indx_deq = indx_deq+1;
    if indx_deq+1 == length(x_ly_cm(1:indx_last)) %Stop while if the daimeter of the drop reachs holder outer diameter
        exit = 1; 
        break
    end
end
        
%Find the position corresponding to second interception of the Laplace curve with the holder radius
indx_second = indx_deq;
while abs(x_ly_cm(indx_second) - r_h_cm) > tol && exit==0
    indx_second = indx_second + 1;
end

%Trim the simulated data to the second interception of th Laplace curve with the holder radius
x_ly_sim_cm = x_ly_cm(1:indx_second);   
z_ly_sim_cm = z_ly_cm(1:indx_second);
S_ly_sim_cm = S_ly_cm(1:indx_second);
Phi_ly_sim = Phi_ly(1:indx_second);
V_ly_sim_uL = V_ly_uL(1:indx_second);
A_ly_sim_cm2 = A_ly_cm2(1:indx_second);

%Conversion to pixels
x_ly_sim_pxs = x_ly_sim_cm*scale_cm;   % Drop width in pixels
z_ly_sim_pxs = z_ly_sim_cm*scale_cm;   % Drop height in pixels
%S_ly_sim_pxs = S_ly_sim_cm*scale_cm;   % Arc length in pixels

%Drop properties
d_eq_cm = 2*x_ly_sim_cm(indx_deq); %drop maximum diameter in cm
h_eq_cm = z_ly_sim_cm(indx_deq); %distance from the drop apex to the horizontal position of maximum diameter in cm
h_cm = z_ly_sim_cm(end); %distance from the drop apex to the needle tip in cm
%indx_ds = 0;
if d_eq_cm > h_cm %drop height greater than its maximum diameter
    d_s_cm = 0; %It is not possible to calculate d_s_cm
else
    indx_ds = find((abs(z_ly_sim_cm-d_eq_cm)) == min(abs(z_ly_sim_cm-d_eq_cm)),1,'last'); %Find the index of rh desired
    d_s_cm = 2*x_ly_sim_cm(indx_ds); %drop diameter from a d_eq distance from the drop apex
end
A_cm2 = A_ly_sim_cm2(end);
V_uL = V_ly_sim_uL(end);

%Construction of synthetical image
res_v = resolution(2);
res_h = resolution(1);
xc = round(res_h/2); %x coordinate of the drop center in pixels
M = ones(res_v,res_h)*255; %Array with image resolution size
I = uint8(M); %Image creation
x_plot_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)

if needle %Construction of pendant drop with needle
    zc = round(res_v/2);
    z_plot = -z_ly_sim_pxs+z_ly_sim_pxs(end)+zc-(z_ly_sim_pxs(end)/2);   
    needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
    I = insertShape(I,'FilledRectangle',[(xc-needle_diameter_pxs/2) 1 needle_diameter_pxs (zc-(z_ly_sim_pxs(end)/2))],'Color',color,'Opacity',opacity);%Insert rectangle that corresponds to the needle
else %Construction of pendant drop without needle
    z_plot = -z_ly_sim_pxs + z_ly_sim_pxs(end);
end
z_plot_lap_full = [wrev(z_plot);z_plot]; %Creation of drop profile (z coordinate)
ly_polygon = [x_plot_lap_full z_plot_lap_full];
I = insertShape(I,'FilledPolygon',ly_polygon,'Color',color,'Opacity',opacity); %Insert polygon that corresponds to the drop
end

