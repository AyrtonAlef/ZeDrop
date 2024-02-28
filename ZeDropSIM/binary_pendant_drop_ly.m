function [Ibw] = binary_pendant_drop_ly(x_ly_sim_cm,z_ly_sim_cm,r_h_cm,resolution,scale_cm,needle)
%BINARY_PENDANT_DROP_LY %Creation of a syhtetical binary pendant drop image
%   
%   INPUT:
% x_ly_sim_cm - x coordinates of drop profile in cm
% z_ly_sim_cm - z coordinates of drop profile in cm
% r_h_cm - holder outer radius in cm
% resolution - image resolution [res_h res_v]
% scale_cm - image scale in cm
% needle - generation of synthetical pendant drop image with needle (=1) or without (=0)

%   OUTPUT
%Ibw - binary image

%Conversion to pixels
x_ly_sim_pxs = x_ly_sim_cm*scale_cm;   % Drop width in pixels
z_ly_sim_pxs = z_ly_sim_cm*scale_cm;   % Drop height in pixels

%Construction of synthetical image
color = 'black'; %Color of the inserted elements (drop and solid substrate)
opacity = 1; %Opacity of the inserted elements (drop and solid substrate)
res_v = resolution(2);
res_h = resolution(1);
xc = round(res_h/2); %x coordinate of the image center in pixels
M = ones(res_v,res_h)*255; %Array with image resolution size
I = uint8(M); %Criação da imagem
x_plot_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)
if needle %Creation of pendant drop with needle
    zc = round(res_v/2);
    z_plot = -z_ly_sim_pxs+z_ly_sim_pxs(end)+zc-(z_ly_sim_pxs(end)/2);   
    needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
    I = insertShape(I,'FilledRectangle',[(xc-needle_diameter_pxs/2) 1 needle_diameter_pxs (zc-(z_ly_sim_pxs(end)/2))],'Color',color,'Opacity',opacity);%Insert rectangle that corresponds to the needle
else %Creation of pendant drop without needle
    z_plot = -z_ly_sim_pxs + z_ly_sim_pxs(end);
end
z_plot_lap_full = [wrev(z_plot);z_plot]; %Creation of drop profile (z coordinate)
ly_polygon = [x_plot_lap_full z_plot_lap_full];
I = insertShape(I,'FilledPolygon',ly_polygon,'Color',color,'Opacity',opacity); %Insert polygon that corresponds to the drops

Igray = rgb2gray(I); %Convert image to grayscale
Ibw = imbinarize(Igray); %Binarize image
end

