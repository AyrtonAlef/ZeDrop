function [Ibw] = binary_sessile_drop_ly(x_ly_sim_cm,z_ly_sim_cm,hsubstrato,resolution,scale_cm,needle_diameter_mm,needle_in)
%BINARY_PENDANT_DROP_LY %Creation of a synthetical binary sessile drop image
%   
%   INPUT:
% x_ly_sim_cm - x coordinates of drop profile in cm
% z_ly_sim_cm - z coordinates of drop profile in cm
% hsubstrato - height of solid susbtrate in pixels
% resolution - image resolution [res_h res_v]
% scale_cm - image scale in cm
% needle_in - generation of synthetical sessile drop image with needle in (=1) or without needle (=0)

%   OUTPUT:
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
I = uint8(M); %Image creation
I = insertShape(I,'FilledRectangle',[1 (res_v-hsubstrato) res_h hsubstrato],'Color',color,'Opacity',opacity); %Insert rectangle that corresponds to the solid substrate
x_plot_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)
%z_plot = z_ly_sim_pxs - (hsubstrato + h_pxs - res_v);
z_plot = z_ly_sim_pxs - (hsubstrato + z_ly_sim_pxs(end) - res_v);
z_plot_lap_full = [wrev(z_plot);z_plot]; %Creation of drop profile (z coordinate)
ly_polygon = [x_plot_lap_full z_plot_lap_full];
I = insertShape(I,'FilledPolygon',ly_polygon,'Color',color,'Opacity',opacity); %Insert polygon that corresponds to the drop
if needle_in %Creation of needle
    needle_diameter_pxs = (needle_diameter_mm/10)*scale_cm; %Needle outer diameter in pixels
    I = insertShape(I,'FilledRectangle',[(xc-needle_diameter_pxs/2) 1 needle_diameter_pxs (z_plot(end)-z_ly_sim_pxs(end)/2)],'Color',color,'Opacity',opacity);%Insert rectangle that corresponds to the needle
end
Igray = rgb2gray(I); %Convert image to grayscale
Ibw = imbinarize(Igray); %Binarize image
end

