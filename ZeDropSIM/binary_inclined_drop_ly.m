function [Ibw] = binary_inclined_drop_ly(x_ly_sim_cm,z_ly_sim_cm,hsubstrato,resolution,scale_cm,inclined,alpha_deg)
%BINARY_INCLINED_DROP_LY %Creation of a synthetical binary inclined drop image
%   INPUT:
% x_ly_sim_cm - x coordinates of drop profile in cm
% z_ly_sim_cm - z coordinates of drop profile in cm
% hsubstrato - height of the solid substrate in pixels
% resolution - image resolution [res_h res_v]
% scale_cm - image scale in cm
% inclined - show drop and the baseline inclined (=1)

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
zc = round(res_v/2); %z coordinate of the image center in pixels
M = ones(res_v,res_h)*255; %Array with image resolution size
I = uint8(M); %Image creation
x_plot = x_ly_sim_pxs + xc - (x_ly_sim_pxs(1)+(x_ly_sim_pxs(end)-x_ly_sim_pxs(1))/2);
if inclined
    %Rotating inclined drop profile
    z_plot = z_ly_sim_pxs + res_v/2 - (max(abs(z_ly_sim_pxs))/2);
    [x_plot_rot,z_plot_rot] = rotate_data(x_plot,z_plot,xc,zc,alpha_deg);
    ly_polygon = [x_plot_rot' z_plot_rot'];
    %Generating and rotating solid substrate polygon 
    z_subst_sup = z_plot(end);
    z_subst_inf = z_plot(end) + 2*res_h;
    x_subst_l = -res_h/2;
    x_subst_r = 3*res_h/2;
    [x_subst_l_sup_rot, z_subst_l_sup_rot] = rotate_data(x_subst_l,z_subst_sup,xc,zc,alpha_deg);
    [x_subst_r_sup_rot, z_subst_r_sup_rot] = rotate_data(x_subst_r,z_subst_sup,xc,zc,alpha_deg);
    [x_subst_r_inf_rot, z_subst_r_inf_rot] = rotate_data(x_subst_r,z_subst_inf,xc,zc,alpha_deg);
    [x_subst_l_inf_rot, z_subst_l_inf_rot] = rotate_data(x_subst_l,z_subst_inf,xc,zc,alpha_deg);
    x_subst_rot = [x_subst_l_sup_rot x_subst_r_sup_rot x_subst_r_inf_rot x_subst_l_inf_rot];
    z_subst_rot = [z_subst_l_sup_rot z_subst_r_sup_rot z_subst_r_inf_rot z_subst_l_inf_rot];
    subst_polygon = [x_subst_rot' z_subst_rot'];
    I = insertShape(I,'FilledPolygon',subst_polygon,'Color',color,'Opacity',opacity); %Insert polygon that corresponds to the solid substrate
else
    z_plot = z_ly_sim_pxs + res_v - (max(abs(z_ly_sim_pxs))+hsubstrato);
    ly_polygon = [x_plot z_plot];
    I = insertShape(I,'FilledRectangle',[1 (res_v-hsubstrato) res_h hsubstrato],'Color',color,'Opacity',opacity); %Insert rectangle that corresponds to the solid substrate
end
I = insertShape(I,'FilledPolygon',ly_polygon,'Color',color,'Opacity',opacity); %Insert polygon that corresponds to the drop
Igray = rgb2gray(I); %Convert image to grayscale
Ibw = imbinarize(Igray); %Binarize image
%{
plot(x_plot_rot,z_plot_rot)
xlim([0 1920])
ylim([0 1080])
hold on
plot(x_subst_rot,z_subst_rot)
%}
end

