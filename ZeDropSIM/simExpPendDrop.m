function simExpPendDrop()
%SIMEXPPENDDROP Simulation of quasi-static pendant drop experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of quasi-static pendant drop experiment. Generation of
%   grayscale pendant drop images. Possibility of export 3D drop
%   surface/mesh. The routine requires the user to enter: liquid capillary 
%   constant [cm-2](c), minimum curvature at the apex possible [cm-1] 
%   (bmin), maximum curvature at the apex possbile [cm-1](bmax), needle 
%   diameter [mm] (needleDiam), possibility of simulating images with and 
%   without needle (needleON), image resolution [pxs] (res_h and res_v), 
%   image scale [pxs/mm] (scale), edge width (W_pxs), maximum and minimum 
%   image intensity levels (g1 and g2), initial drop volume [uL] (Vini), 
%   final drop volume (Vend), volume deposition rate [uL/s](Vrate), camera 
%   framerate [fps](framerate) and number of images (Nimages). The needle 
%   position is fixed. The routine also allows the user to choose the 
%   directory where the images are going to be exported. An excel file 
%   (.xlsx) is created containing important properties of the pendant 
%   drops. During the computation of pendant drop profiles errors may occur
%   due to poor choices of input parameters. For certain combinations of 
%   curavture at the apex (bmin and bmax) and needle diameter, the pendant 
%   drop may not be simulated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ PENDANT DROP - SIMULATION OF QUASI-STATIC EXPERIMENT --------------------\n');
    resp = input("- Do you want to continue (Y/N) [Y]? ", "s");
    if isempty(resp) || (resp ~= 'Y' && resp ~= 'N')
        resp = 'Y';
    end
    if resp == 'N'
        break;
    end
    %------------------------------------------------------------------
    %----------------------- INPUT PARAMETERS -------------------------
    %------------------------------------------------------------------
    %- Input parameters
    fprintf('Input parameters: \n');
    fprintf('1) Liquid propertes: \n');
    c = input("- Liquid capillary constant in cm-2 [e.g. water = 13.55]: "); %Cappilary constant in cm-2
    if isempty(c)
        c = 13.55; %Capillary constant of water in cm-2
    end
    bmin = input("- Minimum curvature at the apex in cm-1 [7.2]: "); %Minimum curvature at the apex possible in cm-1
    if isempty(bmin)
        bmin = 7.2; %For water
    end
    while 1
        bmax = input("- Maximum curvature at the apex in cm-1 [9.7]: "); %Maximum curvature at the apex possible in cm-1;
        if isempty(bmax)
            bmax = 9.7; %For water
        end
        if bmax < bmin
            fprintf('Invalid maximum curavture at the apex. Please, enter a lower value. \n');
        else
            break;
        end
    end
    interval_b = 50; %Number of b values created between bmin and bmax for volume search
    tolV = 1e-3; %Tolerance for volume determination in uL
    fprintf('2) Setup configuration propertes: \n');
    needleDiam = input("- Needle diameter in mm [0.9088]: "); %Needle diameter in mm
    if isempty(needleDiam)
        needleDiam = 0.9088;
    end
    r_h_cm = (needleDiam/2)/10; %holder outer radius in cm
    needleON = input("- Needle ON (Y/N)? [Y]: ", "s"); %Generate images with (needle == 'Y') or without needle (needle == 'N')
    if isempty(needleON) || (needleON ~= 'Y' && needleON ~= 'N')
        needleON = 'Y';
    end
    if needleON == 'N'
        needle = 0;
    else
        needle = 1;
    end
    tol = 10^-5; %tolerance to find the intersection points of the drop with the holder radius
    %- Numerical integration options
    S_span=(0:.0001:8); % S_span is the step variable for ode45 solver
    initialCo = [0 1e-100 0 0 0]; %Initial conditions for integration
    %- Determination of Vmin and Vmax
    [x_ly_cm_min,z_ly_cm_min,x_ly_sim_cm_min,z_ly_sim_cm_min,indx_last,indx_second,indx_deq,indx_ds,d_eqmin,h_eqmin,hmin,d_smin,Amin,Vmin] = profile_pendant_drop_ly(bmax,c,r_h_cm,S_span,initialCo,tol);
    [x_ly_cm_max,z_ly_cm_max,x_ly_sim_cm_max,z_ly_sim_cm_max,indx_last,indx_second,indx_deq,indx_ds,d_eqmax,h_eqmax,hmax,d_smax,Amax,Vmax] = profile_pendant_drop_ly(bmin,c,r_h_cm,S_span,initialCo,tol);
    fprintf('3) Image propertes: \n');
    res_h = input("- Horizontal image resolution in pxs [1920]: ");
    if isempty(res_h)
        res_h = 1920;
    end
    res_v = input("- Vertical image resolution in pxs [1080]: ");
    if isempty(res_v)
        res_v = 1080;
    end
    resolution = [res_h res_v];
    z_needle_max = res_v/3; %z coordinate of the needle tip in the image
    scale = input("- Image scale in pxs/mm [185]: ");
    if isempty(scale)
        scale = 185;
    end
    scale_cm = scale*10; %Scale in pixels/cm
    imFmt = input("- Image format [.tif]: ", "s");
    if isempty(imFmt)
        imFmt = ".tif";
    end
    fprintf('4) Illumination propertes: \n');
    W_pxs = input("- Mesasure of edges width in pxs [1.8]: ");
    if isempty(W_pxs)
        W_pxs = 1.8;
    end
    n_normal = 31; %Cabezas et al 2007 %Number of normal edge points
    dist_pxs = W_pxs/2; %Cabezas et al 2007 %Distance between normal edge points in pixels
    while 1
        g1 = input("- Minimum grayscale image intensity level (0 to 255) [60]: ");
        if isempty(g1)
            g1 = 60;
        end
        if g1 > 255
            fprintf('Invalid maximum intensity level. Please, enter a lower value. \n');
        else
            break;
        end
    end
    while 1
        g2 = input("- Maximum grayscale image intensity level (0 to 255) [220]: ");
        if isempty(g2)
            g2 = 220;
        end
        if g2 < g1 || g2 > 255
            fprintf('Invalid maximum intensity level. Please, enter a lower value. \n');
        else
            break;
        end
    end
    fprintf('5) Quasi-static experiment propertes: \n');
    promptVini = strcat("- Initial drop volume in uL (min. ",num2str(Vmin,'%.2f')," uL/ max. ",num2str(Vmax,'%.2f')," uL): ");
    while 1
        Vini = input(promptVini);
        if Vini < Vmin || Vini > Vmax
            fprintf('Invalid initial drop volume. Please, enter another value. \n');
        else
            break;
        end
    end
    promptVend = strcat("- Final drop volume in uL (min. ",num2str(Vmin,'%.2f')," uL/ max. ",num2str(Vmax,'%.2f')," uL): ");
    while 1
        Vend = input(promptVend);
        if Vend < Vmin || Vend > Vmax
            fprintf('Invalid find drop volume. Please, enter another value. \n');
        elseif Vend < Vini
            fprintf('Final drop volume is lower than initial drop volume. Please, enter a higher value. \n');
        else
            break;
        end
    end
    Vrate = input("- Flow rate in uL/s [0.5]: ");
    if isempty(Vrate)
        Vrate = 0.5;
    end
    framerate = input("- Camera framerate in fps [30]: ");
    if isempty(framerate)
        framerate = 30;
    end
    totalframes = floor(framerate*(Vend-Vini)/Vrate); %Total number of frames acquired by the camera during the experiment
    promptnImages = strcat("- Number of images (max. ",num2str(totalframes),") [100]: ");
    while 1
        Nimages = input(promptnImages);
        if isempty(Nimages)
            Nimages = 100;
        end
        if Nimages > totalframes
            fprintf('The total number of frames is lower than the number of images. Please, enter a lower value. \n');
        else
            break;
        end
    end
    V = linspace(Vini,Vend,Nimages); %Vector containing the volume of each drop image
    t_total = abs(Vmax - Vmin)/Vrate; %Total time of the experiment
    time = linspace(0,t_total,Nimages); %Vector containing the time of each drop image
    %- 3D Surface/mesh properties
    fprintf('6) 3D surface/mesh propertes: \n');
    exp3DDropSurf = input("- Do you want to export 3D drop surface/mesh (Y/N) [Y]? ", "s");
    if isempty(exp3DDropSurf) || (exp3DDropSurf ~= 'Y' && exp3DDropSurf ~= 'N')
        exp3DDropSurf = 'Y';
    end
    if exp3DDropSurf == 'Y'
        n = input("- Number of points (mesh) around droplet circunference [100]: ");
        if isempty(n)
            n = 100; %number of points around the cilinder circunference of drop profile
        end
    end
    fprintf('Selection of the directory to export files ... \n');
    save_directory = uigetdir('C:\'); %Save directory

    %- Variables initialization
    b = zeros(length(V),1); %%Curvature at the apex in cm-1
    r0 = zeros(length(V),1); %Radius of curvature at the in cm
    d_eq_cm = zeros(length(V),1); %Equatorial(maximum) diameter in cm
    d_eq_pxs = zeros(length(V),1); %Equatorial(maximum) diameter in pxs
    h_eq_cm = zeros(length(V),1); %Distance from the drop apex to the equatorial (max.) diameter in cm
    h_eq_pxs = zeros(length(V),1); %Distance from the drop apex to the equatorial (max.) diameter in pxs
    d_s_cm = zeros(length(V),1); %%Diameter at a d_eq distance from the drop apex in cm
    d_s_pxs = zeros(length(V),1); %Diameter at a d_eq distance from the drop apex in pxs
    h_cm = zeros(length(V),1); %Distance from the drop apex to the needle in cm
    h_pxs = zeros(length(V),1); %Distance from the drop apex to the needle in pxs
    A_cm2 = zeros(length(V),1); %Drop surface area in cm²
    V_uL = zeros(length(V),1); %Drop volume in uL calculated concomitantly with the Laplace curve generation

    %--------------------------------------------------------------------------
    %------------ SIMULATION OF GRAYSCALE PENDANT DROP IMAGES -----------------
    %--------------------------------------------------------------------------
    fprintf("--------------------------------------------------------------------------------\n")
    fprintf("Simulation of quasi-static pendant drop experiment... \n")
    for k=1:length(V)
        %Determination of pendant drop profile
        [x_ly_sim_cm,z_ly_sim_cm,b(k),r0(k),d_eq_cm(k),h_eq_cm(k),h_cm(k),d_s_cm(k),A_cm2(k),V_uL(k)] = pendant_findbfromV(c,r_h_cm,tol,V(k),bmin,bmax,interval_b,tolV,S_span,initialCo);

        %Conversion to pixels
        x_ly_sim_pxs = x_ly_sim_cm*scale_cm;   % Drop width in pixels
        z_ly_sim_pxs = z_ly_sim_cm*scale_cm;   % Drop height in pixels
        d_eq_pxs(k) = d_eq_cm(k)*scale_cm;
        h_eq_pxs(k) = h_eq_cm(k)*scale_cm;
        d_s_pxs(k) = d_s_cm(k)*scale_cm;
        h_pxs(k) = h_cm(k)*scale_cm;

        %----------------------------------------------------------------------
        %------------------ EXPORT 3D DROP SURFACE (.OBJ) ---------------------
        %----------------------------------------------------------------------
        if exp3DDropSurf == 'Y'
            %Create 3D drop profile
            if needle
                needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
                pos_x_needle = needle_diameter_pxs/2;
                z_needle_max3D = z_ly_sim_pxs(end)+res_v; %maximum z coordinate of the needle
                z_needle = linspace(z_ly_sim_pxs(end),z_needle_max3D,ceil(res_v)); %z coordinates of the needle profile
                x_needle = ones(1,length(z_needle))*pos_x_needle; %x coordinates of the right needle profile
                x_surf_lap = [x_ly_sim_pxs' x_needle];
                z_surf_lap = [z_ly_sim_pxs' z_needle];
            else
                x_surf_lap = x_ly_sim_pxs';
                z_surf_lap = z_ly_sim_pxs';
            end
            [X,Y,Z,gof] = dropsurface(x_surf_lap',z_surf_lap',n); %Generate 3D drop surface
            %figure
            mesh(X,Y,Z)
            axis equal
            %Saving MATLAB figure
            filename_tif = strcat('3D_Surface_Pendant_drop_Needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_scale(pxs_cm)_',num2str(scale_cm),'_Rsquare_',num2str(gof.rsquare,'%.2f'),'.tif');
            fullfilename = fullfile(save_directory,filename_tif); %Object fullfilename
            saveas(gcf,fullfilename); %Saving Matlab figure
            %Export to .obj format
            filename_obj = strcat('3D_Surface_Pendant_drop_Needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_scale(pxs_cm)_',num2str(scale_cm),'_Rsquare_',num2str(gof.rsquare,'%.2f'),'.obj');
            fullfilename = fullfile(save_directory,filename_obj); %Object fullfilename
            saveobjmesh(fullfilename,X,Y,Z); %Exporting 3D drop surface as .obj
        end
        %----------------------------------------------------------------------
        %------------------ GENERATION OF GRAYSCALE IMAGE ---------------------
        %----------------------------------------------------------------------
        %Creation of needle and pendandt drop edge and normal profiles
        res_v = resolution(2);
        res_h = resolution(1);
        xc = round(res_h/2); %x coordinate of the drop center in pixels
        x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %x coordinates of the full drop profile in the final image
        if needle %Construction of pendant drop with needle
            zc = round(res_v/2);
            %z_drop = -z_ly_sim_pxs+z_ly_sim_pxs(end)+zc-(z_ly_sim_pxs(end)/2);  %z coordinates of the left drop profile in the final image with needle
            z_drop = -z_ly_sim_pxs+z_ly_sim_pxs(end)+z_needle_max;  %z coordinates of the left drop profile in the final image with needle
            % Creation of needle drop edge and normal profiles
            needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
            pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
            pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
            %z_needle_max = zc-(z_ly_sim_pxs(end)/2); %maximum z coordinate of the needle
            z_needle = linspace(0,z_needle_max,ceil(z_needle_max)); %z coordinates of the needle profile
            x_needle_l = ones(1,length(z_needle))*pos_x_needle_l; %x coordinates of the left needle profile
            x_needle_r = ones(1,length(z_needle))*pos_x_needle_r; %x coordinates of the right needle profile
            [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge(x_needle_l,z_needle,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left needle profile
            [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge(x_needle_r,z_needle,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right needle profile
            g_needle_r = wrev(g_needle_r); %Fix grayscale order in the right profile of the needle
            xn_needle = [xn_needle_l xn_needle_r]; %x coordinates of edge and normal points of the needle profile
            zn_needle = [zn_needle_l zn_needle_r]; %z coordinates of edge and normal points of the needle profile
            g_needle = [g_needle_l g_needle_r]; %grayscale value of edge and normal points of the needle profile
        else %Construction of pendant drop without needle
            z_drop = -z_ly_sim_pxs + z_ly_sim_pxs(end); %z coordinates of the left drop profile in the final image without needle
        end
        z_drop_lap_full = [wrev(z_drop);z_drop]; %z coordinates of the full drop profile in the final image)
        %[xn,zn,g] = grayscale_edge(x_ly_sim_cm,z_ly_sim_cm,n_normal,g1,g2,W_lx);
        [xn_drop,zn_drop,g_drop] = grayscale_edge(x_drop_lap_full,z_drop_lap_full,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the drop profile
        %- Combine needle and drop edge and normal points
        if needle
            xn = [xn_needle xn_drop];
            zn = [zn_needle zn_drop];
            g = [g_needle g_drop];
        else
            xn = [xn_drop];
            zn = [zn_drop];
            g = [g_drop];
        end
        %Generation of synthetical illuminated (grayscale) pendant drop images
        M = zeros(res_v,res_h); %Array with image resolution size
        I = uint8(M); %Image creation
        %- Plotting needle and drop profile edge
        for i=1:length(g)
            if (xn(i)>=0) && (xn(i)<res_h) && (zn(i)>=0) && (zn(i)<res_v)
                if double(I(floor(zn(i))+1,floor(xn(i))+1)) == 0
                    I(floor(zn(i))+1,floor(xn(i))+1) = round(g(i));
                else
                    I(floor(zn(i))+1,floor(xn(i))+1) = round((double(I(floor(zn(i))+1,floor(xn(i))+1)) + g(i))/2);
                end
            end
        end
        %- Plotting background
        for i = 1:res_v %Search from left to right
            for j = 1:res_h
                if I(i,j) == 0
                    I(i,j) = g2;
                else
                    break;
                end
            end
        end
        for i = 1:res_v %Search from right to left
            for j = res_h:-1:1
                if I(i,j) == 0
                    I(i,j) = g2;
                else
                    break;
                end
            end
        end
        %- Filling drop and needle
        for i = 1:res_v %Search in all the image
            for j = 1:res_h
                if I(i,j) == 0
                    I(i,j) = g1;
                end
            end
        end

        %Save image
        imagename = strcat('Grayscale_Pendant_drop_Needle_',num2str(needle),'_c(cm-2)_',num2str(c),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
        fullimagename = strcat(imagename,imFmt); %Image name + image format
        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
        imwrite(I,fullfilename) %Save Image
    end

    %- Write excel file with important pendant drop profile properties
    varnames = {'time[s]','V[uL]','b[cm-1]','r0[cm]','d_eq[cm]','d_eq[pxs]','h_eq[cm]','h_eq[pxs]','d_s[cm]','d_s[pxs]','h[cm]','h[pxs]','A[cm2]'}; %Variables name in the excel file
    T = table(time',V_uL,b,r0,d_eq_cm,d_eq_pxs,h_eq_cm,h_eq_pxs,d_s_cm,d_s_pxs,h_cm,h_pxs,A_cm2,'VariableNames',varnames); %Table creation
    sheet_name = strcat('c_',num2str(c,'%.2f'),'_cm-2','_rh_',num2str(r_h_cm,'%.2f'),'_cm'); %Excel sheet name
    savename_xls = strcat('Dimensions_PendantDrop_c(cm-2)_',num2str(c),'_rh(cm)_',num2str(r_h_cm),'_Vini(uL)_',num2str(Vini),'_Vend(uL)_',num2str(Vend),'_Vrate(uL_s)_',num2str(Vrate),'.xlsx'); %Excel file name
    fullsavename_xls = fullfile(save_directory,savename_xls);
    writetable(T,fullsavename_xls,'Sheet',sheet_name) %Write table to an excel file

    fprintf("Simulation of quasi-static experiment completed! \n")
    pause(2)
end
end