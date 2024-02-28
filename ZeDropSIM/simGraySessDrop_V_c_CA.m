function simGraySessDrop_V_c_CA()
%SIMGRAYSESSDROP_V_C_CA Simulation of grayscale sessile drop images for 
% different drop volume (V) given liquid capillary constant (c) and contact
% angle (CA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of grayscale sessile drop images for different drop volumes 
%   (V) given liquid capillary constant (c) and contact angle (CA). Some 
%   source of  errors, such as edge width, camera vertical misalignment, 
%   random perturbation in the drop profile, lack of contrast, non-uniform 
%   illumination, random noise, Gaussian noise, impulse salt-pepper noise
%   and satellite droplets, can be simulated. The 3D drop surface can also 
%   be exported. Depending on the simulated option, the program requires 
%   the user to enter some specific parameters. As common input parameters,
%   any option requires the user to enter: the minimum drop volume [uL]
%   (Vmin), the maximum drop volume [uL](Vmax), the number of drop volume 
%   levels that is going to be simulated between Vmin and Vmax (nVol), the 
%   contact angle [°](theta), the needle diameter [mm](needleDiam), the 
%   possibility of simulating images with and without needle (needleON), 
%   the image resolution [pixels](res_h and res_v), the image scale 
%   [pixels/mm](scale), edge width (W_pxs), maximum and minimum image 
%   intensity levels (g1 and g2). In order to make the execution of the 
%   algorithm faster, it allows the user to load an excel file containing 
%   important drop properties such as curvature at the apex (b), drop 
%   volume (V) and contact angle (CA). The height of the solid substrate 
%   in the image is fixed. The routine also allows the user to choose the 
%   directory where the images are going to be exported. An excel file 
%   (.xlsx) is created containing important properties of the simulated
%   sessile drops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ SESSILE DROP - SIMULATION OF GRAYSCALE IMAGES --------------------\n');
    % Options
    fprintf('OPTIONS: \n');
    fprintf('1. Generation of grayscale images; \n');
    fprintf('2. Generation of 3D drop surface; \n');
    fprintf('3. Simulation of blur in liqud drop and solid substrate profiles (edge width); \n');
    fprintf('4. Simulation of camera vertical misalignment; \n');
    fprintf('5. Simulation of random perturbations in the drop profile; \n')
    fprintf('6. Simulation of lack of contrast; \n')
    fprintf('7. Simulation of non-uniform ilumination; \n')
    fprintf('8. Simulation of random noise; \n')
    fprintf('9. Simulation of Gaussian noise; \n')
    fprintf('10. Simulation of impulse salt-pepper noise; \n')
    fprintf('11. Simulation of satellite droplets; \n')

    option = input('Enter one of the options above (0 to back to main menu) [0]: ');
    if isempty(option)
        option = 0;
    end
    if option == 0
        break;
    elseif option > 0 && option <= 11
        %------------------------------------------------------------------
        %----------------------- INPUT PARAMETERS -------------------------
        %------------------------------------------------------------------
        %- COMMON INPUT PARAMETERS
        fprintf('Input parameters: \n');
        fprintf('1) Liquid propertes: \n');
        c = input("- Liquid capillary constant in cm-2 [e.g. water = 13.55]: "); %Cappilary constant in cm-2
        if isempty(c)
            c = 13.55; %Capillary constant of water in cm-2
        end
        loadExc = input('- Load excel file containing curvature at the apex (b), drop volume (V) and contact angle (CA) values (Y/N) [N]: ', 's');
        if isempty(loadExc)
            loadExc = 'N';
        end
        if loadExc == 'Y'
            [openname_xls,path] = uigetfile('*.xlsx'); %Get file path to open excel file
            fullopenname_xls = fullfile(path,openname_xls);
            readT = readtable(fullopenname_xls,"VariableNamingRule","preserve"); %Read excel file
            b = readT{:,1};
            %r0 = 1./b; %radius at the apex in cm
            theta = readT{:,3};
            V = readT{:,4};
            nSimImages = length(b); %Number of simulated images
        else
            bmin_ini = 0.05; %Minimum curvature at the apex that is going to be searched in cm-1
            bmax_ini = 10; %Maximum curvature at the apex that is going to be searched in cm-1
            interval_b = 50; %Number of b values created between bmin and bmax for drop volume search
            tolV = 1e-3; %Tolerance for volume determination in uL
            itermax = 5; %Maximum number of iterations for volume search
            Vmin = input("- Minimum drop volume in uL [5]: "); %Minimum drop volume in uL
            if isempty(Vmin)
                Vmin = 5;
            end
            while 1
                Vmax = input("- Maximum drop volume in uL [30]: "); %Maximum drop volume in uL
                if isempty(Vmax)
                    Vmax = 30;
                end
                if Vmax < Vmin
                    fprintf('Invalid maximum drop volume. Please, set a lower value. \n');
                else
                    break;
                end
            end
            nVol = input("- Number of drop volume levels desired [3]: "); %Number of simulated images between minimum and maximum drop volume
            if isempty(nVol)
                nVol = 3;
            end
            V = linspace(Vmin,Vmax,nVol);
            theta = input("- Contact angle in degree [70]: "); %Contact angle in degrees
            if isempty(theta)
                theta = 70;
            end
            nSimImages = length(V); %Number of simulated images
            b = zeros(nSimImages,1);
        end
        fprintf('2) Drop configuration propertes: \n');
        needleDiam = input("- Needle diameter in mm [0.9088]: "); %Needle diameter in mm
        if isempty(needleDiam)
            needleDiam = 0.9088;
        end
        r_h_cm = (needleDiam/2)/10; %Holder outer radius in cm
        needleON = input("- Needle ON (Y/N)? [Y]: ", "s"); %Generate images with (needle == 'Y') or without needle (needle == 'N')
        if isempty(needleON) || (needleON ~= 'Y' && needleON ~= 'N')
            needleON = 'Y';
        end
        if needleON == 'N'
            needle_in = 0;
        else
            needle_in = 1;
        end
        hsubstrato = 100; %Substrate height in pixels;
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
        %Numerical integration options
        S_span = (0:.0001:8); % S_span is the step variable for ode45 solver
        initialCo = [0 1e-100 0 0 0]; %Initial conditions for integration

        %- SPECIFIC INPUT PARAMETERS
        if option == 2 %Generation of 3D drop surface
            %- 3D Surface/mesh properties
            fprintf('5) 3D Surface/mesh properties: \n');
            n = input("- Number of points (mesh) around droplet circunference [100]: ");
            if isempty(n)
                n = 100; %number of points around the cilinder circunference of drop profile
            end
        elseif option == 3 %Simulation of edge width
            %- Edge width properties
            fprintf('5) Edge width properties: \n');
            W_pxs_min = input("- Minimum edge width in pxs [0]: ");
            if isempty(W_pxs_min)
                W_pxs_min = 0;
            end
            while 1
                W_pxs_max = input("- Maximum edge width in pxs [2]: ");
                if isempty(W_pxs_max)
                    W_pxs_max = 2;
                end
                if W_pxs_max < W_pxs_min
                    fprintf('Invalid maximum edge width. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            nEdgeWidth = input("- Number of edge width simulated [11]: ");
            if isempty(nEdgeWidth)
                nEdgeWidth = 11;
            end
            W_pxs_array = linspace(W_pxs_min,W_pxs_max,nEdgeWidth);
        elseif option == 4 %Simulation of camera vertical misalignment
            %- Camera vertical misalignment properties
            fprintf('5) Camera vertical misalignment properties: \n');
            phi_min = input("- Minimum camera vertical misalignment in degree [-2]: ");
            if isempty(phi_min)
                phi_min = -2;
            end
            while 1
                phi_max = input("- Maximum camera vertical misalignment in degree [2]: ");
                if isempty(phi_max)
                    phi_max = 2;
                end
                if phi_max < phi_min
                    fprintf('Invalid maximum camera vertical misalignment. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            nVertMisalign = input("- Number of vertical misalignment simulated [11]: ");
            if isempty(nVertMisalign)
                nVertMisalign = 11;
            end
            phi = linspace(phi_min,phi_max,nVertMisalign);
            xc_rot = round(res_h/2); %x coordinate of the center of rotation
            zc_rot = round(res_v/2); %z coordinate of the center of rotation
        elseif option == 5 %Simulation of random perturbations in the drop profile
            %- Random perturbations properties
            fprintf('5) Random perturbation in drop profile properties: \n');
            pert_min = input("- Minimum normal perturbation in pixels [0]: ");
            if isempty(pert_min)
                pert_min = 0;
            end
            while 1
                pert_max = input("- Maximum normal perturbation in pixels [5]: ");
                if isempty(pert_max)
                    pert_max = 5;
                end
                if pert_max < pert_min
                    fprintf('Invalid maximum normal perturbation. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            nRandPert = input("- Number of random perturbations simulated [11]: ");
            if isempty(nRandPert)
                nRandPert = 11;
            end
            pert = linspace(pert_min,pert_max,nRandPert);
            n_images = input("- Number of images simulated for each level [5]: ");
            if isempty(n_images)
                n_images = 5;
            end
        elseif option == 6 %Simulation of lack of contrast
            %- Lack of contrast properties
            fprintf('5) Lack of contrast properties: \n');
            while 1
                diff_g_min = input("- Minimum difference between min.and max.image intensity levels [20]: ");
                if isempty(diff_g_min)
                    diff_g_min = 20;
                end
                if diff_g_min > 255
                    fprintf('Invalid minimum difference. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            while 1
                diff_g_max = input("- Maximum difference between min. and max. image intensity levels [190]: ");
                if isempty(diff_g_max)
                    diff_g_max = 190;
                end
                if diff_g_max > 255 || diff_g_max < diff_g_min
                    fprintf('Invalid maximum difference. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            nLackCont = input("- Number of lack of contrast simulated [11]: ");
            if isempty(nLackCont)
                nLackCont = 11;
            end
            diff_g = linspace(diff_g_min,diff_g_max,nLackCont);
        elseif option == 7 %Simulation of non-uniform ilumination
            %- Non-uniform illumination properties
            fprintf('5) Non-uniform illumination properties: \n');
            while 1
                deltag_min = input("- Minimum variation of image intensity [0]: ");
                if isempty(deltag_min)
                    deltag_min = 0;
                end
                if deltag_min > 255
                    fprintf('Invalid minimum variation. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            while 1
                deltag_max = input("- Maximum variation of image intensity [60]: ");
                if isempty(deltag_max)
                    deltag_max = 60;
                end
                if deltag_max < deltag_min
                    fprintf('Invalid maximum variation of image intensity. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            nNonUniIll = input("- Number of non-uniform illumination levels simulated [11]: ");
            if isempty(nNonUniIll)
                nNonUniIll = 11;
            end
            deltag = linspace(0,deltag_max,nNonUniIll);
        elseif option == 8 %Simulation of random noise
            %- Random noise properties
            fprintf('5) Random noise propertes: \n');
            noise_max = 0.4*(g2-g1); %Maximum level of random noise
            nRandNoise = input("- Number of random noise levels simulated [11]: ");
            if isempty(nRandNoise)
                nRandNoise = 11;
            end
            noise = linspace(0,noise_max,nRandNoise);
            n_images = input("- Number of images simulated for each level [5]: ");
            if isempty(n_images)
                n_images = 5;
            end
        elseif option == 9 %Simulation of Gaussian noise
            %- Gaussian noise properties
            fprintf('5) Gaussian noise propertes: \n');
            g_noise_mean = input("- Mean of Gaussian noise [0]: ");
            if isempty(g_noise_mean)
                g_noise_mean = 0;
            end
            g_noise_var_min = input("- Minimum variance of Gaussian noise [0]: ");
            if isempty(g_noise_var_min)
                g_noise_var_min = 0;
            end
            while 1
                g_noise_var_max = input("- Maximum variance of Gaussian noise [0.05]: ");
                if isempty(g_noise_var_max)
                    g_noise_var_max = 0.05;
                end
                if g_noise_var_max < g_noise_var_min
                    fprintf('Invalid maximum variance of Gaussian noise. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            nGaussNoise = input("- Number of Gaussian noise levels simulated [11]: ");
            if isempty(nGaussNoise)
                nGaussNoise = 11;
            end
            g_noise_var = linspace(g_noise_var_min,g_noise_var_max,nGaussNoise);
            n_images = input("- Number of images simulated for each level [5]: ");
            if isempty(n_images)
                n_images = 5;
            end
        elseif option == 10 %Simulation of impulse salt-pepper noise
            %- Impulse salt-pepper noise properties
            fprintf('5) Impulse salt-pepper noise propertes: \n');
            sp_noise_density_min = input("- Minimum salt-pepper density [0]: ");
            if isempty(sp_noise_density_min)
                sp_noise_density_min = 0;
            end
            while 1
                sp_noise_density_max = input("- Maximum salt-pepper density [0.2]: ");
                if isempty(sp_noise_density_max)
                    sp_noise_density_max = 0.2;
                end
                if sp_noise_density_max < sp_noise_density_min
                    fprintf('Invalid maximum salt-pepper density. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            nSaltPepNoise = input("- Number of impulse salt-pepper noise levels simulated [11]: ");
            if isempty(nSaltPepNoise)
                nSaltPepNoise = 11;
            end
            sp_noise_density = linspace(sp_noise_density_min,sp_noise_density_max,nSaltPepNoise);
            n_images = input("- Number of images simulated for each level [5]: ");
            if isempty(n_images)
                n_images = 5;
            end
        elseif option == 11 %Simulation of satellite droplets
            %- Satellite droplets properties
            fprintf('5) Satellite droplets propertes: \n');
            d_sat_cm = input("- Diameter of satellite droplets in cm [0.006]: ");
            if isempty(d_sat_cm)
                d_sat_cm = 0.006;
            end
            d_sat_pxs = d_sat_cm*scale_cm; %Diameter of satellite droplets in pixels
            a_min = input("- Minimum density of satellite droplets [0]: ");
            if isempty(a_min)
                a_min = 0;
            end
            while 1
                a_max = input("- Maximum density of satellite droplets [0.06]: ");
                if isempty(a_max)
                    a_max = 0.06;
                end
                if a_max < a_min
                    fprintf('Invalid maximum density of satellite droplets. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            nSatDrop = input("- Number of satellite droplets density levels simulated [11]: ");
            if isempty(nSatDrop)
                nSatDrop = 11;
            end
            a = linspace(a_min,a_max,nSatDrop);
            n_images = input("- Number of images simulated for each level [5]: ");
            if isempty(n_images)
                n_images = 5;
            end
        end
        %- Save options
        fprintf('Selection of the directory to export files ... \n');
        save_directory = uigetdir('C:\'); %Save directory
   
        %- Variables initialization
        theta_sim = zeros(nSimImages,1); %Contact angle (theta) effectively used for drop simulation
        r_w_cm = zeros(nSimImages,1); %Wetting radius in cm
        r_w_pxs = zeros(nSimImages,1); %Wetting radius in pxs
        r_eq_cm = zeros(nSimImages,1); %rEquatorial radius (maximum) of the drop in cm
        r_eq_pxs = zeros(nSimImages,1); %Equatorial radius (maximum) of the drop in pxs
        h_cm = zeros(nSimImages,1); %Drop height in cm
        h_pxs = zeros(nSimImages,1); %Drop height in pxs
        A_cm2 = zeros(nSimImages,1); %Drop surface area in cm²
        V_uL = zeros(nSimImages,1); %Drop volume in uL calculated concomitantly with the Laplace curve generation
        xCPleft_pxs = zeros(nSimImages,1); %x coordinate of the left triple contact point
        yCPleft_pxs = zeros(nSimImages,1); %y coordinate of the left triple contact point
        xCPright_pxs = zeros(nSimImages,1); %x coordinate of the right triple contact point
        yCPright_pxs = zeros(nSimImages,1); %y coordinate of the right triple contact point
        
    %------------------------------------------------------------------------
    %----------------- SIMULATION OF SESSILE DROP IMAGES --------------------
    %------------------------------------------------------------------------
    %- Generation of grayscale sessile drop images
    %fprintf('Generating binary sessile drop images... \n');
    for k = 1:nSimImages
        fprintf('----------------------------------------------------------------------------------------- \n');     
        fprintf('Simulating sessile drop profile for c = %.2f cm-2, CA = %.2f° and V = %.2f uL ... \n',c,mean(theta),V(k));
        if loadExc == 'Y'   
            %Determination of sessile drop profile by numerical integration of the Young-Laplace equation - given c, b and CA values
            [x_ly_sim_cm,z_ly_sim_cm,theta_sim(k),r_w_cm(k),r_eq_cm(k),h_cm(k),A_cm2(k),V_uL(k)] = profile_sessile_drop_ly(b(k),c,theta(k),S_span,initialCo);
        else
            %Determination of sessile drop profile by numerical integration of the Young-Laplace equation - given c, V and CA values
            [x_ly_sim_cm,z_ly_sim_cm,theta_sim(k),b(k),r_w_cm(k),r_eq_cm(k),h_cm(k),A_cm2(k),V_uL(k)] = profile_sessile_drop_ly_V(V(k),c,theta,S_span,initialCo,bmin_ini,bmax_ini,interval_b,tolV,itermax);
        end
        
        %Conversion to pixels
        x_ly_sim_pxs = x_ly_sim_cm*scale_cm;   % Drop width in pixels
        z_ly_sim_pxs = z_ly_sim_cm*scale_cm;   % Drop height in pixels
        %S_ly_sim_pxs = S_ly_sim_cm*scale_cm;   % Arc length in pixels
        %Construction of synthetical image
        res_v = resolution(2);
        res_h = resolution(1);
        xc = round(res_h/2); %x coordinate of the drop center in pixels

        %--------------------------------------------------------------
        %-------------- GENERATION OF GRAYSCALE IMAGES ----------------
        %--------------------------------------------------------------
        if option == 1 || option == 8 || option == 9 || option == 10 || option == 11 % Generation of grayscale images
            %Creation of sessile drop and needle edge normal profile
            x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)
            z_drop = z_ly_sim_pxs - (hsubstrato + z_ly_sim_pxs(end) - res_v);
            z_drop_lap_full = [wrev(z_drop);z_drop]; %Creation of drop profile (z coordinate)
            xCPleft_pxs(k) = x_drop_lap_full(1);
            xCPright_pxs(k) = x_drop_lap_full(end);
            yCPleft_pxs(k) = z_drop_lap_full(1);
            yCPright_pxs(k) = z_drop_lap_full(end);

            if needle_in
                %Construction of the needle profile
                needle_diameter_pxs = (needleDiam/10)*scale_cm;
                pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                indx_l = find((abs(x_drop_lap_full-pos_x_needle_l)) == min(abs(x_drop_lap_full-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                indx_r = find((abs(x_drop_lap_full-pos_x_needle_r)) == min(abs(x_drop_lap_full-pos_x_needle_r)),1,'last'); %Find the index of the instersection between the right side of the needle and the drop
                z_needle_max_l = z_drop_lap_full(indx_l);
                z_needle_max_r = z_drop_lap_full(indx_r);
                %z_needle_max = z_drop(end)-z_ly_sim_pxs(end)/2;
                z_needle_l = linspace(0,z_needle_max_l,ceil(z_needle_max_l)); %z coordinates of the needle profile
                z_needle_r = linspace(0,z_needle_max_r,ceil(z_needle_max_r)); %z coordinates of the needle profile
                %z_needle = linspace(0,z_needle_max,ceil(z_needle_max)); %z coordinates of the needle profile
                x_needle_l = ones(1,length(z_needle_l))*pos_x_needle_l; %x coordinates of the left needle profile
                x_needle_r = ones(1,length(z_needle_r))*pos_x_needle_r; %x coordinates of the right needle profile
                [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge(x_needle_l,z_needle_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left needle profile
                [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge(x_needle_r,z_needle_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right needle profile
                g_needle_r = wrev(g_needle_r); %Fix grayscale order in the right profile of the needle
                xn_needle = [xn_needle_l xn_needle_r]; %x coordinates of edge and normal points of the needle profile
                zn_needle = [zn_needle_l zn_needle_r]; %z coordinates of edge and normal points of the needle profile
                g_needle = [g_needle_l g_needle_r]; %grayscale value of edge and normal points of the needle profile

                %Trim drop profile based on the needle profile
                x_drop_l = x_drop_lap_full(1:indx_l);
                z_drop_l = z_drop_lap_full(1:indx_l);
                x_drop_r = x_drop_lap_full(indx_r:end);
                z_drop_r = z_drop_lap_full(indx_r:end);
                [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the drop
                g_drop_l = wrev(g_drop_l);
                [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the drop
                g_drop_r = wrev(g_drop_r);
                xn_drop = [xn_drop_l xn_drop_r];
                zn_drop = [zn_drop_l zn_drop_r];
                g_drop = [g_drop_l g_drop_r];
            else
                [xn_drop,zn_drop,g_drop] = grayscale_edge(x_drop_lap_full,z_drop_lap_full,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the drop profile
                g_drop = wrev(g_drop);
            end

            %Creation of solid subratre edge normal profile
            x_sol_l = linspace(0,x_drop_lap_full(1),ceil(x_drop_lap_full(1)));
            z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
            x_sol_r = linspace(x_drop_lap_full(end),res_h,ceil(res_h-x_drop_lap_full(end)));
            z_sol_r = ones(1,length(x_sol_r))*(res_v-hsubstrato);
            [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the solid substrate surface
            g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
            [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the solid substrate surface
            g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
            xn_sol = [xn_sol_l xn_sol_r];
            zn_sol = [zn_sol_l zn_sol_r];
            g_sol = [g_sol_l g_sol_r];

            %Trim drop normal edge according to solid substrate
            ind_z = find(zn_drop <= (res_v-hsubstrato));
            xn_drop = xn_drop(ind_z);
            zn_drop = zn_drop(ind_z);
            g_drop = g_drop(ind_z);

            %Combine solid, drop and needle edge and normal points
            if needle_in
                xn = [xn_sol xn_needle xn_drop];
                zn = [zn_sol zn_needle zn_drop];
                g = [g_sol g_needle g_drop];
            else
                xn = [xn_sol xn_drop];
                zn = [zn_sol zn_drop];
                g = [g_sol g_drop];
            end
            %Generation of synthetical illuminated (grayscale) pendant drop images
            M = zeros(res_v,res_h); %Array with image resolution size
            I = uint8(M); %Image creation
            %- Plotting solid substrate, drop profile and needle edge
            for i=1:length(g)
                if (xn(i)>=0) && (xn(i)<res_h) && (zn(i)>=0) && (zn(i)<res_v)
                    if double(I(floor(zn(i))+1,floor(xn(i))+1)) == 0
                        I(floor(zn(i))+1,floor(xn(i))+1) = round(g(i));
                    else
                        I(floor(zn(i))+1,floor(xn(i))+1) = round((double(I(floor(zn(i))+1,floor(xn(i))+1)) + g(i))/2);
                    end
                end
            end
            %- Filling background
            for i = 1:(res_v-hsubstrato) %Search from left to right
                for j = 1:res_h
                    if I(i,j) == 0
                        I(i,j) = g2;
                    else
                        break;
                    end
                end
            end
            for i = 1:(res_v-hsubstrato) %Search from right to left
                for j = res_h:-1:1
                    if I(i,j) == 0
                        I(i,j) = g2;
                    else
                        break;
                    end
                end
            end
            %- Filling solid substrate,drop and needle
            for i = 1:res_v %Search in all the image
                for j = 1:res_h
                    if I(i,j) == 0
                        I(i,j) = g1;
                    end
                end
            end
            if option == 1
                %Save images
                imagename = strcat('Grayscale_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_theta(°)_',num2str(theta_sim(k),'%.1f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
                fullimagename = strcat(imagename,imFmt); %Image name + image format
                fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                imwrite(I,fullfilename) %Save Image
            %----------------------------------------------------------
            %------------------ RANDOM NOISE ADDITION -----------------
            %----------------------------------------------------------
            elseif option == 8 %Simulation of random noise
                %- Addition of random noise
                for q = 1:length(noise)
                    for p = 1:n_images
                        Inoise = I;
                        for m = 1:res_v
                            for n = 1:res_h
                                Inoise(m,n) = Inoise(m,n) + randi([round(-noise(q)) round(noise(q))],1);
                            end
                        end
                        %Save noise images
                        imagename = strcat('Noise_Grayscale_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_theta(°)_',num2str(theta_sim(k),'%.1f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_noise_',num2str(round(noise(q))),'_n_images_',num2str(p)); %Image filename
                        fullimagename = strcat(imagename,imFmt); %Image name + image format
                        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                        imwrite(Inoise,fullfilename) %Save Image
                    end
                end
            %----------------------------------------------------------
            %------------------ GAUSSIAN NOISE ADDITION ---------------
            %----------------------------------------------------------
            elseif option == 9 %Simulation of Gaussian noise
                %- Addition of gaussian noise
                for q = 1:length(g_noise_var)
                    for p = 1:n_images
                        Inoise = imnoise(I,'gaussian',g_noise_mean,g_noise_var(q));
                        %Save noise images
                        imagename = strcat('Guassian_Noise_Grayscale_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_theta(°)_',num2str(theta_sim(k),'%.1f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_g_noise_mean_',num2str(g_noise_mean),'_g_noise_var_',num2str(g_noise_var(q),'%.3f'),'_n_images_',num2str(p)); %Image filename
                        fullimagename = strcat(imagename,imFmt); %Image name + image format
                        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                        imwrite(Inoise,fullfilename) %Save Image
                    end
                end
            %----------------------------------------------------------
            %---------- IMPULSE SALT-PEPPER NOISE ADDITION ----------
            %----------------------------------------------------------
            elseif option == 10 %Addition of impulse salt-pepper noise
                %- Addition of impulse salt-pepper noise
                for q = 1:length(sp_noise_density)
                    for p = 1:n_images
                        Inoise = imnoise(I,'salt & pepper',sp_noise_density(q));
                        %Save noise images
                        imagename = strcat('Salt_Pepper_Noise_Grayscale_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_theta(°)_',num2str(theta_sim(k),'%.1f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_salt_pepper_density_',num2str(sp_noise_density(q),'%.2f'),'_n_images_',num2str(p)); %Image filename
                        fullimagename = strcat(imagename,imFmt); %Image name + image format
                        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                        imwrite(Inoise,fullfilename) %Save Image
                    end
                end
            %----------------------------------------------------------
            %-------------- SATELLITE DROPLETS ADDITION ---------------
            %----------------------------------------------------------
            elseif option == 11 %Addition of satellite droplets
                %- Addition of satellite droples
                for p=1:length(a)
                    N = round((4*a(p)*res_h*res_v)/(pi*(d_sat_pxs^2))); %Number of satellite droplets
                    fprintf('Running a = %.3f and N = %d ... \n',a(p),N);
                    for m = 1:n_images
                        Isat = I;
                        for q = 1:N
                            [columnsInImage rowsInImage] = meshgrid(1:res_h, 1:res_v);
                            x_sat_center = randi([1 res_h],1);
                            z_sat_center = randi([1 res_v],1);
                            circlePixels = (rowsInImage-z_sat_center).^2 + (columnsInImage-x_sat_center).^2 <= (d_sat_pxs/2).^2;
                            indx_circle = find(circlePixels);
                            Isat(indx_circle) = g1; %Plotiing satellite droplets on the image
                        end
                        %Save satellite droplet images
                        imagename = strcat('Satellite_Grayscale_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_theta(°)_',num2str(theta_sim(k),'%.1f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_d_sat(cm)_',num2str(d_sat_cm,'%.3f'),'_a_',num2str(a(p),'%.4f'),'_n_images_',num2str(m)); %Image filename
                        fullimagename = strcat(imagename,imFmt); %Image name + image format
                        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                        imwrite(Isat,fullfilename) %Save Image
                    end
                end
            end
        %--------------------------------------------------------------
        %--------------- GENERATION OF 3D DROP SURFACE ----------------
        %--------------------------------------------------------------
        elseif option == 2 % GENERATION OF 3D DROP SURFACE
            %Create 3D drop profile
            x_surf_lap =  x_ly_sim_pxs';
            z_surf_lap = -z_ly_sim_pxs';
            if needle_in
                %Construction of the needle profile
                needle_diameter_pxs = (needleDiam/10)*scale_cm;
                pos_x_needle = round(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                %if theta(j)>90
                %indx_needle = find((abs(x_surf_lap(round(length(x_surf_lap)/2:end))-pos_x_needle)) == min(abs(x_surf_lap(round(length(x_surf_lap)/2:end))-pos_x_needle)),1,'first'); %Find the index of the instersection between the drop profile and the needle
                %indx_needle = length(x_surf_lap)/2 + indx_needle;
                %else
                indx_needle = find((abs(x_surf_lap'-pos_x_needle)) == min(abs(x_surf_lap'-pos_x_needle)),1,'first'); %Find the index of the instersection between the drop profile and the needle
                %end
                z_needle_max = z_surf_lap(indx_needle);
                z_needle = linspace(z_needle_max,z_needle_max+res_v,ceil(res_v)); %z coordinates of the needle profile
                x_needle = ones(1,length(z_needle))*pos_x_needle; %x coordinates of the needle profile
                %Trim drop profile based on the needle profile
                x_drop_lap = x_surf_lap(indx_needle:end);
                z_drop_lap = z_surf_lap(indx_needle:end);
                x_surf_lap = [wrev(x_drop_lap) x_needle];
                z_surf_lap = [wrev(z_drop_lap) z_needle] - z_drop_lap(end);
            end
            [X,Y,Z,gof] = dropsurface(x_surf_lap',z_surf_lap',n); %Generate 3D drop surface
            mesh(X,Y,Z)
            axis equal
            %Saving MATLAB figure
            filename_tif = strcat('3D_Surface_Sessile_drop_Needle_In_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_theta(°)_',num2str(theta_sim(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_Rsquare_',num2str(gof.rsquare,'%.2f'),'.tif');
            fullfilename = fullfile(save_directory,filename_tif); %Object fullfilename
            saveas(gcf,fullfilename); %Saving Matlab figure
            %Export to .obj format
            filename_obj = strcat('3D_Surface_Sessile_drop_Needle_In_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_theta(°)_',num2str(theta_sim(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_Rsquare_',num2str(gof.rsquare,'%.2f'),'.obj');
            fullfilename = fullfile(save_directory,filename_obj); %Object fullfilename
            saveobjmesh(fullfilename,X,Y,Z); %Exporting 3D drop surface as .obj
        %--------------------------------------------------------------
        %------------ SIMULATION OF DIFFERENT EDGE WIDTHS -------------
        %--------------------------------------------------------------
        elseif option == 3 %SIMULATION OF BLUR IN LIQUID DROP AND SOLID SUBSTRATE PROFILES (EDGE WIDTH)
            for q = 1:length(W_pxs_array)
                %Creation of sessile drop and needle edge normal profile
                x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)
                z_drop = z_ly_sim_pxs - (hsubstrato + z_ly_sim_pxs(end) - res_v);
                z_drop_lap_full = [wrev(z_drop);z_drop]; %Creation of drop profile (z coordinate)
                xCPleft_pxs(k) = x_drop_lap_full(1);
                xCPright_pxs(k) = x_drop_lap_full(end);
                yCPleft_pxs(k) = z_drop_lap_full(1);
                yCPright_pxs(k) = z_drop_lap_full(end);
                if needle_in
                    %Construction of the needle profile
                    needle_diameter_pxs = (needleDiam/10)*scale_cm;
                    pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                	pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                    indx_l = find((abs(x_drop_lap_full-pos_x_needle_l)) == min(abs(x_drop_lap_full-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                    indx_r = find((abs(x_drop_lap_full-pos_x_needle_r)) == min(abs(x_drop_lap_full-pos_x_needle_r)),1,'last'); %Find the index of the instersection between the right side of the needle and the drop
                    z_needle_max_l = z_drop_lap_full(indx_l);
                    z_needle_max_r = z_drop_lap_full(indx_r);
                    %z_needle_max = z_drop(end)-z_ly_sim_pxs(end)/2;
                    z_needle_l = linspace(0,z_needle_max_l,ceil(z_needle_max_l)); %z coordinates of the needle profile
                    z_needle_r = linspace(0,z_needle_max_r,ceil(z_needle_max_r)); %z coordinates of the needle profile
                    %z_needle = linspace(0,z_needle_max,ceil(z_needle_max)); %z coordinates of the needle profile
                    x_needle_l = ones(1,length(z_needle_l))*pos_x_needle_l; %x coordinates of the left needle profile
                    x_needle_r = ones(1,length(z_needle_r))*pos_x_needle_r; %x coordinates of the right needle profile
                    [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge(x_needle_l,z_needle_l,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the left needle profile
                    [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge(x_needle_r,z_needle_r,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the right needle profile
                    g_needle_r = wrev(g_needle_r); %Fix grayscale order in the right profile of the needle
                    xn_needle = [xn_needle_l xn_needle_r]; %x coordinates of edge and normal points of the needle profile
                    zn_needle = [zn_needle_l zn_needle_r]; %z coordinates of edge and normal points of the needle profile
                    g_needle = [g_needle_l g_needle_r]; %grayscale value of edge and normal points of the needle profile

                    %Trim drop profile based on the needle profile
                    x_drop_l = x_drop_lap_full(1:indx_l);
                    z_drop_l = z_drop_lap_full(1:indx_l);
                    x_drop_r = x_drop_lap_full(indx_r:end);
                    z_drop_r = z_drop_lap_full(indx_r:end);
                    [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the left portion of the drop
                    g_drop_l = wrev(g_drop_l);
                    [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the right portion of the drop
                    g_drop_r = wrev(g_drop_r);
                    xn_drop = [xn_drop_l xn_drop_r];
                    zn_drop = [zn_drop_l zn_drop_r];
                    g_drop = [g_drop_l g_drop_r];
                else
                    [xn_drop,zn_drop,g_drop] = grayscale_edge(x_drop_lap_full,z_drop_lap_full,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the drop profile
                    g_drop = wrev(g_drop);
                end

                %Creation of solid subratre edge normal profile
                x_sol_l = linspace(0,x_drop_lap_full(1),ceil(x_drop_lap_full(1)));
                z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                x_sol_r = linspace(x_drop_lap_full(end),res_h,ceil(res_h-x_drop_lap_full(end)));
                z_sol_r = ones(1,length(x_sol_r))*(res_v-hsubstrato);
                [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the left portion of the solid substrate surface
                g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
                [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the right portion of the solid substrate surface
                g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
                xn_sol = [xn_sol_l xn_sol_r];
                zn_sol = [zn_sol_l zn_sol_r];
                g_sol = [g_sol_l g_sol_r];

                %Trim drop normal edge according to solid substrate
                ind_z = find(zn_drop <= (res_v-hsubstrato));
                xn_drop = xn_drop(ind_z);
                zn_drop = zn_drop(ind_z);
                g_drop = g_drop(ind_z);

                %Combine solid, drop and needle edge and normal points
                if needle_in
                    xn = [xn_sol xn_needle xn_drop];
                    zn = [zn_sol zn_needle zn_drop];
                    g = [g_sol g_needle g_drop];
                else
                    xn = [xn_sol xn_drop];
                    zn = [zn_sol zn_drop];
                    g = [g_sol g_drop];
                end
                %Generation of synthetical illuminated (grayscale) pendant drop images
                M = zeros(res_v,res_h); %Array with image resolution size
                I = uint8(M); %Image creation
                %- Plotting solid substrate, drop profile and needle edge
                for i=1:length(g)
                    if (xn(i)>=0) && (xn(i)<res_h) && (zn(i)>=0) && (zn(i)<res_v)
                        if double(I(floor(zn(i))+1,floor(xn(i))+1)) == 0
                            I(floor(zn(i))+1,floor(xn(i))+1) = round(g(i));
                        else
                            I(floor(zn(i))+1,floor(xn(i))+1) = round((double(I(floor(zn(i))+1,floor(xn(i))+1)) + g(i))/2);
                        end
                    end
                end
                %- Filling background
                for i = 1:(res_v-hsubstrato) %Search from left to right
                    for j = 1:res_h
                        if I(i,j) == 0
                            I(i,j) = g2;
                        else
                            break;
                        end
                    end
                end
                for i = 1:(res_v-hsubstrato) %Search from right to left
                    for j = res_h:-1:1
                        if I(i,j) == 0
                            I(i,j) = g2;
                        else
                            break;
                        end
                    end
                end
                %- Filling solid substrate,drop and needle
                for i = 1:res_v %Search in all the image
                    for j = 1:res_h
                        if I(i,j) == 0
                            I(i,j) = g1;
                        end
                    end
                end
                %Save images
                imagename = strcat('Edge_Width_Grayscale_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_theta(°)_',num2str(theta_sim(k),'%.1f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs_array(q),'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
                fullimagename = strcat(imagename,imFmt); %Image name + image format
                fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                imwrite(I,fullfilename) %Save Image
            end
        %--------------------------------------------------------------
        %------- SIMULATION OF CAMERA VERTICAL MISALIGNMENT -----------
        %--------------------------------------------------------------
        elseif option == 4 % SIMULATION OF CAMERA VERTICAL MISALIGNMENT
            for q = 1:length(phi)
                x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)
                z_drop = z_ly_sim_pxs - (hsubstrato + z_ly_sim_pxs(end) - res_v);
                z_drop_lap_full = [wrev(z_drop);z_drop]; %Creation of drop profile (z coordinate)
                xCPleft_pxs(k) = x_drop_lap_full(1);
                xCPright_pxs(k) = x_drop_lap_full(end);
                yCPleft_pxs(k) = z_drop_lap_full(1);
                yCPright_pxs(k) = z_drop_lap_full(end);
                if needle_in
                    %Construction of the needle profile
                    needle_diameter_pxs = (needleDiam/10)*scale_cm;
                    pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                    pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                    indx_l = find((abs(x_drop_lap_full-pos_x_needle_l)) == min(abs(x_drop_lap_full-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                    indx_r = find((abs(x_drop_lap_full-pos_x_needle_r)) == min(abs(x_drop_lap_full-pos_x_needle_r)),1,'last'); %Find the index of the instersection between the right side of the needle and the drop
                    z_needle_max_l = z_drop_lap_full(indx_l);
                    z_needle_max_r = z_drop_lap_full(indx_r);
                    z_needle_l = linspace(-100,z_needle_max_l,ceil(z_needle_max_l)+100); %z coordinates of the needle profile
                    z_needle_r = linspace(-100,z_needle_max_r,ceil(z_needle_max_r)+100); %z coordinates of the needle profile
                    x_needle_l = ones(1,length(z_needle_l))*pos_x_needle_l; %x coordinates of the left needle profile
                    x_needle_r = ones(1,length(z_needle_r))*pos_x_needle_r; %x coordinates of the right needle profile
                    [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge(x_needle_l,z_needle_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left needle profile
                    [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge(x_needle_r,z_needle_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right needle profile
                    g_needle_r = wrev(g_needle_r); %Fix grayscale order in the right profile of the needle
                    xn_needle = [xn_needle_l xn_needle_r]; %x coordinates of edge and normal points of the needle profile
                    zn_needle = [zn_needle_l zn_needle_r]; %z coordinates of edge and normal points of the needle profile
                    g_needle = [g_needle_l g_needle_r]; %grayscale value of edge and normal points of the needle profile

                    %Trim drop profile based on the needle profile
                    x_drop_l = x_drop_lap_full(1:indx_l);
                    z_drop_l = z_drop_lap_full(1:indx_l);
                    x_drop_r = x_drop_lap_full(indx_r:end);
                    z_drop_r = z_drop_lap_full(indx_r:end);
                    [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the drop
                    g_drop_l = wrev(g_drop_l);
                    [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the drop
                    g_drop_r = wrev(g_drop_r);
                    xn_drop = [xn_drop_l xn_drop_r];
                    zn_drop = [zn_drop_l zn_drop_r];
                    g_drop = [g_drop_l g_drop_r];
                else
                    [xn_drop,zn_drop,g_drop] = grayscale_edge(x_drop_lap_full,z_drop_lap_full,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the drop profile
                    g_drop = wrev(g_drop);
                end

                %Creation of solid substrate edge normal profile
                x_sol_l = linspace(-100,x_drop_lap_full(1),ceil(x_drop_lap_full(1))+100);
                z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                x_sol_r = linspace(x_drop_lap_full(end),res_h+100,ceil(res_h-x_drop_lap_full(end)+100));
                z_sol_r = ones(1,length(x_sol_r))*(res_v-hsubstrato);
                [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the solid substrate surface
                g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
                [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the solid substrate surface
                g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
                xn_sol = [xn_sol_l xn_sol_r];
                zn_sol = [zn_sol_l zn_sol_r];
                g_sol = [g_sol_l g_sol_r];

                %Trim drop normal edge according to solid substrate
                ind_z = find(zn_drop <= (res_v-hsubstrato));
                xn_drop = xn_drop(ind_z);
                zn_drop = zn_drop(ind_z);
                g_drop = g_drop(ind_z);

                %Combine solid, drop and needle edge and normal points/ Rotation of data
                if needle_in
                    xn = [xn_sol xn_needle xn_drop];
                    zn = [zn_sol zn_needle zn_drop];
                    g = [g_sol g_needle g_drop];
                    [x_drop_l_rot,z_drop_l_rot] = rotate_data(x_drop_l,z_drop_l,xc_rot,zc_rot,phi(q));
                    [x_drop_r_rot,z_drop_r_rot] = rotate_data(x_drop_r,z_drop_r,xc_rot,zc_rot,phi(q));
                    [x_needle_l_rot,z_needle_l_rot] = rotate_data(x_needle_l,z_needle_l,xc_rot,zc_rot,phi(q));
                    [x_needle_r_rot,z_needle_r_rot] = rotate_data(x_needle_r,z_needle_r,xc_rot,zc_rot,phi(q));
                else
                    xn = [xn_sol xn_drop];
                    zn = [zn_sol zn_drop];
                    g = [g_sol g_drop];
                    [x_drop_lap_full_rot,z_drop_lap_full_rot] = rotate_data(x_drop_lap_full,z_drop_lap_full,xc_rot,zc_rot,phi(q));
                end
                [xn_rot,zn_rot] = rotate_data(xn,zn,xc_rot,zc_rot,phi(q)); %all data


                %Generation of synthetical rotated grayscale sessile drop images
                M = zeros(res_v,res_h); %Array with image resolution size
                I = uint8(M); %Image creation
                %- Plotting solid substrate, drop profile and needle edge
                for n=1:length(g)
                    if (xn_rot(n)>=0) && (xn_rot(n)<res_h) && (zn_rot(n)>=0) && (zn_rot(n)<res_v)
                        if double(I(floor(zn_rot(n))+1,floor(xn_rot(n))+1)) == 0
                            I(floor(zn_rot(n))+1,floor(xn_rot(n))+1) = round(g(n));
                        else
                            I(floor(zn_rot(n))+1,floor(xn_rot(n))+1) = round((double(I(floor(zn_rot(n))+1,floor(xn_rot(n))+1)) + g(n))/2);
                        end
                    end
                end
                %Sessile drop with needle
                if needle_in
                    if phi(q)<0 %Counterclockwise rotation
                        %Filling image inside needle
                        for n = 1:length(z_needle_l_rot)
                            if round(z_needle_l_rot(n)) > 0
                                for m = round(x_needle_l_rot(n)):round(res_h)
                                    if I(round(z_needle_l_rot(n)),m) == 0
                                        I(round(z_needle_l_rot(n)),m) = g1;
                                    elseif I(round(z_needle_l_rot(n)),m) == g2
                                        break;
                                    end
                                end
                            end
                        end
                        %Filling image inside drop
                        for n = 1:length(z_drop_l_rot)
                            for m = round(x_drop_l_rot(n)):round(res_h)
                                if I(round(z_drop_l_rot(n)),m) == 0
                                    I(round(z_drop_l_rot(n)),m) = g1;
                                elseif I(round(z_drop_l_rot(n)),m) == g2
                                    break;
                                end
                            end
                        end
                    else %Clockwise rotation
                        %Filling image inside needle
                        for n = 1:length(z_needle_r_rot)
                            if round(z_needle_r_rot(n)) > 0
                                for m = round(x_needle_r_rot(n)):-1:1
                                    if I(round(z_needle_r_rot(n)),m) == 0
                                        I(round(z_needle_r_rot(n)),m) = g1;
                                    elseif I(round(z_needle_r_rot(n)),m) == g2
                                        break;
                                    end
                                end
                            end
                        end
                        %Filling image inside drop
                        for n = 1:length(z_drop_r_rot)
                            for m = round(x_drop_r_rot(n)):-1:1
                                if I(round(z_drop_r_rot(n)),m) == 0
                                    I(round(z_drop_r_rot(n)),m) = g1;
                                elseif I(round(z_drop_r_rot(n)),m) == g2
                                    break;
                                end
                            end
                        end
                    end
                    %Sessile drop without needle
                else
                    %Finding drop apex and dividing drop profile
                    indx_z = find(z_drop_lap_full_rot == min(z_drop_lap_full_rot),1,'first'); %Find drop apex
                    x_drop_lap_l_rot = x_drop_lap_full_rot(1:indx_z);
                    z_drop_lap_l_rot = z_drop_lap_full_rot(1:indx_z);
                    x_drop_lap_r_rot = x_drop_lap_full_rot(indx_z:end);
                    z_drop_lap_r_rot = z_drop_lap_full_rot(indx_z:end);
                    if phi(q)<0 %Counterclockwise rotation
                        %Filling image inside drop
                        for n = 1:length(z_drop_lap_l_rot)
                            for m = round(x_drop_lap_l_rot(n)):round(res_h)
                                if I(round(z_drop_lap_l_rot(n)),m) == 0
                                    I(round(z_drop_lap_l_rot(n)),m) = g1;
                                elseif I(round(z_drop_lap_l_rot(n)),m) == g2
                                    break;
                                end
                            end
                        end
                    else %Clockwise rotation
                        %Filling image inside drop
                        for n = 1:length(z_drop_lap_r_rot)
                            for m = round(x_drop_lap_r_rot(n)):-1:1
                                if I(round(z_drop_lap_r_rot(n)),m) == 0
                                    I(round(z_drop_lap_r_rot(n)),m) = g1;
                                elseif I(round(z_drop_lap_r_rot(n)),m) == g2
                                    break;
                                end
                            end
                        end
                    end
                end
                %Filling solid substrate
                for m=1:res_h %Search from bottom to top
                    for n =res_v:-1:1
                        if I(n,m) == 0
                            I(n,m) = g1;
                        elseif I(n,m) == g2
                            break;
                        end
                    end
                end
                %Filling background
                for n = 1:res_v %Search in all the image
                    for m = 1:res_h
                        if I(n,m) == 0
                            I(n,m) = g2;
                        end
                    end
                end
                %Save images
                imagename = strcat('Vertical_Misalignment_Grayscale_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_theta(°)_',num2str(theta_sim(k),'%.1f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_phi(°)_',num2str(phi(q),'%.3f')); %Image filename
                fullimagename = strcat(imagename,imFmt); %Image name + image format
                fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                imwrite(I,fullfilename) %Save Image
            end
        %--------------------------------------------------------------
        %----- SIMULATION OF RANDOM PERTURBATIONS IN DROP PROFILE -----
        %--------------------------------------------------------------
        elseif option == 5 %Simulation of random perturbations in drop profile
            for p = 1:length(pert)
                for q = 1:n_images
                    %Creation of sessile drop and needle edge normal profile
                    x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)
                    z_drop = z_ly_sim_pxs - (hsubstrato + z_ly_sim_pxs(end) - res_v);
                    z_drop_lap_full = [wrev(z_drop);z_drop]; %Creation of drop profile (z coordinate)
                    xCPleft_pxs(k) = x_drop_lap_full(1);
                    xCPright_pxs(k) = x_drop_lap_full(end);
                    yCPleft_pxs(k) = z_drop_lap_full(1);
                    yCPright_pxs(k) = z_drop_lap_full(end);
                    if needle_in
                        %Construction of the needle profile
                        needle_diameter_pxs = (needleDiam/10)*scale_cm;
                        pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                        pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                        indx_l = find((abs(x_drop_lap_full-pos_x_needle_l)) == min(abs(x_drop_lap_full-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                        indx_r = find((abs(x_drop_lap_full-pos_x_needle_r)) == min(abs(x_drop_lap_full-pos_x_needle_r)),1,'last'); %Find the index of the instersection between the right side of the needle and the drop
                        z_needle_max_l = z_drop_lap_full(indx_l);
                        z_needle_max_r = z_drop_lap_full(indx_r);
                        %z_needle_max = z_drop(end)-z_ly_sim_pxs(end)/2;
                        z_needle_l = linspace(0,z_needle_max_l,ceil(z_needle_max_l)); %z coordinates of the needle profile
                        z_needle_r = linspace(0,z_needle_max_r,ceil(z_needle_max_r)); %z coordinates of the needle profile
                        %z_needle = linspace(0,z_needle_max,ceil(z_needle_max)); %z coordinates of the needle profile
                        x_needle_l = ones(1,length(z_needle_l))*pos_x_needle_l; %x coordinates of the left needle profile
                        x_needle_r = ones(1,length(z_needle_r))*pos_x_needle_r; %x coordinates of the right needle profile
                        [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge(x_needle_l,z_needle_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left needle profile
                        [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge(x_needle_r,z_needle_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right needle profile
                        g_needle_r = wrev(g_needle_r); %Fix grayscale order in the right profile of the needle
                        xn_needle = [xn_needle_l xn_needle_r]; %x coordinates of edge and normal points of the needle profile
                        zn_needle = [zn_needle_l zn_needle_r]; %z coordinates of edge and normal points of the needle profile
                        g_needle = [g_needle_l g_needle_r]; %grayscale value of edge and normal points of the needle profile

                        %Trim drop profile based on the needle profile
                        x_drop_l = x_drop_lap_full(1:indx_l);
                        z_drop_l = z_drop_lap_full(1:indx_l);
                        x_drop_r = x_drop_lap_full(indx_r:end);
                        z_drop_r = z_drop_lap_full(indx_r:end);
                        [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the drop
                        g_drop_l = wrev(g_drop_l);
                        [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the drop
                        g_drop_r = wrev(g_drop_r);
                        xn_drop = [xn_drop_l xn_drop_r];
                        zn_drop = [zn_drop_l zn_drop_r];
                        g_drop = [g_drop_l g_drop_r];
                    else
                        [xn_drop,zn_drop,g_drop] = grayscale_edge(x_drop_lap_full,z_drop_lap_full,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the drop profile
                        g_drop = wrev(g_drop);
                    end

                    %Random perturbation of the drop profile
                    [xn_drop,zn_drop,xn_drop_edge,zn_drop_edge] = normal_perturbation_edge(xn_drop,zn_drop,pert(p),n_normal);

                    %Creation of solid subratre edge normal profile
                    x_sol_l = linspace(0,x_drop_lap_full(1),ceil(x_drop_lap_full(1)));
                    z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                    x_sol_r = linspace(x_drop_lap_full(end),res_h,ceil(res_h-x_drop_lap_full(end)));
                    z_sol_r = ones(1,length(x_sol_r))*(res_v-hsubstrato);
                    [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the solid substrate surface
                    g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
                    [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the solid substrate surface
                    g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
                    xn_sol = [xn_sol_l xn_sol_r];
                    zn_sol = [zn_sol_l zn_sol_r];
                    g_sol = [g_sol_l g_sol_r];

                    %Trim drop normal edge according to solid substrate
                    ind_z = find(zn_drop <= (res_v-hsubstrato));
                    xn_drop = xn_drop(ind_z);
                    zn_drop = zn_drop(ind_z);
                    g_drop = g_drop(ind_z);

                    %Combine solid, drop and needle edge and normal points
                    if needle_in
                        xn = [xn_sol xn_needle xn_drop];
                        zn = [zn_sol zn_needle zn_drop];
                        g = [g_sol g_needle g_drop];
                    else
                        xn = [xn_sol xn_drop];
                        zn = [zn_sol zn_drop];
                        g = [g_sol g_drop];
                    end
                    %Generation of synthetical illuminated (grayscale) pendant drop images
                    M = zeros(res_v,res_h); %Array with image resolution size
                    I = uint8(M); %Image creation
                    %- Plotting solid substrate, drop profile and needle edge
                    for t = 1:length(g)
                        if (xn(t)>=0) && (xn(t)<res_h) && (zn(t)>=0) && (zn(t)<res_v)
                            if double(I(floor(zn(t))+1,floor(xn(t))+1)) == 0
                                I(floor(zn(t))+1,floor(xn(t))+1) = round(g(t));
                            else
                                I(floor(zn(t))+1,floor(xn(t))+1) = round((double(I(floor(zn(t))+1,floor(xn(t))+1)) + g(t))/2);
                            end
                        end
                    end
                    if needle_in == 1
                        for t = 1:length(z_needle_r)
                            if round(z_needle_r(t)) ~= 0
                                for m = round(x_needle_r(t)):-1:1
                                    if I(round(z_needle_r(t)),m) == 0
                                        I(round(z_needle_r(t)),m) = g1;
                                    elseif I(round(z_needle_r(t)),m) == g2
                                        break;
                                    end
                                end
                            end
                        end
                        %Find drop apex
                        indx_apex_sup = find(abs(zn_drop_edge) == min(abs(zn_drop_edge)),1,'last');
                        indx_apex_inf = find(abs(zn_drop_edge) == max(abs(zn_drop_edge)),1,'last');
                        %Divide drop profile and determine the right side
                        xn_drop_edge_r = xn_drop_edge(indx_apex_sup:indx_apex_inf);
                        zn_drop_edge_r = zn_drop_edge(indx_apex_sup:indx_apex_inf);
                        %Filling image inside drop
                        for t = 1:length(zn_drop_edge_r)
                            for m = round(xn_drop_edge_r(t)):-1:1
                                if I(round(zn_drop_edge_r(t)),m) == 0
                                    I(round(zn_drop_edge_r(t)),m) = g1;
                                elseif I(round(zn_drop_edge_r(t)),m) == g2
                                    break;
                                end
                            end
                        end
                    else
                        %Find drop apex
                        indx_apex = find(abs(zn_drop_edge) == min(abs(zn_drop_edge)),1,'first');
                        %Divide drop profile and determine the right side
                        xn_drop_edge_r = xn_drop_edge(indx_apex:end);
                        zn_drop_edge_r = zn_drop_edge(indx_apex:end);
                        %Filling image inside drop
                        for t = 1:length(zn_drop_edge_r)
                            for m = round(xn_drop_edge_r(t)):-1:1
                                if I(round(zn_drop_edge_r(t)),m) == 0
                                    I(round(zn_drop_edge_r(t)),m) = g1;
                                elseif I(round(zn_drop_edge_r(t)),m) == g2
                                    break;
                                end
                            end
                        end
                    end
                    %Filling solid substrate
                    for m = 1:res_h %Search from bottom to top - left to right
                        for t = res_v:-1:1
                            if I(t,m) == g2
                                break;
                            elseif I(t,m) == 0
                                I(t,m) = g1;
                            end
                        end
                    end
                    %- Filling background
                    for t = 1:res_v %Search in all the image
                        for m = 1:res_h
                            if I(t,m) == 0
                                I(t,m) = g2;
                            end
                        end
                    end
                    %Save images
                    imagename = strcat('Random_Edge_Perturbation_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_theta(°)_',num2str(theta_sim(k),'%.1f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_pert(pxs)_',num2str(pert(p),'%.2f'),'_image_',num2str(q),'.tif'); %Image filename
                    fullimagename = strcat(imagename,imFmt); %Image name + image format
                    fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                    imwrite(I,fullfilename) %Save Image
                end
            end
        %--------------------------------------------------------------
        %--------------- SIMULATION OF LACK OF CONTRAST ---------------
        %--------------------------------------------------------------
        elseif option == 6 %Simulation of lack of contrast
            for n = 1:length(diff_g)
                g2 = g1 + diff_g(n); %Setting g2 value
                %Creation of sessile drop and needle edge normal profile
                x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)
                z_drop = z_ly_sim_pxs - (hsubstrato + z_ly_sim_pxs(end) - res_v);
                z_drop_lap_full = [wrev(z_drop);z_drop]; %Creation of drop profile (z coordinate)
                xCPleft_pxs(k) = x_drop_lap_full(1);
                xCPright_pxs(k) = x_drop_lap_full(end);
                yCPleft_pxs(k) = z_drop_lap_full(1);
                yCPright_pxs(k) = z_drop_lap_full(end);
                if needle_in
                    %Construction of the needle profile
                    needle_diameter_pxs = (needleDiam/10)*scale_cm;
                    pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                    pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                    indx_l = find((abs(x_drop_lap_full-pos_x_needle_l)) == min(abs(x_drop_lap_full-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                    indx_r = find((abs(x_drop_lap_full-pos_x_needle_r)) == min(abs(x_drop_lap_full-pos_x_needle_r)),1,'last'); %Find the index of the instersection between the right side of the needle and the drop
                    z_needle_max_l = z_drop_lap_full(indx_l);
                    z_needle_max_r = z_drop_lap_full(indx_r);
                    z_needle_l = linspace(0,z_needle_max_l,ceil(z_needle_max_l)); %z coordinates of the needle profile
                    z_needle_r = linspace(0,z_needle_max_r,ceil(z_needle_max_r)); %z coordinates of the needle profile
                    x_needle_l = ones(1,length(z_needle_l))*pos_x_needle_l; %x coordinates of the left needle profile
                    x_needle_r = ones(1,length(z_needle_r))*pos_x_needle_r; %x coordinates of the right needle profile
                    [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge(x_needle_l,z_needle_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left needle profile
                    [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge(x_needle_r,z_needle_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right needle profile
                    g_needle_r = wrev(g_needle_r); %Fix grayscale order in the right profile of the needle
                    xn_needle = [xn_needle_l xn_needle_r]; %x coordinates of edge and normal points of the needle profile
                    zn_needle = [zn_needle_l zn_needle_r]; %z coordinates of edge and normal points of the needle profile
                    g_needle = [g_needle_l g_needle_r]; %grayscale value of edge and normal points of the needle profile

                    %Trim drop profile based on the needle profile
                    x_drop_l = x_drop_lap_full(1:indx_l);
                    z_drop_l = z_drop_lap_full(1:indx_l);
                    x_drop_r = x_drop_lap_full(indx_r:end);
                    z_drop_r = z_drop_lap_full(indx_r:end);
                    [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the drop
                    g_drop_l = wrev(g_drop_l);
                    [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the drop
                    g_drop_r = wrev(g_drop_r);
                    xn_drop = [xn_drop_l xn_drop_r];
                    zn_drop = [zn_drop_l zn_drop_r];
                    g_drop = [g_drop_l g_drop_r];
                else
                    [xn_drop,zn_drop,g_drop] = grayscale_edge(x_drop_lap_full,z_drop_lap_full,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the drop profile
                    g_drop = wrev(g_drop);
                end

                %Creation of solid subratre edge normal profile
                x_sol_l = linspace(0,x_drop_lap_full(1),ceil(x_drop_lap_full(1)));
                z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                x_sol_r = linspace(x_drop_lap_full(end),res_h,ceil(res_h-x_drop_lap_full(end)));
                z_sol_r = ones(1,length(x_sol_r))*(res_v-hsubstrato);
                [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the solid substrate surface
                g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
                [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the solid substrate surface
                g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
                xn_sol = [xn_sol_l xn_sol_r];
                zn_sol = [zn_sol_l zn_sol_r];
                g_sol = [g_sol_l g_sol_r];

                %Trim drop normal edge according to solid substrate
                ind_z = find(zn_drop <= (res_v-hsubstrato));
                xn_drop = xn_drop(ind_z);
                zn_drop = zn_drop(ind_z);
                g_drop = g_drop(ind_z);

                %Combine solid, drop and needle edge and normal points
                if needle_in
                    xn = [xn_sol xn_needle xn_drop];
                    zn = [zn_sol zn_needle zn_drop];
                    g = [g_sol g_needle g_drop];
                else
                    xn = [xn_sol xn_drop];
                    zn = [zn_sol zn_drop];
                    g = [g_sol g_drop];
                end
                %Generation of synthetical illuminated (grayscale) pendant drop images
                M = zeros(res_v,res_h); %Array with image resolution size
                I = uint8(M); %Image creation
                %- Plotting solid substrate, drop profile and needle edge
                for i=1:length(g)
                    if (xn(i)>=0) && (xn(i)<res_h) && (zn(i)>=0) && (zn(i)<res_v)
                        if double(I(floor(zn(i))+1,floor(xn(i))+1)) == 0
                            I(floor(zn(i))+1,floor(xn(i))+1) = round(g(i));
                        else
                            I(floor(zn(i))+1,floor(xn(i))+1) = round((double(I(floor(zn(i))+1,floor(xn(i))+1)) + g(i))/2);
                        end
                    end
                end
                %- Filling background
                for i = 1:(res_v-hsubstrato) %Search from left to right
                    for j = 1:res_h
                        if I(i,j) == 0
                            I(i,j) = g2;
                        else
                            break;
                        end
                    end
                end
                for i = 1:(res_v-hsubstrato) %Search from right to left
                    for j = res_h:-1:1
                        if I(i,j) == 0
                            I(i,j) = g2;
                        else
                            break;
                        end
                    end
                end
                %- Filling solid substrate,drop and needle
                for i = 1:res_v %Search in all the image
                    for j = 1:res_h
                        if I(i,j) == 0
                            I(i,j) = g1;
                        end
                    end
                end
                %Save images
                imagename = strcat('Grayscale_SessilGerar e_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_theta(°)_',num2str(theta_sim(k),'%.1f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
                fullimagename = strcat(imagename,imFmt); %Image name + image format
                fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                imwrite(I,fullfilename) %Save Image
            end
        %--------------------------------------------------------------
        %--------- SIMULATION OF NON-UNIFORM ILLUMINATION -------------
        %--------------------------------------------------------------
        elseif option == 7 %Simulation of non-uniform illumination
            iref = res_v/2; %Coordenada no eixo i do centro da iluminação não-uniforme
            jref = res_h/2; %Coordenada no eixo j do centro da iluminação não-uniforme
            for n = 1:length(deltag)
                %Creation of sessile drop and needle edge normal profile
                x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)
                z_drop = z_ly_sim_pxs - (hsubstrato + z_ly_sim_pxs(end) - res_v);
                z_drop_lap_full = [wrev(z_drop);z_drop]; %Creation of drop profile (z coordinate)
                xCPleft_pxs(k) = x_drop_lap_full(1);
                xCPright_pxs(k) = x_drop_lap_full(end);
                yCPleft_pxs(k) = z_drop_lap_full(1);
                yCPright_pxs(k) = z_drop_lap_full(end);
                if needle_in
                    %Construction of the needle profile
                    needle_diameter_pxs = (needleDiam/10)*scale_cm;
                    pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                    pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                    indx_l = find((abs(x_drop_lap_full-pos_x_needle_l)) == min(abs(x_drop_lap_full-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                    indx_r = find((abs(x_drop_lap_full-pos_x_needle_r)) == min(abs(x_drop_lap_full-pos_x_needle_r)),1,'last'); %Find the index of the instersection between the right side of the needle and the drop
                    z_needle_max_l = z_drop_lap_full(indx_l);
                    z_needle_max_r = z_drop_lap_full(indx_r);
                    z_needle_l = linspace(0,z_needle_max_l,ceil(z_needle_max_l)); %z coordinates of the needle profile
                    z_needle_r = linspace(0,z_needle_max_r,ceil(z_needle_max_r)); %z coordinates of the needle profile
                    x_needle_l = ones(1,length(z_needle_l))*pos_x_needle_l; %x coordinates of the left needle profile
                    x_needle_r = ones(1,length(z_needle_r))*pos_x_needle_r; %x coordinates of the right needle profile
                    [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge_nonuniform_illumination(x_needle_l,z_needle_l,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the left needle profile
                    [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge_nonuniform_illumination(x_needle_r,z_needle_r,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the right needle profile
                    %g_needle_r = wrev(g_needle_r); %Fix grayscale order in the right profile of the needle
                    for p = 1:(length(g_needle_r)/n_normal) %Correcting g_sol_r order
                        prov = g_needle_r(n_normal*(p-1)+1:n_normal*p);
                        for m = 1:floor(n_normal/2)
                            g_needle_r((n_normal*(p-1))+m) = prov((n_normal)-(m-1));
                            g_needle_r((n_normal*p)-(m-1)) = prov(m);
                        end
                    end
                    xn_needle = [xn_needle_l xn_needle_r]; %x coordinates of edge and normal points of the needle profile
                    zn_needle = [zn_needle_l zn_needle_r]; %z coordinates of edge and normal points of the needle profile
                    g_needle = [g_needle_l g_needle_r]; %grayscale value of edge and normal points of the needle profile

                    %Trim drop profile based on the needle profile
                    x_drop_l = x_drop_lap_full(1:indx_l);
                    z_drop_l = z_drop_lap_full(1:indx_l);
                    x_drop_r = x_drop_lap_full(indx_r:end);
                    z_drop_r = z_drop_lap_full(indx_r:end);
                    [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge_nonuniform_illumination(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the left portion of the drop
                    %g_drop_l = wrev(g_drop_l);
                    for p = 1:(length(g_drop_l)/n_normal) %Correcting g_sol_r order
                        prov = g_drop_l(n_normal*(p-1)+1:n_normal*p);
                        for m = 1:floor(n_normal/2)
                            g_drop_l((n_normal*(p-1))+m) = prov((n_normal)-(m-1));
                            g_drop_l((n_normal*p)-(m-1)) = prov(m);
                        end
                    end
                    [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge_nonuniform_illumination(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the right portion of the drop
                    %g_drop_r = wrev(g_drop_r);
                    for p = 1:(length(g_drop_r)/n_normal) %Correcting g_sol_r order
                        prov = g_drop_r(n_normal*(p-1)+1:n_normal*p);
                        for m = 1:floor(n_normal/2)
                            g_drop_r((n_normal*(p-1))+m) = prov((n_normal)-(m-1));
                            g_drop_r((n_normal*p)-(m-1)) = prov(m);
                        end
                    end
                    xn_drop = [xn_drop_l xn_drop_r];
                    zn_drop = [zn_drop_l zn_drop_r];
                    g_drop = [g_drop_l g_drop_r];
                else
                    [xn_drop,zn_drop,g_drop] = grayscale_edge_nonuniform_illumination(x_drop_lap_full,z_drop_lap_full,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the drop profile
                    g_drop = wrev(g_drop);
                end

                %Creation of solid subratre edge normal profile
                x_sol_l = linspace(0,x_drop_lap_full(1),ceil(x_drop_lap_full(1)));
                z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                x_sol_r = linspace(x_drop_lap_full(end),res_h,ceil(res_h-x_drop_lap_full(end)));
                z_sol_r = ones(1,length(x_sol_r))*(res_v-hsubstrato);
                [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge_nonuniform_illumination(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the left portion of the solid substrate surface
                %g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
                for p = 1:(length(g_sol_l)/n_normal) %Correcting g_sol_r order
                    prov = g_sol_l(n_normal*(p-1)+1:n_normal*p);
                    for m = 1:floor(n_normal/2)
                        g_sol_l((n_normal*(p-1))+m) = prov((n_normal)-(m-1));
                        g_sol_l((n_normal*p)-(m-1)) = prov(m);
                    end
                end
                [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge_nonuniform_illumination(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the right portion of the solid substrate surface
                %g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
                for p = 1:(length(g_sol_r)/n_normal) %Correcting g_sol_r order
                    prov = g_sol_r(n_normal*(p-1)+1:n_normal*p);
                    for m = 1:floor(n_normal/2)
                        g_sol_r((n_normal*(p-1))+m) = prov((n_normal)-(m-1));
                        g_sol_r((n_normal*p)-(m-1)) = prov(m);
                    end
                end
                xn_sol = [xn_sol_l xn_sol_r];
                zn_sol = [zn_sol_l zn_sol_r];
                g_sol = [g_sol_l g_sol_r];
                %Trim drop normal edge according to solid substrate
                ind_z = find(zn_drop <= (res_v-hsubstrato));
                xn_drop = xn_drop(ind_z);
                zn_drop = zn_drop(ind_z);
                g_drop = g_drop(ind_z);

                %Combine solid, drop and needle edge and normal points
                if needle_in
                    xn = [xn_sol xn_needle xn_drop];
                    zn = [zn_sol zn_needle zn_drop];
                    g = [g_sol g_needle g_drop];
                else
                    xn = [xn_sol xn_drop];
                    zn = [zn_sol zn_drop];
                    g = [g_sol g_drop];
                end
                %Generation of synthetical illuminated (grayscale) pendant drop images
                M = zeros(res_v,res_h); %Array with image resolution size
                I = uint8(M); %Image creation
                %- Plotting solid substrate, drop profile and needle edge
                for i=1:length(g)
                    if (xn(i)>=0) && (xn(i)<res_h) && (zn(i)>=0) && (zn(i)<res_v)
                        if double(I(floor(zn(i))+1,floor(xn(i))+1)) == 0
                            I(floor(zn(i))+1,floor(xn(i))+1) = round(g(i));
                        else
                            I(floor(zn(i))+1,floor(xn(i))+1) = round((double(I(floor(zn(i))+1,floor(xn(i))+1)) + g(i))/2);
                        end
                    end
                end
                %Iorig = I;
                %- Filling background
                for i = 1:(res_v-hsubstrato) %Search from left to right
                    for j = 1:res_h
                        if I(i,j) == 0
                            r = sqrt((i-iref)^2+(j-jref)^2);
                            I(i,j) = round(g2 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                            %I(i,j) = g2;
                        else
                            break;
                        end
                    end
                end
                for i = 1:(res_v-hsubstrato) %Search from right to left
                    for j = res_h:-1:1
                        if I(i,j) == 0
                            r = sqrt((i-iref)^2+(j-jref)^2);
                            I(i,j) = round(g2 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                            %I(i,j) = g2;
                        else
                            break;
                        end
                    end
                end
                %- Filling solid substrate,drop and needle
                for i = 1:res_v %Search in all the image
                    for j = 1:res_h
                        if I(i,j) == 0
                            r = sqrt((i-iref)^2+(j-jref)^2);
                            I(i,j) = round(g1 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                            %I(i,j) = g1;
                        end
                    end
                end
                %Save images
                imagename = strcat('Non_uniform_Grayscale_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_theta(°)_',num2str(theta_sim(k),'%.1f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_deltag_',num2str(round(deltag(n)))); %Image filename
                fullimagename = strcat(imagename,imFmt); %Image name + image format
                fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                imwrite(I,fullfilename) %Save Image
            end
        end     
    end
    r0 = 1./b; %radius at the apex in cm
    %Write excel file with important pendant drop profile properties
    fprintf('Exporting excel file containing important sessile drop properties... \n');
    varnames = {'b[cm-1]','r0[cm]','theta[°]','xCP_left[pxs]','yCP_left[pxs]','xCP_right[pxs]','yCP_right[pxs]','r_w[cm]','r_eq[cm]','h[cm]','A [cm2]','V[uL]'}; %Variables name in the excel file
    T = table(b,r0,theta_sim,xCPleft_pxs,yCPleft_pxs,xCPright_pxs,yCPright_pxs,r_w_cm,r_eq_cm,h_cm,A_cm2,V_uL,'VariableNames',varnames); %Build a table with important drop properties
    sheet_name = strcat('c = ',num2str(c,'%.2f'),' cm-2'); %Excel sheet name
    savename_xls = strcat('Dimensions_sessiledrop(c(cm-2)_',num2str(c),'_theta(°)_',num2str(mean(theta),'%.2f'),').xlsx'); %Excel file name
    fullsavename_xls = fullfile(save_directory,savename_xls);
    writetable(T,fullsavename_xls,'Sheet',sheet_name) %Write table to an excel file

    fprintf('Simulation completed! \n');
    pause(2)
    end
end
end