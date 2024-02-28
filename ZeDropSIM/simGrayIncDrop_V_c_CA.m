function simGrayIncDrop_V_c_CA()
%SIMGRAYINCDROP_V_C_CA Simulation of grayscale inclined drop images for 
% different drop volume (V) given liquid capillary constant (c), 
% inclination angle (alpha) and contact angle (CA).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of grayscale inclined drop images for different drop volumes 
%   (V) given liquid capillary constant (c), inclination angle (alpha) and 
%   contact angle (CA). Some source of  errors, such as edge width, camera 
%   vertical misalignment, random perturbation in the drop profile, lack of
%   contrast, non-uniform illumination, random noise, Gaussian noise, 
%   impulse salt-pepper noise and satellite droplets, can be simulated. 
%   Depending on the simulated option, the program requires the user to 
%   enter some specific parameters. As common input parameters, any option 
%   requires the user to enter: the liquid capillary constant [cm-2](c), 
%   the inclination angle [°](CA), the possibility of simulating images 
%   inclined or planed solid substrate (inclined), the image resolution 
%   [pixels](res_h and res_v) and the image scale [pixels/mm](scale), edge 
%   width (W_pxs), maximum and minimum image intensity levels (g1 and g2). 
%   In order to make the execution of the algorithm faster, it allows the 
%   user to load an excel file containing essential drop properties: 
%   curvature at the apex (b), drop volume (V) and contact angle (CA). The 
%   height of the solid substrate in the image is fixed. The routine also 
%   allows the user to choose the directory where the images are going to 
%   be exported. An excel file (.xlsx) is created containing important 
% properties of the inclined drops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ INCLINED DROP - SIMULATION OF GRAYSCALE IMAGES --------------------\n');
    % Options
    fprintf('OPTIONS: \n');
    fprintf('1. Generation of grayscale images; \n');
    %fprintf('2. Generation of 3D drop surface; \n');
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
    elseif option == 2

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
        alpha_deg = input("- Inclination angle in degree [15]: ");
        if isempty(alpha_deg)
            alpha_deg = 15;
        end
        %fprintf('- Load excel file containing curvature at the apex (b), drop volume (V) and contact angle (CA) values ...\n ');
        loadExc = input('- Load excel file containing curvature at the apex (b), drop volume (V) and contact angle (CA) values (Y/N) [N]: ', 's');
        if isempty(loadExc)
        loadExc = 'N';
        end
        if loadExc == 'Y'
            [openname_xls,path] = uigetfile('*.xlsx'); %Get file path to open excel file
            fullopenname_xls = fullfile(path,openname_xls);
            readT = readtable(fullopenname_xls,"VariableNamingRule","preserve"); %Read excel file
            b = readT{:,1};
            r0 = 1./b; %radius at the apex in cm
            CA = readT{:,3};
            V = readT{:,9};
        else
            CA = input("- Contact angle in degree [70]: "); %Contact angle in degrees
            if isempty(CA)
                CA = 70;
            end
            Vmin = input("- Minimum drop volume in uL [2]: "); %Minimum drop volume in uL
            if isempty(Vmin)
                Vmin = 2;
            end
            while 1
                Vmax = input("- Maximum drop volume in uL [10]: "); %Maximum drop volume in uL
                if isempty(Vmax)
                    Vmax = 10;
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
            bmin_ini = input("- Minimum searched curvature at the apex in cm-1: "); %Minimum curvature at the apex possible in cm-1
            %{
        if isempty(bmin)
            bmin = 7.2; %For water
        end
            %}
            while 1
                bmax_ini = input("- Maximum searched curvature at the apex in cm-1: "); %Maximum curvature at the apex possible in cm-1;
                %{
            if isempty(bmax)
                bmax = 9.7; %For water
            end
                %}
                if bmax_ini < bmin_ini
                    fprintf('Invalid maximum curavture at the apex. Please, enter a lower value. \n');
                else
                    break;
                end
            end
            interval_b = 50; %Number of b values created between bmin and bmax for drop volume search
            tolV = 1e-3; %Tolerance for volume determination in uL
            itermax = 5; %Maximum number of iterations for volume search
            nSimImages = length(V); %Number of simulated images
            b = zeros(nSimImages,1);
            r0 = zeros(nSimImages,1);
        end
        fprintf('2) Drop configuration propertes: \n');
        inclined = input("- Do you want to show the solid substrate inclined (Y/N)? [Y]: ", "s"); %Generate images with substrate inclined (inclined == 'Y') or planed (inclined == 'N')
        if isempty(inclined) || (inclined ~= 'Y' && inclined ~= 'N')
            inclined = 'Y';
        end
        if inclined == 'N'
            inclined = 0;
        else
            inclined = 1;
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
        dist_pxs = W_pxs/2; %%Cabezas et al 2007 %Distance between normal edge points in pixels
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
        N = 13; %number of elements along theta/ Suggestion: give n odd numbers, so the middle element is given at theta=90°

        %- SPECIFIC INPUT PARAMETERS
        %{
        if option == 2 %Generation of 3D drop surface
            %- 3D Surface/mesh properties
            fprintf('5) 3D Surface/mesh properties: \n');
            n = input("- Number of points (mesh) around droplet circunference [100]: ");
            if isempty(n)
                n = 100; %number of points around the cilinder circunference of drop profile
            end
        elseif option == 3 %Simulation of edge width
        %}
        if option == 3 %Simulation of edge width
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
        CA_sim = zeros(length(b),1); %Intrinsic apparent CA for a W distance from the apex of the drop in degrees
        CA_adv = zeros(length(b),1); %Advancing CA in degrees
        CA_rec = zeros(length(b),1); %Receding CA in degrees
        Xcp_left = zeros(length(b),1); %X coordinate of the left triple contact point in pxs
        Ycp_left = zeros(length(b),1); %Y coordinate of the left triple contact point in pxs
        Xcp_right = zeros(length(b),1); %X coordinate of the right triple contact point in pxs
        Ycp_right = zeros(length(b),1); %Y coordinate of the right triple contact point in pxs
        d_w_cm = zeros(length(b),1); %Distance between contact points in cm
        h_cm = zeros(length(b),1); %Drop height in cm
        A_cm2 = zeros(length(b),1); %Drop surface area in cm²
        V_uL = zeros(length(b),1); %Drop volume in mm³ or uL
        
    %------------------------------------------------------------------------
    %----------------- SIMULATION OF INCLINED DROP IMAGES -------------------
    %------------------------------------------------------------------------
    %- Generation of grayscale inclined drop images
    %fprintf('Generating binary sessile drop images... \n');
    for k = 1:length(b)
        fprintf('----------------------------------------------------------------------------------------- \n');     
        fprintf('Simulating inclined drop profile for c = %.2f cm-2, alpha = %.2f°, CA = %.2f° and V = %.2f uL ... \n',c,alpha_deg,mean(CA),V(k));
        if loadExc == 'Y'
            %Determination of the inclined drop profile by using numerical solution technique for solving Young-Laplace PDE
            [x_ly_sim_cm,z_ly_sim_cm,W,U,theta0,CA_sim(k),CA_adv(k),CA_rec(k),d_w_cm(k),h_cm(k),A_cm2(k),V_uL(k)] = profile_inclined_drop_ly(b(k),c,alpha_deg,CA(k),N);
        else
            % Variable initialization
            x_ly_sim_cm_prov = cell(1,interval_b);
            z_ly_sim_cm_prov = cell(1,interval_b);
            CA_sim_prov = zeros(1,interval_b); %Intrinsic apparent CA for a W distance from the apex of the drop in degrees
            CA_adv_prov = zeros(1,interval_b); %Advancing CA in degrees
            CA_rec_prov = zeros(1,interval_b); %Receding CA in degrees
            d_w_cm_prov = zeros(1,interval_b); %Distance between contact points in cm
            h_cm_prov = zeros(1,interval_b); %Drop height in cm
            A_cm2_prov = zeros(1,interval_b); %Drop surface area in cm²
            V_uL_prov = zeros(1,interval_b); %Drop volume in mm³ or uL

            % Search for b value
            diffV = 1; %reset volume difference
            bmin = bmin_ini; %reset bmin
            bmax = bmax_ini; %reset bmax
            iter = 1;
            while diffV > tolV && iter <= itermax
                fprintf('Searching for V = %.2f uL, between bmin = %.4f and bmax = %.4f (iteration %d)... \n',V(k),bmin,bmax,iter);
                b_prov = linspace(bmin,bmax,interval_b);
                r0_prov = 1./b_prov;
                for u = 1:length(b_prov)
                    %Determination of the inclined drop profile by using numerical solution technique for solving Young-Laplace PDE
                    [x_ly_sim_cm_prov{u},z_ly_sim_cm_prov{u},W,U,theta0,CA_sim_prov(u),CA_adv_prov(u),CA_rec_prov(u),d_w_cm_prov(u),h_cm_prov(u),A_cm2_prov(u),V_uL_prov(u)] = profile_inclined_drop_ly(b_prov(u),c,alpha_deg,CA,N);
                end
                diffV = min(abs(V_uL_prov-V(k)));
                fprintf('diffV = %f\n',diffV);
                indx_V = find(abs(V_uL_prov-V(k)) == diffV,1,'last');
                if (indx_V - 1) <= 0
                    bmin = b_prov(indx_V);
                    bmax = b_prov(indx_V+1);
                elseif (indx_V + 1) >= length(b_prov)
                    bmin = b_prov(indx_V-1);
                    bmax = b_prov(indx_V);
                else
                    bmin = b_prov(indx_V-1);
                    bmax = b_prov(indx_V+1);
                end
                iter = iter +1;
            end
            V_uL(k) = V_uL_prov(indx_V);
            x_ly_sim_cm = x_ly_sim_cm_prov{indx_V};
            z_ly_sim_cm = z_ly_sim_cm_prov{indx_V};
            b(k) = b_prov(indx_V);
            r0(k) = r0_prov(indx_V);
            CA_sim(k) = CA_sim_prov(indx_V);
            CA_adv(k) = CA_adv_prov(indx_V);
            CA_rec(k) = CA_rec_prov(indx_V);
            d_w_cm(k) = d_w_cm_prov(indx_V);
            h_cm(k) = h_cm_prov(indx_V);
            A_cm2(k) = A_cm2_prov(indx_V);
        end
        
        %Determination of the inclined drop profile by using numerical solution technique for solving Young-Laplace PDE
        %[x_ly_sim_cm,z_ly_sim_cm,W,U,theta0,CA_sim(k),CA_adv(k),CA_rec(k),d_w_cm(k),h_cm(k),A_cm2(k),V_uL(k)] = profile_inclined_drop_ly(b(k),c,alpha_deg,CA(k),N);
        
        %Conversion to pixels
        x_ly_sim_pxs = x_ly_sim_cm*scale_cm;   %Drop width in pixels
        z_ly_sim_pxs = z_ly_sim_cm*scale_cm;   %Drop height in pixels
        %Construction of synthetical image
        res_v = resolution(2);
        res_h = resolution(1);
        xc = round(res_h/2); %Coordenada x do centro da imagem em pixels
        zc = round(res_v/2); %Coordenada z do centro da imagem em pixels

        %--------------------------------------------------------------
        %-------------- GENERATION OF GRAYSCALE IMAGES ----------------
        %--------------------------------------------------------------
        if option == 1 || option == 8 || option == 9 || option == 10 || option == 11 % Generation of grayscale images
            %Creation of the inclined drop normal edge and normal profile
            x_plot = x_ly_sim_pxs + xc - (x_ly_sim_pxs(1)+(x_ly_sim_pxs(end)-x_ly_sim_pxs(1))/2);

            %Defining inclined drop profile position in the image
            if inclined
                z_plot = z_ly_sim_pxs + res_v/2 - (max(abs(z_ly_sim_pxs))/2);
                z_subst = z_plot(end);
            else
                z_plot = z_ly_sim_pxs + res_v - (max(abs(z_ly_sim_pxs))+hsubstrato);
                z_subst = res_v-hsubstrato;
            end
            [xn_drop,zn_drop,g_drop] = grayscale_edge(x_plot,z_plot,n_normal,g1,g2,W_pxs); %Calculation of the normal points to inclined drop profile
            g_drop = wrev(g_drop);

            %Creation of solid subratre edge normal profile
            %x_sol_l = linspace(-res_h/2,x_plot(1),ceil(x_plot(1)+res_h/2));
            x_sol_l = linspace(-res_h/2,x_plot(1),10*ceil(x_plot(1)+res_h/2));
            z_sol_l = ones(1,length(x_sol_l))*(z_subst);
            %x_sol_r = linspace(x_plot(end),3*res_h/2,ceil((3*res_h/2)-x_plot(end)));
            x_sol_r = linspace(x_plot(end),3*res_h/2,10*ceil((3*res_h/2)-x_plot(end)));
            z_sol_r = ones(1,length(x_sol_r))*(z_subst);
            [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the solid substrate surface
            g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
            [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the solid substrate surface
            g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
            xn_sol = [xn_sol_l xn_sol_r];
            zn_sol = [zn_sol_l zn_sol_r];
            g_sol = [g_sol_l g_sol_r];

            %Trim drop normal edge according to solid substrate
            %ind_z = find(zn_drop <= (res_v-hsubstrato));
            ind_z = find(zn_drop <= z_subst);
            xn_drop = xn_drop(ind_z);
            zn_drop = zn_drop(ind_z);
            g_drop = g_drop(ind_z);

            %Combine solid and drop edge and normal points
            xn = [xn_sol xn_drop];
            zn = [zn_sol zn_drop];
            g = [g_sol g_drop];

            %Rotation of edge and normal points
            if inclined
                [xn,zn] = rotate_data(xn,zn,xc,zc,alpha_deg);
            end

            %Determination of  triple contact points
            if inclined
                [x_plot_rot,z_plot_rot] = rotate_data(x_plot,z_plot,xc,zc,alpha_deg);
                Xcp_left(k) = x_plot_rot(1);
                Ycp_left(k) = z_plot_rot(1);
                Xcp_right(k) = x_plot_rot(end);
                Ycp_right(k) = z_plot_rot(end);
            else
                Xcp_left(k) = x_plot(1);
                Ycp_left(k) = z_plot(1);
                Xcp_right(k) = x_plot(end);
                Ycp_right(k) = z_plot(end);
            end

            %Generation of synthetical illuminated (grayscale) pendant drop images
            M = zeros(res_v,res_h); %Array with image resolution size
            I = uint8(M); %Image creation
            %- Plotting solid substrate and drop profile edge
            for i=1:length(g)
                if (xn(i)>=0) && (xn(i)<res_h) && (zn(i)>=0) && (zn(i)<res_v)
                    if double(I(floor(zn(i))+1,floor(xn(i))+1)) == 0
                        I(floor(zn(i))+1,floor(xn(i))+1) = round(g(i));
                    else
                        I(floor(zn(i))+1,floor(xn(i))+1) = round((double(I(floor(zn(i))+1,floor(xn(i))+1)) + g(i))/2);
                    end
                end
            end
            if inclined %If drop is showed inclined
                %Find drop apex
                indx_apex = find(x_ly_sim_pxs == 0);
                %Divide drop profile and determine the right side
                x_plot_r = x_plot(indx_apex-1000:end);
                z_plot_r = z_plot(indx_apex-1000:end);
                [x_plot_r_rot,z_plot_r_rot] = rotate_data(x_plot_r,z_plot_r,xc,zc,alpha_deg);
                x_plot_l = x_plot(1:indx_apex);
                z_plot_l = z_plot(1:indx_apex);
                [x_plot_l_rot,z_plot_l_rot] = rotate_data(x_plot_l,z_plot_l,xc,zc,alpha_deg);
                %Filling image inside drop
                for n = 1:length(z_plot_r_rot)
                    for m = round(x_plot_r_rot(n)):-1:1
                        if I(round(z_plot_r_rot(n)),m) == 0
                            I(round(z_plot_r_rot(n)),m) = g1;
                        elseif I(round(z_plot_r_rot(n)),m) == g2
                            break;
                        end
                    end
                end
                for n = 1:length(z_plot_l_rot)
                    for m = round(x_plot_l_rot(n)):res_h
                        if I(round(z_plot_l_rot(n)),m) == 0
                            I(round(z_plot_l_rot(n)),m) = g1;
                        elseif I(round(z_plot_l_rot(n)),m) == g2
                            break;
                        end
                    end
                end
                %Filling solid substrate
                for m=1:res_h %Search from bottom to top
                    if I(res_v,m) == g2
                        break;
                    else
                        for n =res_v:-1:1
                            if I(n,m) == 0
                                I(n,m) = g1;
                            elseif I(n,m) == g2
                                break;
                            end
                        end
                    end
                end
                for n=res_v:-1:1 %Search from left to right
                    if I(n,1) == g2
                        break;
                    else
                        for m =1:res_h
                            if I(n,m) == 0
                                I(n,m) = g1;
                            elseif I(n,m) == g2
                                break;
                            end
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
            else %If drop is not showed inclined
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
                %- Filling solid substrate and drop
                for i = 1:res_v %Search in all the image
                    for j = 1:res_h
                        if I(i,j) == 0
                            I(i,j) = g1;
                        end
                    end
                end
            end
            if option == 1
                %Save images
                imagename = strcat('Grayscale_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_CA(°)_',num2str(CA_sim(k),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_b(cm-1)_',num2str(b(k),'%.3f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_CA_max(°)_',num2str(CA_adv(k),'%.2f'),'_CA_min(°)_',num2str(CA_rec(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
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
                        imagename = strcat('Noise_Grayscale_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_CA(°)_',num2str(CA_sim(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_noise_',num2str(round(noise(q))),'_n_images_',num2str(p)); %Image filename
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
                        imagename = strcat('Gaussian_Grayscale_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_CA(°)_',num2str(CA_sim(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_g_noise_mean_',num2str(g_noise_mean),'_g_noise_var_',num2str(g_noise_var(q),'%.3f'),'_n_images_',num2str(p)); %Image filename
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
                        imagename = strcat('Salt_Pepper_Noise_Grayscale_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_CA(°)_',num2str(CA_sim(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_salt_pepper_density_',num2str(sp_noise_density(q),'%.2f'),'_n_images_',num2str(p)); %Image filename
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
                        imagename = strcat('Satellite_Grayscale_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_CA(°)_',num2str(CA_sim(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_d_sat(cm)_',num2str(d_sat_cm,'%.3f'),'_a_',num2str(a(p),'%.4f'),'_n_images_',num2str(m)); %Image filename
                        fullimagename = strcat(imagename,imFmt); %Image name + image format
                        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                        imwrite(Isat,fullfilename) %Save Image
                    end
                end
            end
        %--------------------------------------------------------------
        %--------------- GENERATION OF 3D DROP SURFACE ----------------
        %--------------------------------------------------------------
        %{
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
        %}
        %--------------------------------------------------------------
        %------------ SIMULATION OF DIFFERENT EDGE WIDTHS -------------
        %--------------------------------------------------------------
        elseif option == 3 %SIMULATION OF BLUR IN LIQUID DROP AND SOLID SUBSTRATE PROFILES (EDGE WIDTH)
            for q = 1:length(W_pxs_array)
                %Creation of the inclined drop normal edge and normal profile
                x_plot = x_ly_sim_pxs + xc - (x_ly_sim_pxs(1)+(x_ly_sim_pxs(end)-x_ly_sim_pxs(1))/2);

                %Defining inclined drop profile position in the image
                if inclined
                    z_plot = z_ly_sim_pxs + res_v/2 - (max(abs(z_ly_sim_pxs))/2);
                    z_subst = z_plot(end);
                else
                    z_plot = z_ly_sim_pxs + res_v - (max(abs(z_ly_sim_pxs))+hsubstrato);
                    z_subst = res_v-hsubstrato;
                end
                [xn_drop,zn_drop,g_drop] = grayscale_edge(x_plot,z_plot,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to inclined drop profile
                g_drop = wrev(g_drop);

                %Creation of solid subratre edge normal profile
                %x_sol_l = linspace(-res_h/2,x_plot(1),ceil(x_plot(1)+res_h/2));
                x_sol_l = linspace(-res_h/2,x_plot(1),10*ceil(x_plot(1)+res_h/2));
                z_sol_l = ones(1,length(x_sol_l))*(z_subst);
                %x_sol_r = linspace(x_plot(end),3*res_h/2,ceil((3*res_h/2)-x_plot(end)));
                x_sol_r = linspace(x_plot(end),3*res_h/2,10*ceil((3*res_h/2)-x_plot(end)));
                z_sol_r = ones(1,length(x_sol_r))*(z_subst);
                [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the left portion of the solid substrate surface
                g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
                [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the right portion of the solid substrate surface
                g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
                xn_sol = [xn_sol_l xn_sol_r];
                zn_sol = [zn_sol_l zn_sol_r];
                g_sol = [g_sol_l g_sol_r];

                %Trim drop normal edge according to solid substrate
                %ind_z = find(zn_drop <= (res_v-hsubstrato));
                ind_z = find(zn_drop <= z_subst);
                xn_drop = xn_drop(ind_z);
                zn_drop = zn_drop(ind_z);
                g_drop = g_drop(ind_z);

                %Combine solid and drop edge and normal points
                xn = [xn_sol xn_drop];
                zn = [zn_sol zn_drop];
                g = [g_sol g_drop];

                %Rotation of edge and normal points
                if inclined
                    [xn,zn] = rotate_data(xn,zn,xc,zc,alpha_deg);
                end

                %Determination of  triple contact points
                if inclined
                    [x_plot_rot,z_plot_rot] = rotate_data(x_plot,z_plot,xc,zc,alpha_deg);
                    Xcp_left(k) = x_plot_rot(1);
                    Ycp_left(k) = z_plot_rot(1);
                    Xcp_right(k) = x_plot_rot(end);
                    Ycp_right(k) = z_plot_rot(end);
                else
                    Xcp_left(k) = x_plot(1);
                    Ycp_left(k) = z_plot(1);
                    Xcp_right(k) = x_plot(end);
                    Ycp_right(k) = z_plot(end);
                end

                %Generation of synthetical illuminated (grayscale) pendant drop images
                M = zeros(res_v,res_h); %Array with image resolution size
                I = uint8(M); %Image creation
                %- Plotting solid substrate and drop profile edge
                for i=1:length(g)
                    if (xn(i)>=0) && (xn(i)<res_h) && (zn(i)>=0) && (zn(i)<res_v)
                        if double(I(floor(zn(i))+1,floor(xn(i))+1)) == 0
                            I(floor(zn(i))+1,floor(xn(i))+1) = round(g(i));
                        else
                            I(floor(zn(i))+1,floor(xn(i))+1) = round((double(I(floor(zn(i))+1,floor(xn(i))+1)) + g(i))/2);
                        end
                    end
                end
                if inclined %If drop is showed inclined
                    %Find drop apex
                    indx_apex = find(x_ly_sim_pxs == 0);
                    %Divide drop profile and determine the right side
                    x_plot_r = x_plot(indx_apex-1000:end);
                    z_plot_r = z_plot(indx_apex-1000:end);
                    [x_plot_r_rot,z_plot_r_rot] = rotate_data(x_plot_r,z_plot_r,xc,zc,alpha_deg);
                    x_plot_l = x_plot(1:indx_apex);
                    z_plot_l = z_plot(1:indx_apex);
                    [x_plot_l_rot,z_plot_l_rot] = rotate_data(x_plot_l,z_plot_l,xc,zc,alpha_deg);
                    %Filling image inside drop
                    for n = 1:length(z_plot_r_rot)
                        for m = round(x_plot_r_rot(n)):-1:1
                            if I(round(z_plot_r_rot(n)),m) == 0
                                I(round(z_plot_r_rot(n)),m) = g1;
                            elseif I(round(z_plot_r_rot(n)),m) == g2
                                break;
                            end
                        end
                    end
                    for n = 1:length(z_plot_l_rot)
                        for m = round(x_plot_l_rot(n)):res_h
                            if I(round(z_plot_l_rot(n)),m) == 0
                                I(round(z_plot_l_rot(n)),m) = g1;
                            elseif I(round(z_plot_l_rot(n)),m) == g2
                                break;
                            end
                        end
                    end
                    %Filling solid substrate
                    for m=1:res_h %Search from bottom to top
                        if I(res_v,m) == g2
                            break;
                        else
                            for n =res_v:-1:1
                                if I(n,m) == 0
                                    I(n,m) = g1;
                                elseif I(n,m) == g2
                                    break;
                                end
                            end
                        end
                    end
                    for n=res_v:-1:1 %Search from left to right
                        if I(n,1) == g2
                            break;
                        else
                            for m =1:res_h
                                if I(n,m) == 0
                                    I(n,m) = g1;
                                elseif I(n,m) == g2
                                    break;
                                end
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
                else %If drop is not showed inclined
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
                    %- Filling solid substrate and drop
                    for i = 1:res_v %Search in all the image
                        for j = 1:res_h
                            if I(i,j) == 0
                                I(i,j) = g1;
                            end
                        end
                    end
                end
                %Save images
                imagename = strcat('Edge_Width_Grayscale_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_CA(°)_',num2str(CA_sim(k),'%.2f'),'_CA_max(°)_',num2str(CA_adv(k),'%.2f'),'_CA_min(°)_',num2str(CA_rec(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs_array(q),'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
                fullimagename = strcat(imagename,imFmt); %Image name + image format
                fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                imwrite(I,fullfilename) %Save Image
            end
        %--------------------------------------------------------------
        %------- SIMULATION OF CAMERA VERTICAL MISALIGNMENT -----------
        %--------------------------------------------------------------
        elseif option == 4 % SIMULATION OF CAMERA VERTICAL MISALIGNMENT
            for q = 1:length(phi)
                %Creation of the inclined drop normal edge and normal profile
                x_plot = x_ly_sim_pxs + xc - (x_ly_sim_pxs(1)+(x_ly_sim_pxs(end)-x_ly_sim_pxs(1))/2);

                %Defining inclined drop profile position in the image
                if inclined
                    z_plot = z_ly_sim_pxs + res_v/2 - (max(abs(z_ly_sim_pxs))/2);
                    z_subst = z_plot(end);
                else
                    z_plot = z_ly_sim_pxs + res_v - (max(abs(z_ly_sim_pxs))+hsubstrato);
                    z_subst = res_v-hsubstrato;
                end
                [xn_drop,zn_drop,g_drop] = grayscale_edge(x_plot,z_plot,n_normal,g1,g2,W_pxs); %Calculation of the normal points to inclined drop profile
                g_drop = wrev(g_drop);

                %Creation of solid subratre edge normal profile
                %x_sol_l = linspace(-res_h/2,x_plot(1),ceil(x_plot(1)+res_h/2));
                x_sol_l = linspace(-res_h/2,x_plot(1),10*ceil(x_plot(1)+res_h/2));
                z_sol_l = ones(1,length(x_sol_l))*(z_subst);
                %x_sol_r = linspace(x_plot(end),3*res_h/2,ceil((3*res_h/2)-x_plot(end)));
                x_sol_r = linspace(x_plot(end),3*res_h/2,10*ceil((3*res_h/2)-x_plot(end)));
                z_sol_r = ones(1,length(x_sol_r))*(z_subst);
                [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the solid substrate surface
                g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
                [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the solid substrate surface
                g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
                xn_sol = [xn_sol_l xn_sol_r];
                zn_sol = [zn_sol_l zn_sol_r];
                g_sol = [g_sol_l g_sol_r];

                %Trim drop normal edge according to solid substrate
                %ind_z = find(zn_drop <= (res_v-hsubstrato));
                ind_z = find(zn_drop <= z_subst);
                xn_drop = xn_drop(ind_z);
                zn_drop = zn_drop(ind_z);
                g_drop = g_drop(ind_z);

                %Combine solid and drop edge and normal points
                xn = [xn_sol xn_drop];
                zn = [zn_sol zn_drop];
                g = [g_sol g_drop];

                %Rotation of edge and normal points
                if inclined
                    [xn,zn] = rotate_data(xn,zn,xc,zc,alpha_deg);
                end

                %Determination of  triple contact points
                if inclined
                    [x_plot_rot,z_plot_rot] = rotate_data(x_plot,z_plot,xc,zc,alpha_deg);
                    Xcp_left(k) = x_plot_rot(1);
                    Ycp_left(k) = z_plot_rot(1);
                    Xcp_right(k) = x_plot_rot(end);
                    Ycp_right(k) = z_plot_rot(end);
                else
                    Xcp_left(k) = x_plot(1);
                    Ycp_left(k) = z_plot(1);
                    Xcp_right(k) = x_plot(end);
                    Ycp_right(k) = z_plot(end);
                end

                %Rotation of data according to vertical misalignment
                if inclined
                    [x_plot_rot,z_plot_rot] = rotate_data(x_plot_rot,z_plot_rot,xc_rot,zc_rot,phi(q));
                else
                    [x_plot_rot,z_plot_rot] = rotate_data(x_plot,z_plot,xc_rot,zc_rot,phi(q));
                end
                [xn_rot,zn_rot] = rotate_data(xn,zn,xc_rot,zc_rot,phi(q)); %all data

                %Generation of synthetical illuminated (grayscale) pendant drop images
                M = zeros(res_v,res_h); %Array with image resolution size
                I = uint8(M); %Image creation
                %- Plotting solid substrate and drop profile edge
                for i=1:length(g)
                    if (xn_rot(i)>=0) && (xn_rot(i)<res_h) && (zn_rot(i)>=0) && (zn_rot(i)<res_v)
                        if double(I(floor(zn_rot(i))+1,floor(xn_rot(i))+1)) == 0
                            I(floor(zn_rot(i))+1,floor(xn_rot(i))+1) = round(g(i));
                        else
                            I(floor(zn_rot(i))+1,floor(xn_rot(i))+1) = round((double(I(floor(zn_rot(i))+1,floor(xn_rot(i))+1)) + g(i))/2);
                        end
                    end
                end
                %Finding drop apex and dividing drop profile
                indx_z = find(z_plot_rot == min(z_plot_rot),1,'first'); %Find drop apex
                x_plot_l_rot = x_plot_rot(1:indx_z);
                z_plot_l_rot = z_plot_rot(1:indx_z);
                x_plot_r_rot = x_plot_rot(indx_z:end);
                z_plot_r_rot = z_plot_rot(indx_z:end);
                if (alpha_deg + phi(q)) < 0 %Rotação anti-horária
                    %Filling image inside drop
                    for n = 1:length(z_plot_l_rot)
                        for m = round(x_plot_l_rot(n)):round(res_h)
                            if I(round(z_plot_l_rot(n)),m) == 0
                                I(round(z_plot_l_rot(n)),m) = g1;
                            elseif I(round(z_plot_l_rot(n)),m) == g2
                                break;
                            end
                        end
                    end
                else %Rotação horária
                    %Filling image inside drop
                    for n = 1:length(z_plot_r_rot)
                        for m = round(x_plot_r_rot(n)):-1:1
                            if I(round(z_plot_r_rot(n)),m) == 0
                                I(round(z_plot_r_rot(n)),m) = g1;
                            elseif I(round(z_plot_r_rot(n)),m) == g2
                                break;
                            end
                        end
                    end
                end
                %Filling solid substrate
                for m=1:res_h %Search from bottom to top
                    if I(res_v,m) == g2
                        break;
                    else
                        for n =res_v:-1:1
                            if I(n,m) == 0
                                I(n,m) = g1;
                            elseif I(n,m) == g2
                                break;
                            end
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
                imagename = strcat('Vertical_Misalignment_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_CA(°)_',num2str(CA_sim(k),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_CA_max(°)_',num2str(CA_adv(k),'%.2f'),'_CA_min(°)_',num2str(CA_rec(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_phi(°)_',num2str(phi(q),'%.3f')); %Image filename
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
                    %Creation of the inclined drop normal edge and normal profile
                    x_plot = x_ly_sim_pxs + xc - (x_ly_sim_pxs(1)+(x_ly_sim_pxs(end)-x_ly_sim_pxs(1))/2);

                    %Defining inclined drop profile position in the image
                    if inclined
                        z_plot = z_ly_sim_pxs + res_v/2 - (max(abs(z_ly_sim_pxs))/2);
                        z_subst = z_plot(end);
                    else
                        z_plot = z_ly_sim_pxs + res_v - (max(abs(z_ly_sim_pxs))+hsubstrato);
                        z_subst = res_v-hsubstrato;
                    end
                    [xn_drop,zn_drop,g_drop] = grayscale_edge(x_plot,z_plot,n_normal,g1,g2,W_pxs); %Calculation of the normal points to inclined drop profile
                    g_drop = wrev(g_drop);

                    %Random perturbation of the drop profile
                    [xn_drop,zn_drop,xn_drop_edge,zn_drop_edge] = normal_perturbation_edge(xn_drop,zn_drop,pert(p),n_normal);

                    %Creation of solid subratre edge normal profile
                    %x_sol_l = linspace(-res_h/2,x_plot(1),ceil(x_plot(1)+res_h/2));
                    x_sol_l = linspace(-res_h/2,x_plot(1),10*ceil(x_plot(1)+res_h/2));
                    z_sol_l = ones(1,length(x_sol_l))*(z_subst);
                    %x_sol_r = linspace(x_plot(end),3*res_h/2,ceil((3*res_h/2)-x_plot(end)));
                    x_sol_r = linspace(x_plot(end),3*res_h/2,10*ceil((3*res_h/2)-x_plot(end)));
                    z_sol_r = ones(1,length(x_sol_r))*(z_subst);
                    [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the solid substrate surface
                    g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
                    [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the solid substrate surface
                    g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
                    xn_sol = [xn_sol_l xn_sol_r];
                    zn_sol = [zn_sol_l zn_sol_r];
                    g_sol = [g_sol_l g_sol_r];

                    %Trim drop normal edge according to solid substrate
                    %ind_z = find(zn_drop <= (res_v-hsubstrato));
                    ind_z = find(zn_drop <= z_subst);
                    xn_drop = xn_drop(ind_z);
                    zn_drop = zn_drop(ind_z);
                    g_drop = g_drop(ind_z);

                    %Combine solid and drop edge and normal points
                    xn = [xn_sol xn_drop];
                    zn = [zn_sol zn_drop];
                    g = [g_sol g_drop];

                    %Rotation of edge and normal points
                    if inclined
                        [xn,zn] = rotate_data(xn,zn,xc,zc,alpha_deg);
                    end

                    %Determination of  triple contact points
                    if inclined
                        [x_plot_rot,z_plot_rot] = rotate_data(x_plot,z_plot,xc,zc,alpha_deg);
                        Xcp_left(k) = x_plot_rot(1);
                        Ycp_left(k) = z_plot_rot(1);
                        Xcp_right(k) = x_plot_rot(end);
                        Ycp_right(k) = z_plot_rot(end);
                    else
                        Xcp_left(k) = x_plot(1);
                        Ycp_left(k) = z_plot(1);
                        Xcp_right(k) = x_plot(end);
                        Ycp_right(k) = z_plot(end);
                    end

                    %Generation of synthetical illuminated (grayscale) pendant drop images
                    M = zeros(res_v,res_h); %Array with image resolution size
                    I = uint8(M); %Image creation
                    %- Plotting solid substrate and drop profile edge
                    for i=1:length(g)
                        if (xn(i)>=0) && (xn(i)<res_h) && (zn(i)>=0) && (zn(i)<res_v)
                            if double(I(floor(zn(i))+1,floor(xn(i))+1)) == 0
                                I(floor(zn(i))+1,floor(xn(i))+1) = round(g(i));
                            else
                                I(floor(zn(i))+1,floor(xn(i))+1) = round((double(I(floor(zn(i))+1,floor(xn(i))+1)) + g(i))/2);
                            end
                        end
                    end
                    if inclined %If drop is showed inclined
                        [xn_drop_edge_rot,zn_drop_edge_rot] = rotate_data(xn_drop_edge,zn_drop_edge,xc,zc,alpha_deg);
                        %Finding drop apex and dividing drop profile
                        indx_z = find(abs(zn_drop_edge_rot) == min(abs(zn_drop_edge_rot)),1,'first'); %Find drop apex
                        xn_drop_edge_rot_l = xn_drop_edge_rot(1:indx_z);
                        zn_drop_edge_rot_l = zn_drop_edge_rot(1:indx_z);
                        xn_drop_edge_rot_r = xn_drop_edge_rot(indx_z:end);
                        zn_drop_edge_rot_r = zn_drop_edge_rot(indx_z:end);

                        %Filling image inside drop
                        for n = 1:length(zn_drop_edge_rot_l)
                            for m = round(xn_drop_edge_rot_l(n)):round(res_h)
                                if I(round(zn_drop_edge_rot_l(n)),m) == 0
                                    I(round(zn_drop_edge_rot_l(n)),m) = g1;
                                elseif I(round(zn_drop_edge_rot_l(n)),m) == g2
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
                    %{
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
                    %}
                    %Filling solid substrate
                    for m=1:res_h %Search from bottom to top
                        if I(res_v,m) == g2
                            break;
                        else
                            for n =res_v:-1:1
                                if I(n,m) == 0
                                    I(n,m) = g1;
                                elseif I(n,m) == g2
                                    break;
                                end
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
                    imagename = strcat('Random_Edge_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_CA(°)_',num2str(CA_sim(k),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_CA_max(°)_',num2str(CA_adv(k),'%.2f'),'_CA_min(°)_',num2str(CA_rec(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_pert(pxs)_',num2str(pert(p),'%.2f'),'_image_',num2str(q)); %Image filename
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

                %Creation of the inclined drop normal edge and normal profile
                x_plot = x_ly_sim_pxs + xc - (x_ly_sim_pxs(1)+(x_ly_sim_pxs(end)-x_ly_sim_pxs(1))/2);

                %Defining inclined drop profile position in the image
                if inclined
                    z_plot = z_ly_sim_pxs + res_v/2 - (max(abs(z_ly_sim_pxs))/2);
                    z_subst = z_plot(end);
                else
                    z_plot = z_ly_sim_pxs + res_v - (max(abs(z_ly_sim_pxs))+hsubstrato);
                    z_subst = res_v-hsubstrato;
                end
                [xn_drop,zn_drop,g_drop] = grayscale_edge(x_plot,z_plot,n_normal,g1,g2,W_pxs); %Calculation of the normal points to inclined drop profile
                g_drop = wrev(g_drop);

                %Creation of solid subratre edge normal profile
                %x_sol_l = linspace(-res_h/2,x_plot(1),ceil(x_plot(1)+res_h/2));
                x_sol_l = linspace(-res_h/2,x_plot(1),10*ceil(x_plot(1)+res_h/2));
                z_sol_l = ones(1,length(x_sol_l))*(z_subst);
                %x_sol_r = linspace(x_plot(end),3*res_h/2,ceil((3*res_h/2)-x_plot(end)));
                x_sol_r = linspace(x_plot(end),3*res_h/2,10*ceil((3*res_h/2)-x_plot(end)));
                z_sol_r = ones(1,length(x_sol_r))*(z_subst);
                [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the solid substrate surface
                g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
                [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the solid substrate surface
                g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
                xn_sol = [xn_sol_l xn_sol_r];
                zn_sol = [zn_sol_l zn_sol_r];
                g_sol = [g_sol_l g_sol_r];

                %Trim drop normal edge according to solid substrate
                %ind_z = find(zn_drop <= (res_v-hsubstrato));
                ind_z = find(zn_drop <= z_subst);
                xn_drop = xn_drop(ind_z);
                zn_drop = zn_drop(ind_z);
                g_drop = g_drop(ind_z);

                %Combine solid and drop edge and normal points
                xn = [xn_sol xn_drop];
                zn = [zn_sol zn_drop];
                g = [g_sol g_drop];

                %Rotation of edge and normal points
                if inclined
                    [xn,zn] = rotate_data(xn,zn,xc,zc,alpha_deg);
                end

                %Determination of  triple contact points
                if inclined
                    [x_plot_rot,z_plot_rot] = rotate_data(x_plot,z_plot,xc,zc,alpha_deg);
                    Xcp_left(k) = x_plot_rot(1);
                    Ycp_left(k) = z_plot_rot(1);
                    Xcp_right(k) = x_plot_rot(end);
                    Ycp_right(k) = z_plot_rot(end);
                else
                    Xcp_left(k) = x_plot(1);
                    Ycp_left(k) = z_plot(1);
                    Xcp_right(k) = x_plot(end);
                    Ycp_right(k) = z_plot(end);
                end

                %Generation of synthetical illuminated (grayscale) pendant drop images
                M = zeros(res_v,res_h); %Array with image resolution size
                I = uint8(M); %Image creation
                %- Plotting solid substrate and drop profile edge
                for i=1:length(g)
                    if (xn(i)>=0) && (xn(i)<res_h) && (zn(i)>=0) && (zn(i)<res_v)
                        if double(I(floor(zn(i))+1,floor(xn(i))+1)) == 0
                            I(floor(zn(i))+1,floor(xn(i))+1) = round(g(i));
                        else
                            I(floor(zn(i))+1,floor(xn(i))+1) = round((double(I(floor(zn(i))+1,floor(xn(i))+1)) + g(i))/2);
                        end
                    end
                end
                if inclined %If drop is showed inclined
                    %Find drop apex
                    indx_apex = find(x_ly_sim_pxs == 0);
                    %Divide drop profile and determine the right side
                    x_plot_r = x_plot(indx_apex-1000:end);
                    z_plot_r = z_plot(indx_apex-1000:end);
                    [x_plot_r_rot,z_plot_r_rot] = rotate_data(x_plot_r,z_plot_r,xc,zc,alpha_deg);
                    x_plot_l = x_plot(1:indx_apex);
                    z_plot_l = z_plot(1:indx_apex);
                    [x_plot_l_rot,z_plot_l_rot] = rotate_data(x_plot_l,z_plot_l,xc,zc,alpha_deg);
                    %Filling image inside drop
                    for i = 1:length(z_plot_r_rot)
                        for j = round(x_plot_r_rot(i)):-1:1
                            if I(round(z_plot_r_rot(i)),j) == 0
                                I(round(z_plot_r_rot(i)),j) = g1;
                            elseif I(round(z_plot_r_rot(i)),j) == g2
                                break;
                            end
                        end
                    end
                    for i = 1:length(z_plot_l_rot)
                        for j = round(x_plot_l_rot(i)):res_h
                            if I(round(z_plot_l_rot(i)),j) == 0
                                I(round(z_plot_l_rot(i)),j) = g1;
                            elseif I(round(z_plot_l_rot(i)),j) == g2
                                break;
                            end
                        end
                    end
                    %Filling solid substrate
                    for i=1:res_h %Search from bottom to top
                        if I(res_v,i) == g2
                            break;
                        else
                            for j =res_v:-1:1
                                if I(j,i) == 0
                                    I(j,i) = g1;
                                elseif I(j,i) == g2
                                    break;
                                end
                            end
                        end
                    end
                    for i=res_v:-1:1 %Search from left to right
                        if I(i,1) == g2
                            break;
                        else
                            for j =1:res_h
                                if I(i,j) == 0
                                    I(i,j) = g1;
                                elseif I(i,j) == g2
                                    break;
                                end
                            end
                        end
                    end
                    %Filling background
                    for i = 1:res_v %Search in all the image
                        for j = 1:res_h
                            if I(i,j) == 0
                                I(i,j) = g2;
                            end
                        end
                    end
                else %If drop is not showed inclined
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
                    %- Filling solid substrate and drop
                    for i = 1:res_v %Search in all the image
                        for j = 1:res_h
                            if I(i,j) == 0
                                I(i,j) = g1;
                            end
                        end
                    end
                end
                %Save images
                imagename = strcat('Grayscale_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_CA(°)_',num2str(CA_sim(k),'%.2f'),'_CA_max(°)_',num2str(CA_adv(k),'%.2f'),'_CA_min(°)_',num2str(CA_rec(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
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
                %Creation of the inclined drop normal edge and normal profile
                x_plot = x_ly_sim_pxs + xc - (x_ly_sim_pxs(1)+(x_ly_sim_pxs(end)-x_ly_sim_pxs(1))/2);
                if inclined
                    z_plot = z_ly_sim_pxs + res_v/2 - (max(abs(z_ly_sim_pxs))/2);
                    z_subst = z_plot(end);
                    [x_plot_rot,z_plot_rot] = rotate_data(x_plot,z_plot,xc,zc,alpha_deg);
                    [xn_drop,zn_drop,g_drop] = grayscale_edge_nonuniform_illumination(x_plot_rot,z_plot_rot,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the left portion of the solid substrate surface
                else
                    z_plot = z_ly_sim_pxs + res_v - (max(abs(z_ly_sim_pxs))+hsubstrato);
                    z_subst = res_v-hsubstrato;
                    [xn_drop,zn_drop,g_drop] = grayscale_edge_nonuniform_illumination(x_plot,z_plot,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the left portion of the solid substrate surface
                end
                g_drop = wrev(g_drop);

                %Creation of solid subratre edge normal profile
                %x_sol_l = linspace(-res_h/2,x_plot(1),ceil(x_plot(1)+res_h/2));
                x_sol_l = linspace(-res_h/2,x_plot(1),10*ceil(x_plot(1)+res_h/2));
                z_sol_l = ones(1,length(x_sol_l))*(z_subst);
                %x_sol_r = linspace(x_plot(end),3*res_h/2,ceil((3*res_h/2)-x_plot(end)));
                x_sol_r = linspace(x_plot(end),3*res_h/2,10*ceil((3*res_h/2)-x_plot(end)));
                z_sol_r = ones(1,length(x_sol_r))*(z_subst);
                if inclined
                    [x_sol_l,z_sol_l] = rotate_data(x_sol_l,z_sol_l,xc,zc,alpha_deg);
                    [x_sol_r,z_sol_r] = rotate_data(x_sol_r,z_sol_r,xc,zc,alpha_deg);
                end
                x_sol_l = wrev(x_sol_l);
                z_sol_l = wrev(z_sol_l);
                [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge_nonuniform_illumination(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the left portion of the solid substrate surface
                [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge_nonuniform_illumination(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the right portion of the solid substrate surface
                for i = 1:(length(g_sol_r)/n_normal) %Correcting g_sol_r order
                    prov = g_sol_r(n_normal*(i-1)+1:n_normal*i);
                    for j = 1:floor(n_normal/2)
                        g_sol_r((n_normal*(i-1))+j) = prov((n_normal)-(j-1));
                        g_sol_r((n_normal*i)-(j-1)) = prov(j);
                    end
                end
                xn_sol = [wrev(xn_sol_l) xn_sol_r];
                zn_sol = [wrev(zn_sol_l) zn_sol_r];
                g_sol = [wrev(g_sol_l) g_sol_r];

                %Trim drop normal edge according to solid substrate
                if inclined
                    if alpha_deg < 90
                        p_coef = polyfit(x_sol_l,z_sol_l,1);
                        indx_subst = find(zn_drop < (p_coef(1)*xn_drop+p_coef(2)));
                        xn_drop = xn_drop(indx_subst);
                        zn_drop = zn_drop(indx_subst);
                        g_drop = g_drop(indx_subst);
                    else
                        indx_subst = find(xn_drop > x_sol_r(1));
                        xn_drop = xn_drop(indx_subst);
                        zn_drop = zn_drop(indx_subst);
                        g_drop = g_drop(indx_subst);
                    end
                else
                    ind_subst = find(zn_drop <= (res_v-hsubstrato));
                    xn_drop = xn_drop(ind_subst);
                    zn_drop = zn_drop(ind_subst);
                    g_drop = g_drop(ind_subst);
                end

                %Combine solid and drop edge and normal points
                xn = [xn_sol xn_drop];
                zn = [zn_sol zn_drop];
                g = [g_sol g_drop];

                %Generation of synthetical illuminated (grayscale) pendant drop images
                M = zeros(res_v,res_h); %Matriz com a resolução indicada
                I = uint8(M); %Criação da imagem
                %- Plotting solid substrate and drop profile edge
                for i=1:length(g)
                    if (xn(i)>=0) && (xn(i)<res_h) && (zn(i)>=0) && (zn(i)<res_v)
                        if double(I(floor(zn(i))+1,floor(xn(i))+1)) == 0
                            I(floor(zn(i))+1,floor(xn(i))+1) = round(g(i));
                        else
                            I(floor(zn(i))+1,floor(xn(i))+1) = round((double(I(floor(zn(i))+1,floor(xn(i))+1)) + g(i))/2);
                        end
                    end
                end
                if inclined %If drop is showed inclined
                    if alpha_deg < 90
                        %- 1st methodd - fit using cubic interpolation
                        %{
                %Find drop apex
                indx_apex = find(x_ly_sim_pxs == 0);
                %Divide drop profile and determine the left and right side
                x_plot_l = x_plot(1:indx_apex);
                z_plot_l = z_plot(1:indx_apex);
                [x_plot_l_rot,z_plot_l_rot] = rotate_data(x_plot_l,z_plot_l,xc,zc,alpha_deg);
                x_plot_r = x_plot(indx_apex+1:end);
                z_plot_r = z_plot(indx_apex+1:end);
                [x_plot_r_rot,z_plot_r_rot] = rotate_data(x_plot_r,z_plot_r,xc,zc,alpha_deg);
                %Return drop profile data but with no repetitions
                x_prof_l = [wrev(x_sol_l) x_plot_l_rot];
                x_prof_l = unique(x_prof_l,'stable');
                z_prof_l = [wrev(z_sol_l) z_plot_l_rot];
                z_prof_l = unique(z_prof_l,'stable');
                x_prof_r = [x_plot_r_rot x_sol_r];
                x_prof_r = unique(x_prof_r,'stable');
                z_prof_r = [z_plot_r_rot z_sol_r];
                z_prof_r = unique(z_prof_r,'stable');
                x_prof = [x_prof_l x_prof_r];
                z_prof = [z_prof_l z_prof_r];
                %Fit a pieciwise cubic interpolation on data
                f_prof = fit(x_prof',z_prof','cubicinterp');
            
                %Filling image drop, substrate and background
                dim_j = linspace(1,res_h,res_h);
                icalc = f_prof(dim_j);
                for i=1:res_v
                    for j=1:res_h
                        if i<icalc(j) && I(i,j) == 0
                            r = sqrt((i-iref)^2+(j-jref)^2);
                            I(i,j) = round(g2 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                        elseif i>icalc(j) && I(i,j) == 0
                            r = sqrt((i-iref)^2+(j-jref)^2);
                            I(i,j) = round(g1 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                        end
                    end
                end
                        %}
                        %- 2nd method - filling each part separately
                        %Find drop apex
                        indx_apex = find(x_ly_sim_pxs == 0);
                        %Divide drop profile and determine the left and right side
                        x_plot_l = x_plot(1:indx_apex);
                        z_plot_l = z_plot(1:indx_apex);
                        [x_plot_l_rot,z_plot_l_rot] = rotate_data(x_plot_l,z_plot_l,xc,zc,alpha_deg);
                        x_plot_r = x_plot(indx_apex+1:end);
                        z_plot_r = z_plot(indx_apex+1:end);
                        [x_plot_r_rot,z_plot_r_rot] = rotate_data(x_plot_r,z_plot_r,xc,zc,alpha_deg);
                        %Filling image inside drop
                        for i = 1:length(x_plot_l_rot) %Left profile - from top to bottom
                            for j = round(z_plot_l_rot(i)):res_v
                                if j == round(z_plot_r_rot(end))
                                    break;
                                elseif I(j,round(x_plot_l_rot(i))) == 0
                                    r = sqrt((j-iref)^2+(x_plot_l_rot(i)-jref)^2);
                                    I(j,round(x_plot_l_rot(i))) = round(g1 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                elseif I(j+2,round(x_plot_l_rot(i))) >= g2
                                    break;
                                end
                            end
                        end
                        for i = 1:length(x_plot_r_rot) %Right profile - from top to bottom
                            for j = round(z_plot_r_rot(i)):res_v
                                if j == round(z_plot_r_rot(end))
                                    break;
                                elseif I(j,round(x_plot_r_rot(i))) == 0
                                    r = sqrt((j-iref)^2+(x_plot_r_rot(i)-jref)^2);
                                    I(j,round(x_plot_r_rot(i))) = round(g1 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                elseif I(j+2,round(x_plot_r_rot(i))) >= g2
                                    break;
                                end
                            end
                        end
                        %Filling backgorund
                        gmax_l = max(I(1,:)); %Search from top to bottom(right to the left)
                        for j=res_h:-1:1
                            if I(1,j) == gmax_l && gmax_l~=0
                                break;
                            end
                            gmax = max(I(:,j));
                            for i=1:res_v
                                if I(i,j) == gmax && gmax~=0
                                    break;
                                elseif I(i,j) == 0
                                    r = sqrt((i-iref)^2+(j-jref)^2);
                                    I(i,j) = round(g2 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                end
                            end
                        end
                        for i = 1:length(z_plot_r_rot) %Search from right drop profile - from left to right
                            for j = round(x_plot_r_rot(i)):res_h
                                if I(round(z_plot_r_rot(i)),j) == 0
                                    r = sqrt((z_plot_r_rot(i)-iref)^2+(j-jref)^2);
                                    I(round(z_plot_r_rot(i)),j) = round(g2 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                end
                            end
                        end
                        for i = 1:length(z_sol_r) %%Search from right substrate - from left to right
                            if z_sol_r(i) > res_v
                                break;
                            else
                                for j = round(x_sol_r(i)):res_h
                                    if I(round(z_sol_r(i)),j) == 0
                                        r = sqrt((z_sol_r(i)-iref)^2+(j-jref)^2);
                                        I(round(z_sol_r(i)),j) = round(g2 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                        %elseif I(j+2,round(x_plot_r_rot(i))) >= g2
                                        %    break;
                                    end
                                end
                            end
                        end
                        %Filling solid substrate
                        for j=1:res_h %Search from bottom to top (left to right)
                            gmax = max(I(:,j));
                            for i=res_v:-1:1
                                if I(i,j) == gmax && gmax~=0
                                    break;
                                elseif I(i,j) == 0
                                    r = sqrt((i-iref)^2+(j-jref)^2);
                                    I(i,j) = round(g1 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                end
                            end
                        end
                    else %alpha_deg = 90°
                        %Find drop apex
                        indx_apex = find(x_ly_sim_pxs == 0);
                        %Divide drop profile and determine the left and right side
                        x_plot_l = x_plot(1:indx_apex);
                        z_plot_l = z_plot(1:indx_apex);
                        [x_plot_l_rot,z_plot_l_rot] = rotate_data(x_plot_l,z_plot_l,xc,zc,alpha_deg);
                        x_plot_r = x_plot(indx_apex+1:end);
                        z_plot_r = z_plot(indx_apex+1:end);
                        [x_plot_r_rot,z_plot_r_rot] = rotate_data(x_plot_r,z_plot_r,xc,zc,alpha_deg);
                        %Filling image inside drop
                        for i = 1:length(x_plot_r_rot)
                            for j = round(z_plot_r_rot(i)):-1:1
                                if I(j,round(x_plot_r_rot(i))) == 0
                                    r = sqrt((j-iref)^2+(x_plot_r_rot(i)-jref)^2);
                                    I(j,round(x_plot_r_rot(i))) = round(g1 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                elseif I(j,round(x_plot_r_rot(i))) >= g2
                                    break;
                                end
                            end
                        end
                        %Filling solid substrate
                        for i=res_v:-1:1 %Search from left to right
                            gmax = max(I(i,:));
                            for j =1:res_h
                                if I(i,j) == 0
                                    r = sqrt((i-iref)^2+(j-jref)^2);
                                    I(i,j) = round(g1 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                elseif I(i,j) == gmax
                                    break;
                                end
                            end
                        end
                        %Filling background
                        for i=res_v:-1:1 %Search from right to left
                            gmax = max(I(i,:));
                            for j =res_h:-1:1
                                if I(i,j) == 0
                                    r = sqrt((i-iref)^2+(j-jref)^2);
                                    I(i,j) = round(g2 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                elseif I(i,j) == gmax
                                    break;
                                end
                            end
                        end
                    end
                else %Inclined = 0
                    %- Filling background
                    for i = 1:(res_v-hsubstrato) %Search from left to right
                        for j = 1:res_h
                            if I(i,j) == 0
                                r = sqrt((i-iref)^2+(j-jref)^2);
                                I(i,j) = round(g2 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
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
                            else
                                break;
                            end
                        end
                    end
                    %- Filling solid substrate and drop
                    for i = 1:res_v %Search in all the image
                        for j = 1:res_h
                            if I(i,j) == 0
                                r = sqrt((i-iref)^2+(j-jref)^2);
                                I(i,j) = round(g1 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                            end
                        end
                    end
                end
                %Save images
                imagename = strcat('Non_uniform_Grayscale_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_CA(°)_',num2str(CA_sim(k),'%.2f'),'_CA_max(°)_',num2str(CA_adv(k),'%.2f'),'_CA_min(°)_',num2str(CA_rec(k),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_deltag_',num2str(round(deltag(n)))); %Image filename
                fullimagename = strcat(imagename,imFmt); %Image name + image format
                fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                imwrite(I,fullfilename) %Save Image
            end
        end     
    end
    r0 = 1./b; %radius at the apex in cm
    %Write excel file with important pendant drop profile properties
    varnames = {'b[cm-1]','r0[cm]','CA[°]','CA_max[°]','CA_min[°]','Xcp_left[pxs]','Ycp_left[pxs]','Xcp_right[pxs]','Ycp_right[pxs]','d_w[cm]','h[cm]','A[cm2]','V[uL]'}; %Variables name in the excel file
    T = table(b,r0,CA_sim,CA_adv,CA_rec,Xcp_left,Ycp_left,Xcp_right,Ycp_right,d_w_cm,h_cm,A_cm2,V_uL,'VariableNames',varnames); %Build a table with important drop properties
    sheet_name = strcat('c = ',num2str(c),' cm^-2',' alpha = ',num2str(alpha_deg),' °'); %Excel sheet name
    savename_xls = strcat('Grayscale_Inclined_drop_',num2str(inclined),'(c(cm-2)_',num2str(c),'_CA(°)_',num2str(mean(CA_sim(k)),'%.2f'),'_alpha(°)_',num2str(alpha_deg),').xlsx'); %Excel file name
    fullsavename_xls = fullfile(save_directory,savename_xls);
    writetable(T,fullsavename_xls,'Sheet',sheet_name) %Write table to an excel file

    fprintf('Simulation completed! \n');
    pause(2)
    end
end
end