function simGrayPendDrop_b_c_rh()
%SIMGRAYPENDDROP_B_C_RH Simulation of grayscale images of pendant drop.
%Some source of errros can be simulated and the 3D drop surface can also be
%exported
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of grayscale images of pendant drop. Some source of errors, 
%   such as edge width, camera vertical misalignment, random perturbation 
%   in the drop profile, lack of contrast, non-uniform illumination, random
%   noise, Gaussian noise, impulse salt-pepper noise and satellite
%   droplets, can be simulated. The 3D drop surface can also be exported.
%   Depending on the simulated option, the program requires the user to
%   enter some parameters. As common input parameters, any option requires
%   the user to enter: the liquid capillary constant [cm-2](c), the minimum
%   curvature at the apex [cm-1](bmin), the maximum curvature at the apex 
%   [cm-1](bmax), the number of curvature at the apex levels simulated 
%   (nb), the needle diameter [mm](needleDiam), the possibility of 
%   simulating images with and without needle (needleON), image resolution 
%   [pxs] (res_h and res_v), image scale [pxs/mm] (scale), edge width 
%   (W_pxs), maximum and minimum image intensity levels (g1 and g2). The 
%   simulated pendant drops are centered in the image. The routine also 
%   allows the user to choose the directory where the images are going to 
%   be exported. An excel file (.xlsx) is created containing important 
%   properties of the pendant drops. During the computation of pendant drop
%   profiles errors may occur due to poor choices of input parameters. For 
%   certain combinations of curavture at the apex (bmin and bmax) and 
%   needle diameter, the pendant drop may not be simulated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while 1
    fprintf('------------------ PENDANT DROP - SIMULATION OF GRAYSCALE IMAGES --------------------\n');
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
        bmin = input("- Minimum curvature at the apex in cm-1 [7.2]: "); %Minimum curvature at the apex in cm-1
        if isempty(bmin)
            bmin = 7.2; %For water
        end
        while 1
            bmax = input("- Maximum curvature at the apex in cm-1 [9.7]: "); %Maximum curvature at the apex in cm-1;
            if isempty(bmax)
                bmax = 9.7; %For water
            end
            if bmax < bmin
                fprintf('Invalid maximum curavture at the apex. Please, enter a lower value. \n');
            else
                break;
            end
        end
        nb = input("- Number of curvature at the apex levels simulated [3]: "); %Number of simulated images between minimum and maximum curvature at the apex
        if isempty(nb)
            nb = 3;
        end
        b = linspace(bmin,bmax,nb);
        r0 = 1./b; %radius at the apex in cm
        fprintf('2) Drop configuration propertes: \n');
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
        %- Numerical integration options
        S_span=(0:.0001:8); % S_span is the step variable for ode45 solver
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
                diff_g_max = input("- Maximum difference between min.and max.image intensity levels [190]: ");
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
        fprintf('Selection of the directory to export files ... \n');
        save_directory = uigetdir('C:\'); %Save directory

        %- VARIABLES INITIALIZATION
        d_eq_cm = zeros(length(b),1); %Maximum drop diameter in cm
        h_eq_cm = zeros(length(b),1); %Distance from the apex of the drop to its maximum diameter (dop height) in cm
        d_s_cm = zeros(length(b),1); %Drop diameter from a distance equivalent to the maximum drop diameter from the apex of the drop in cm
        h_cm = zeros(length(b),1); %Distance from the apex of the drop to the needle tip in cm
        A_cm2 = zeros(length(b),1); %Drop surface area in cm²
        V_uL = zeros(length(b),1); %Drop volume in uL calculated concomitantly with the Laplace curve generation

        %-----------------------------------------------------------------
        %-------------- SIMULATION OF PENDANT DROP IMAGES -----------------
        %------------------------------------------------------------------
        fprintf("--------------------------------------------------------------------------------\n")
        fprintf("Simulation of grayscale pendant drop images... \n")
        for k=1:length(b)
            fprintf('Running b = %.4f cm-1 and rh = %.3f cm...\n',b(k),r_h_cm);
            %Determination of the pendant drop profile
            [x_ly_cm,z_ly_cm,x_ly_sim_cm,z_ly_sim_cm,indx_last,indx_second,indx_deq,indx_ds,d_eq_cm(k),h_eq_cm(k),h_cm(k),d_s_cm(k),A_cm2(k),V_uL(k)] = profile_pendant_drop_ly(b(k),c,r_h_cm,S_span,initialCo,tol);

            %Conversion to pixels
            x_ly_sim_pxs = x_ly_sim_cm*scale_cm;   % Drop width in pixels
            z_ly_sim_pxs = z_ly_sim_cm*scale_cm;   % Drop height in pixels
            %S_ly_sim_pxs = S_ly_sim_cm*scale_cm;   % Arc length in pixels

            %--------------------------------------------------------------
            %-------------- GENERATION OF GRAYSCALE IMAGES ----------------
            %--------------------------------------------------------------
            if option == 1 || option == 8 || option == 9 || option == 10 || option == 11 % Generation of grayscale images
                %Creation of needle and pendandt drop edge and normal profiles
                res_v = resolution(2);
                res_h = resolution(1);
                xc = round(res_h/2); %x coordinate of the drop center in pixels
                x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %x coordinates of the full drop profile in the final image
                if needle %Construction of pendant drop with needle
                    zc = round(res_v/2);
                    z_drop = -z_ly_sim_pxs+z_ly_sim_pxs(end)+zc-(z_ly_sim_pxs(end)/2);  %z coordinates of the left drop profile in the final image with needle
                    % Creation of needle drop edge and normal profiles
                    needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
                    pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                    pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                    z_needle_max = zc-(z_ly_sim_pxs(end)/2); %maximum z coordinate of the needle
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
                if option == 1
                    % - Save image
                    imagename = strcat('Grayscale_Pendant_drop_Needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_scale(pxs_cm)_',num2str(scale_cm),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_W(pxs)_',num2str(W_pxs,'.%4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
                    fullimagename = strcat(imagename,imFmt); %Image name + image format
                    fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                    imwrite(I,fullfilename) %Save Image
                %----------------------------------------------------------
                %------------------ RANDOM NOISE ADDITION ----------------
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
                            imagename = strcat('Noise_Grayscale_Pendant_drop_Needle_',num2str(needle),'_c(cm-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_noise_',num2str(round(noise(q))),'_n_images_',num2str(p)); %Image filename
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
                            imagename = strcat('Gaussian_Noise_Grayscale_Pendant_drop_Needle_',num2str(needle),'_c(cm-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_g_noise_mean_',num2str(g_noise_mean),'_g_noise_var_',num2str(g_noise_var(q),'%.3f'),'_n_images_',num2str(p)); %Image filename
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
                            imagename = strcat('Salt_Pepper_Noise_Pendant_drop_Needle_',num2str(needle),'_c(cm-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_salt_pepper_density_',num2str(sp_noise_density(q),'%.2f'),'_n_images_',num2str(p)); %Image filename
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
                                 [columnsInImage, rowsInImage] = meshgrid(1:res_h, 1:res_v);
                                 x_sat_center = randi([1 res_h],1);
                                 z_sat_center = randi([1 res_v],1);
                                 circlePixels = (rowsInImage-z_sat_center).^2 + (columnsInImage-x_sat_center).^2 <= (d_sat_pxs/2).^2;
                                 indx_circle = find(circlePixels);
                                 Isat(indx_circle) = g1; %Ploting satellite droplets on the image
                             end
                             %Save satellite droplet images
                             imagename = strcat('Satellite_Grayscale_Pendant_drop_Needle_',num2str(needle),'_c(cm-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_d_sat(cm)_',num2str(d_sat_cm,'%.3f'),'_a_',num2str(a(p),'%.4f'),'_n_images_',num2str(m)); %Image filename
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
                %- Create 3D drop profile
                %x_surf_lap =  x_ly_sim_pxs;
                %z_surf_lap = -z_ly_sim_pxs;
                if needle
                    needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
                    pos_x_needle = needle_diameter_pxs/2;
                    z_needle_max = z_ly_sim_pxs(end)+res_v; %maximum z coordinate of the needle
                    z_needle = linspace(z_ly_sim_pxs(end),z_needle_max,ceil(res_v)); %z coordinates of the needle profile
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
                filename_tif = strcat('3D_Surface_Pendant_drop_Needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_scale(pxs_cm)_',num2str(scale_cm),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_n_',num2str(n),'_Rsquare_',num2str(gof.rsquare,'%.2f'),'.tif');
                fullfilename = fullfile(save_directory,filename_tif); %Object fullfilename
                saveas(gcf,fullfilename); %Saving Matlab figure
                %Export to .obj format
                filename_obj = strcat('3D_Surface_Pendant_drop_Needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_scale(pxs_cm)_',num2str(scale_cm),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_n_',num2str(n),'_Rsquare_',num2str(gof.rsquare,'%.2f'),'.obj');
                fullfilename = fullfile(save_directory,filename_obj); %Object fullfilename
                saveobjmesh(fullfilename,X,Y,Z); %Exporting 3D drop surface as .obj
            %--------------------------------------------------------------
            %------------ SIMULATION OF DIFFERENT EDGE WIDTHS -------------
            %--------------------------------------------------------------
            elseif option == 3 %SIMULATION OF BLUR IN LIQUID DROP AND SOLID SUBSTRATE PROFILES (EDGE WIDTH)
                %- Simulation of drops with various edge widths
                for q = 1:length(W_pxs_array)
                    %Creation of needle and pendandt drop edge and normal profiles
                    res_v = resolution(2);
                    res_h = resolution(1);
                    xc = round(res_h/2); %x coordinate of the drop center in pixels
                    x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %x coordinates of the full drop profile in the final image
                    if needle %Construction of pendant drop with needle
                        zc = round(res_v/2);
                        z_drop = -z_ly_sim_pxs+z_ly_sim_pxs(end)+zc-(z_ly_sim_pxs(end)/2);  %z coordinates of the left drop profile in the final image with needle
                        % Creation of needle drop edge and normal profiles
                        needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
                        pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                        pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                        z_needle_max = zc-(z_ly_sim_pxs(end)/2); %maximum z coordinate of the needle
                        z_needle = linspace(0,z_needle_max,ceil(z_needle_max)); %z coordinates of the needle profile
                        x_needle_l = ones(1,length(z_needle))*pos_x_needle_l; %x coordinates of the left needle profile
                        x_needle_r = ones(1,length(z_needle))*pos_x_needle_r; %x coordinates of the right needle profile
                        [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge(x_needle_l,z_needle,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the left needle profile
                        [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge(x_needle_r,z_needle,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the right needle profile
                        g_needle_r = wrev(g_needle_r); %Fix grayscale order in the right profile of the needle
                        xn_needle = [xn_needle_l xn_needle_r]; %x coordinates of edge and normal points of the needle profile
                        zn_needle = [zn_needle_l zn_needle_r]; %z coordinates of edge and normal points of the needle profile
                        g_needle = [g_needle_l g_needle_r]; %grayscale value of edge and normal points of the needle profile
                    else %Construction of pendant drop without needle
                        z_drop = -z_ly_sim_pxs + z_ly_sim_pxs(end); %z coordinates of the left drop profile in the final image without needle
                    end
                    z_drop_lap_full = [wrev(z_drop);z_drop]; %z coordinates of the full drop profile in the final image)
                    %[xn,zn,g] = grayscale_edge(x_ly_sim_cm,z_ly_sim_cm,n_normal,g1,g2,W_lx);
                    [xn_drop,zn_drop,g_drop] = grayscale_edge(x_drop_lap_full,z_drop_lap_full,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the drop profile
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
                    %%{
                    %Save image
                    imagename = strcat('Edge_Width_Grayscale_Pendant_drop_Needle_',num2str(needle),'_c(cm-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_scale(pxs_cm)_',num2str(scale_cm),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_W(pxs)_',num2str(W_pxs_array(q),'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
                    fullimagename = strcat(imagename,imFmt); %Image name + image format
                    fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                    imwrite(I,fullfilename) %Save Image
                end
            %--------------------------------------------------------------
            %------- SIMULATION OF CAMERA VERTICAL MISALIGNMENT -----------
            %--------------------------------------------------------------
            elseif option == 4 % SIMULATION OF CAMERA VERTICAL MISALIGNMENT
                %- Simulation of camera vertical misalignment
                if needle == 0 % Camera verical misalignment can only be simulated with needle in the image
                    fprintf('Camera vertical misalignment cannot be simulated. Please, confirm needleON option. \n')
                    break;
                else
                    for q = 1:length(phi)
                        % Creation of needle and pendant drop edge and normal profiles
                        res_v = resolution(2);
                        res_h = resolution(1);
                        xc = round(res_h/2); %x coordinate of the drop center in pixels
                        x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %x coordinates of the full drop profile in the final image
                        zc = round(res_v/2);
                        z_drop = -z_ly_sim_pxs+z_ly_sim_pxs(end)+zc-(z_ly_sim_pxs(end)/2);  %z coordinates of the left drop profile in the final image with needle
            
                        % Creation of needle drop edge and normal profiles
                        needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
                        pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                        pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                        z_needle_max = zc-(z_ly_sim_pxs(end)/2); %maximum z coordinate of the needle
                        z_needle = linspace(-100,z_needle_max,ceil(z_needle_max)+100); %z coordinates of the needle profile
                        x_needle_l = ones(1,length(z_needle))*pos_x_needle_l; %x coordinates of the left needle profile
                        x_needle_r = ones(1,length(z_needle))*pos_x_needle_r; %x coordinates of the right needle profile
                        [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge(x_needle_l,z_needle,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left needle profile
                        [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge(x_needle_r,z_needle,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right needle profile
                        g_needle_r = wrev(g_needle_r); %Fix grayscale order in the right profile of the needle
                        xn_needle = [xn_needle_l xn_needle_r]; %x coordinates of edge and normal points of the needle profile
                        zn_needle = [zn_needle_l zn_needle_r]; %z coordinates of edge and normal points of the needle profile
                        g_needle = [g_needle_l g_needle_r]; %grayscale value of edge and normal points of the needle profile
                        z_drop_lap_full = [wrev(z_drop);z_drop]; %z coordinates of the full drop profile in the final image)
                        [xn_drop,zn_drop,g_drop] = grayscale_edge(x_drop_lap_full,z_drop_lap_full,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the drop profile
                        %- Combine needle and drop edge and normal points / Rotation of data
                        xn = [xn_needle xn_drop];
                        zn = [zn_needle zn_drop];
                        g = [g_needle g_drop];
                        [xn_rot,zn_rot] = rotate_data(xn,zn,xc_rot,zc_rot,phi(q)); %all data
        
                        % Generation of synthetical illuminated (grayscale) pendant drop images
                        M = zeros(res_v,res_h); %Array with image resolution size
                        I = uint8(M); %Image creation
                        %- Plotting needle and drop profile edge
                        for i=1:length(g)
                            if (xn_rot(i)>=0) && (xn_rot(i)<res_h) && (zn_rot(i)>=0) && (zn_rot(i)<res_v)
                                if double(I(floor(zn_rot(i))+1,floor(xn_rot(i))+1)) == 0
                                    I(floor(zn_rot(i))+1,floor(xn_rot(i))+1) = round(g(i));
                                else                
                                    I(floor(zn_rot(i))+1,floor(xn_rot(i))+1) = round((double(I(floor(zn_rot(i))+1,floor(xn_rot(i))+1)) + g(i))/2);
                                end
                            end
                        end
                        %- Filling background
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
                        %- Save image
                        imagename = strcat('Vertical_Misalignment_Grayscale_Pendant_drop_Needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_phi(°)_',num2str(phi(q),'%.3f')); %Image filename
                        fullimagename = strcat(imagename,imFmt); %Image name + image format
                        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                        imwrite(I,fullfilename) %Save Image
                    end
                end
            %--------------------------------------------------------------
            %----- SIMULATION OF RANDOM PERTURBATIONS IN DROP PROFILE -----
            %--------------------------------------------------------------
            elseif option == 5 %Simulation of random perturbations in drop profile
                %- Simulation of random perturbation in drop profile
                for p = 1:length(pert)
                    for q = 1:n_images
                        %Creation of needle and pendandt drop edge and normal profiles
                        res_v = resolution(2);
                        res_h = resolution(1);
                        xc = round(res_h/2); %x coordinate of the drop center in pixels
                        x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %x coordinates of the full drop profile in the final image
                        if needle %Construction of pendant drop with needle
                            zc = round(res_v/2);
                            z_drop = -z_ly_sim_pxs+z_ly_sim_pxs(end)+zc-(z_ly_sim_pxs(end)/2);  %z coordinates of the left drop profile in the final image with needle
                            % Creation of needle drop edge and normal profiles
                            needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
                            pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                            pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                            z_needle_max = zc-(z_ly_sim_pxs(end)/2); %maximum z coordinate of the needle
                            z_needle = linspace(0,z_needle_max,ceil(z_needle_max)+1); %z coordinates of the needle profile
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

                        %Random perturbation of the drop profile
                        [xn_drop,zn_drop,xn_drop_edge,zn_drop_edge] = normal_perturbation_edge(xn_drop,zn_drop,pert(p),n_normal);

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
                        %- Filling needle and drop profile
                        for i = 1:length(z_needle)
                            if round(z_needle(i)) ~= 0
                                for m = round(x_needle_r(i)):-1:1
                                    if I(round(z_needle(i)),m) == 0
                                        I(round(z_needle(i)),m) = g1;
                                    elseif I(round(z_needle(i)),m) == g2
                                        break;
                                    end
                                end
                            end
                        end
                        %Find drop apex
                        indx_apex = find(abs(zn_drop_edge) == max(abs(zn_drop_edge)),1,'first');
                        %Divide drop profile and determine the right side
                        xn_drop_edge_r = xn_drop_edge(indx_apex:end);
                        zn_drop_edge_r = zn_drop_edge(indx_apex:end);
                        %Filling image inside drop
                        for i = 1:length(zn_drop_edge_r)
                            for m = round(xn_drop_edge_r(i)):-1:1
                                if I(round(zn_drop_edge_r(i)),m) == 0
                                    I(round(zn_drop_edge_r(i)),m) = g1;
                                elseif I(round(zn_drop_edge_r(i)),m) == g2
                                    break;
                                end
                            end
                        end
                        %- Filling background
                        for i = 1:res_v %Search in all the image
                            for m = 1:res_h
                                if I(i,m) == 0
                                    I(i,m) = g2;
                                end
                            end
                        end
                        %%{
                        %Save image
                        imagename = strcat('Random_Edge_Perturbation_Pendant_drop_Needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_scale(pxs_cm)_',num2str(scale_cm),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_pert(pxs)_',num2str(pert(p),'%.2f'),'_image_',num2str(q)); %Image filename
                        fullimagename = strcat(imagename,imFmt); %Image name + image format
                        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                        imwrite(I,fullfilename) %Save Image
                    end
                end
            %--------------------------------------------------------------
            %--------------- SIMULATION OF LACK OF CONTRAST ---------------
            %--------------------------------------------------------------
            elseif option == 6 %Simulation of lack of contrast
                %- Simulation of lack of contrast
                for n = 1:length(diff_g)
                    g2 = g1 + diff_g(n); %Set g2 value
                    %Creation of needle and pendant drop edge and normal profiles
                    res_v = resolution(2);
                    res_h = resolution(1);
                    xc = round(res_h/2); %x coordinate of the drop center in pixels
                    x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %x coordinates of the full drop profile in the final image
                    if needle %Construction of pendant drop with needle
                        zc = round(res_v/2);
                        z_drop = -z_ly_sim_pxs+z_ly_sim_pxs(end)+zc-(z_ly_sim_pxs(end)/2);  %z coordinates of the left drop profile in the final image with needle
                        % Creation of needle drop edge and normal profiles
                        needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
                        pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                        pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                        z_needle_max = zc-(z_ly_sim_pxs(end)/2); %maximum z coordinate of the needle
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
                    %- Save image
                    imagename = strcat('Grayscale_Pendant_drop_Needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_scale(pxs_cm)_',num2str(scale_cm),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
                    fullimagename = strcat(imagename,imFmt); %Image name + image format
                    fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                    imwrite(I,fullfilename) %Save Image
                end
            %--------------------------------------------------------------
            %--------- SIMULATION OF NON-UNIFORM ILLUMINATION -------------
            %--------------------------------------------------------------
            elseif option == 7 %Simulation of non-uniform illumination
                %- Simulation of non-uniform illumination
                iref = res_v/2; %x-coordinate of non-uniform illumination center
                jref = res_h/2; %y-coordinate of non-uniform illumination center
                for n = 1:length(deltag)
                    %Creation of needle and pendant drop edge and normal profiles
                    res_v = resolution(2);
                    res_h = resolution(1);
                    xc = round(res_h/2); %x coordinate of the drop center in pixels
                    x_drop_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %x coordinates of the full drop profile in the final image
                    if needle %Construction of pendant drop with needle
                        zc = round(res_v/2);
                        z_drop = -z_ly_sim_pxs+z_ly_sim_pxs(end)+zc-(z_ly_sim_pxs(end)/2);  %z coordinates of the left drop profile in the final image with needle
                        % Creation of needle drop edge and normal profiles
                        needle_diameter_pxs = 2*r_h_cm*scale_cm; %Needle outer diameter in pixels
                        pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                        pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                        z_needle_max = zc-(z_ly_sim_pxs(end)/2); %maximum z coordinate of the needle
                        z_needle = linspace(0,z_needle_max,ceil(z_needle_max)); %z coordinates of the needle profile
                        x_needle_l = ones(1,length(z_needle))*pos_x_needle_l; %x coordinates of the left needle profile
                        x_needle_r = ones(1,length(z_needle))*pos_x_needle_r; %x coordinates of the right needle profile
                        [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge_nonuniform_illumination(x_needle_l,z_needle,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the left needle profile
                        [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge_nonuniform_illumination(x_needle_r,z_needle,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the right needle profile
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
                    else %Construction of pendant drop without needle
                        z_drop = -z_ly_sim_pxs + z_ly_sim_pxs(end); %z coordinates of the left drop profile in the final image without needle
                    end
                    z_drop_lap_full = [wrev(z_drop);z_drop]; %z coordinates of the full drop profile in the final image)
                    [xn_drop,zn_drop,g_drop] = grayscale_edge_nonuniform_illumination(x_drop_lap_full,z_drop_lap_full,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the drop profile
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
                    %- Filling background
                    for i = 1:res_v %Search from left to right
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
                    for i = 1:res_v %Search from right to left
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
                    %- Filling drop and needle
                    for i = 1:res_v %Search in all the image
                        for j = 1:res_h
                            if I(i,j) == 0
                                r = sqrt((i-iref)^2+(j-jref)^2);
                                I(i,j) = round(g1 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                %I(i,j) = g1;
                            end
                        end
                    end
                    %Save image
                    imagename = strcat('Non-uniform_Grayscale_Pendant_drop_Needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(k),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_V(uL)_',num2str(V_uL(k),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_deltag_',num2str(round(deltag(n)))); %Image filename
                    fullimagename = strcat(imagename,imFmt); %Image name + image format
                    fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                    imwrite(I,fullfilename) %Save Image
                    %imwrite(I,'pendantdrop_ilum.tif') %Save Image
                end
            end
        end
        %Write excel file with important pendant drop profile properties
        fprintf('Exporting excel file containing important pendant drop properties... \n');
        varnames = {'b[cm-1]','r0[cm]','De[cm]','h_De[cm]','h[cm]','Ds[cm]','A[cm²]','V[uL]'}; %Variable names in Excel file
        T = table(b',r0',d_eq_cm,h_eq_cm,h_cm,d_s_cm,A_cm2,V_uL,'VariableNames',varnames);
        sheet_name = strcat('c = ',num2str(c,'%.2f'),' cm-2',' r_h = ',num2str(r_h_cm,'%.4f'),' cm'); %Sheet name in Excel
        savename_xls = strcat('Grayscale_Pendant_drop_(Needle_',num2str(needle),'_c_cm-2_',num2str(c),'_r_h_cm_',num2str(r_h_cm,'%.4f'),').xlsx'); %Excel filename
        fullsavename_xls = fullfile(save_directory,savename_xls);
        writetable(T,fullsavename_xls,'Sheet',sheet_name) %Write table to an excel file

        fprintf('Simulation completed! \n');
        pause(2)
    end
    clc
    clear all
    close all
end