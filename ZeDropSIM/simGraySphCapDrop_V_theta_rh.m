function simGraySphCapDrop_V_theta_rh()
%SIMGRAYSPHCAPDROP_V_THETA_RH Simulation of grayscale images of spherical 
% cap drops. Some source of errros can be simulated and the 3D drop surface
% can also be exported.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of grayscale images of spherical cap drops. Some source of 
%   errors, such as edge width, camera vertical misalignment, random 
%   perturbation in the drop profile, lack of contrast, non-uniform 
%   illumination, random noise, Gaussian noise, impulse salt-pepper noise
%   and satellite droplets, can be simulated. The 3D drop surface can also 
%   be exported. Depending on the simulated option, the program requires 
%   the user to enter some specific parameters. As common input parameters,
%   any option requires the user to enter: the minimum drop volume [uL]
%   (Vmin), the maximum drop volume [uL](Vmax), the number of drop volume 
%   levels that is going to be simulated between Vmin and Vmax (nVol), the 
%   minimum contact angle [°](thetamin), the maximum contact angle [°]
%   (thetamax), the number of contact angle levels that is going to be 
%   simulated between thetamin and thetamax (ntheta), the needle diameter 
%   [mm](needleDiam), the possibility of simulating images with and without
%   needle (needleON), the image resolution [pixels](res_h and res_v), the 
%   image scale [pixels/mm](scale), edge width (W_pxs), maximum and minimum
%   image intensity levels (g1 and g2). The height of the solid substrate 
%   in the image is fixed. The routine also allows the user to choose the 
%   directory where the images are going to be exported. An excel file 
%   (.xlsx) is created containing important properties of the spherical cap
%   drops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while 1
    fprintf('------------------ SPHERICAL CAP DROP - SIMULATION OF GRAYSCALE IMAGES --------------------\n');
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
        nVol = input("- Number of drop volume levels simulated [3]: "); %Number of simulated images between minimum and maximum drop volume
        if isempty(nVol)
            nVol = 3;
        end
        V_uL = linspace(Vmin,Vmax,nVol);
        thetamin = input("- Minimum contact angle in degree [30]: "); %Minimum contact angle in degree
        if isempty(thetamin)
            thetamin = 30;
        end
        while 1
            thetamax = input("- Maximum contact angle in degree [160]: "); %Maximum contact angle in degree
            if isempty(thetamax)
                thetamax = 160;
            end
            if thetamax < thetamin
                fprintf('Invalid maximum contact angle. Please, set a lower value. \n');
            else
                break;
            end
        end
        nTheta = input("- Number of contact angle levels simulated [4]: "); %Number of contact angle levels simulated
        if isempty(nTheta)
            nTheta = 4;
        end
        theta = linspace(thetamin,thetamax,nTheta);
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
        scale_cm = scale*10; %Image scale in pixels/cm
        imFmt = input("- Image format [.tif]: ", "s");
        if isempty(imFmt)
            imFmt = ".tif";
        end
        xc = round(res_h/2); %x coordinate of the drop center in pixels
        Ncircle = 10000; %Number of points of circle edge
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
        fprintf('Selection of the directory to export files ... \n');
        save_directory = uigetdir('C:\'); %Save directory

        %- VARIABLES INITIALIZATION
        r_pxs = zeros(length(theta),length(V_uL)); %Radius of the spherical drop in pixels
        a_pxs = zeros(length(theta),length(V_uL)); %Wetted radius in pixels
        b_pxs = zeros(length(theta),length(V_uL)); %Distance between the center of the spherical drop and the cut level in pixels
        h_pxs = zeros(length(theta),length(V_uL)); %Drop height in pixels
        r_cm = zeros(length(theta),length(V_uL)); %Radius of the spherical drop in cm
        a_cm = zeros(length(theta),length(V_uL)); %Wetted radius in pixels in cm
        b_cm = zeros(length(theta),length(V_uL)); %Distance between the center of the spherical drop and the cut level in cm
        h_cm = zeros(length(theta),length(V_uL)); %Drop height in cm
        xCPleft_pxs = zeros(length(theta),length(V_uL)); %x coordinate of the left triple contact point
        yCPleft_pxs = zeros(length(theta),length(V_uL)); %y coordinate of the left triple contact point
        xCPright_pxs = zeros(length(theta),length(V_uL)); %x coordinate of the right triple contact point
        yCPright_pxs = zeros(length(theta),length(V_uL)); %y coordinate of the right triple contact point

        %------------------------------------------------------------------------
        %-------------- SIMULATION OF SPHERICAL CAP DROP IMAGES -----------------
        %------------------------------------------------------------------------
        fprintf("--------------------------------------------------------------------------------\n")
        fprintf("Simulation of grayscale spherical cap drop images... \n")
        for i = 1:length(V_uL)
            for j = 1:length(theta)
                fprintf('Running V = %.2f uL / CA = %.2f °... \n',V_uL(i),theta(j));
                %Drop dimensions
                [r_mm,a_mm,b_mm,h_mm] = spherical_cap_dim(V_uL(i),theta(j));
                r_cm(j,i) = r_mm/10;
                a_cm(j,i) = a_mm/10;
                b_cm(j,i) = b_mm/10;
                h_cm(j,i) = h_mm/10;

                %Conversion to pixels
                r_pxs(j,i) = r_cm(j,i)*scale_cm; %em pixels
                a_pxs(j,i) = a_cm(j,i)*scale_cm; %em pixels
                b_pxs(j,i) = b_cm(j,i)*scale_cm; %em pixels
                h_pxs(j,i) = h_cm(j,i)*scale_cm; %em pixels

                xCPleft_pxs(j,i) = xc - a_pxs(j,i);
                xCPright_pxs(j,i) = xc + a_pxs(j,i);

                %--------------------------------------------------------------
                %-------------- GENERATION OF GRAYSCALE IMAGES ----------------
                %--------------------------------------------------------------
                if option == 1 || option == 8 || option == 9 || option == 10 || option == 11 % Generation of grayscale images
                    %Creation of edge and normal profiles
                    if theta(j) <= 90 %theta<=90°
                        zc = res_v - hsubstrato + b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        yCPleft_pxs(j,i) = zc - b_pxs(j,i);
                        yCPright_pxs(j,i) = zc - b_pxs(j,i);
                    else %theta>90°
                        zc = res_v - hsubstrato - b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        yCPleft_pxs(j,i) = zc + b_pxs(j,i);
                        yCPright_pxs(j,i) = zc + b_pxs(j,i);
                    end
                    %zc = res_v-hsubstrato-b_pxs(j,i); %Coordenada y do centro da gota em pixels
                    [x,z] = circle(xc,zc,r_pxs(j,i),theta(j),Ncircle); %Spherical cap coordinates in pixels
                    %Flip drop around x-axis that pass through zc
                    x_circ = x;
                    z_circ = 2*zc-z;

                    %Creation of spherical cap drop and needle edge normal profile
                    if needle_in
                        %Construction of the needle profile
                        needle_diameter_pxs = (needleDiam/10)*scale_cm;
                        pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                        pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                        if theta(j)>90
                            indx_l = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                            indx_l = round(length(x_circ)/4 + indx_l);
                            indx_r = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                            indx_r = round(length(x_circ)/4 + indx_r);
                        else
                            indx_l = find((abs(x_circ-pos_x_needle_l)) == min(abs(x_circ-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                            indx_r = find((abs(x_circ-pos_x_needle_r)) == min(abs(x_circ-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                        end
                        z_needle_max_l = z_circ(indx_l);
                        z_needle_max_r = z_circ(indx_r);
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
                        x_drop_l = x_circ(indx_l:end);
                        z_drop_l = z_circ(indx_l:end);
                        x_drop_r = x_circ(1:indx_r);
                        z_drop_r = z_circ(1:indx_r);
                        [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the drop
                        %g_drop_l = wrev(g_drop_l);
                        [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the drop
                        %g_drop_r = wrev(g_drop_r);
                        xn_drop = [xn_drop_l xn_drop_r];
                        zn_drop = [zn_drop_l zn_drop_r];
                        g_drop = [g_drop_l g_drop_r];
                    else
                        [xn_drop,zn_drop,g_drop] = grayscale_edge(x_circ,z_circ,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the drop profile
                        %g_drop = wrev(g_drop);
                    end
                    %Creation of solid subratre edge normal profile
                    x_sol_l = linspace(0,x_circ(end),ceil(x_circ(end)));
                    z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                    x_sol_r = linspace(x_circ(1),res_h,ceil(res_h-x_circ(1)));
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
                    for k=1:length(g)
                        if (xn(k)>=0) && (xn(k)<res_h) && (zn(k)>=0) && (zn(k)<res_v)
                            if double(I(floor(zn(k))+1,floor(xn(k))+1)) == 0
                                I(floor(zn(k))+1,floor(xn(k))+1) = round(g(k));
                            else
                                I(floor(zn(k))+1,floor(xn(k))+1) = round((double(I(floor(zn(k))+1,floor(xn(k))+1)) + g(k))/2);
                            end
                        end
                    end
                    %Iorig = I;
                    %- Filling background
                    for k = 1:(res_v-hsubstrato) %Search from left to right
                        for m = 1:res_h
                            if I(k,m) == 0
                                I(k,m) = g2;
                            else
                                break;
                            end
                        end
                    end
                    for k = 1:(res_v-hsubstrato) %Search from right to left
                        for m = res_h:-1:1
                            if I(k,m) == 0
                                I(k,m) = g2;
                            else
                                break;
                            end
                        end
                    end
                    %- Filling solid substrate,drop and needle
                    for k = 1:res_v %Search in all the image
                        for m = 1:res_h
                            if I(k,m) == 0
                                I(k,m) = g1;
                            end
                        end
                    end
                    if option == 1
                        %Save images
                        imagename = strcat('Grayscale_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
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
                                imagename = strcat('Noise_Grayscale_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_noise_',num2str(round(noise(q))),'_n_images_',num2str(p)); %Image filename
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
                                imagename = strcat('Gaussian_Noise_Grayscale_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_g_noise_mean_',num2str(g_noise_mean),'_g_noise_var_',num2str(g_noise_var(q)),'_n_images_',num2str(p)); %Image filename
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
                                imagename = strcat('Salt_Pepper_Noise_Grayscale_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_salt_pepper_density_',num2str(sp_noise_density(q)),'_n_images_',num2str(p)); %Image filename
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
                        for p = 1:length(a)
                            N = round((4*a(p)*res_h*res_v)/(pi*(d_sat_pxs^2))); %Number of satellite droplets
                            fprintf('Running a = %.3f and N = %d ... \n',a(p),N);
                            for k = 1:n_images
                                Isat = I;
                                for q = 1:N
                                    [columnsInImage, rowsInImage] = meshgrid(1:res_h, 1:res_v);
                                    x_sat_center = randi([1 res_h],1);
                                    z_sat_center = randi([1 res_v],1);
                                    circlePixels = (rowsInImage-z_sat_center).^2 + (columnsInImage-x_sat_center).^2 <= (d_sat_pxs/2).^2;
                                    indx_circle = find(circlePixels);
                                    Isat(indx_circle) = g1; %Plotiing satellite droplets on the image
                                end
                                %Save satellite droplet images
                                imagename = strcat('Satellite_Grayscale_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_d_sat(cm)_',num2str(d_sat_cm,'%.3f'),'_a_',num2str(a(p),'%.4f'),'_n_images_',num2str(k)); %Image filename
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
                    [x_surf,z_surf] = circle(0,-r_pxs(j,i),r_pxs(j,i),theta(j),Ncircle); %Creating the spherical drop with apex centered at (0,0)
                    %ind_x = find(abs(x_surf)==min(abs(x_surf)),1,'first');
                    ind_x = round(length(x_surf)/2); %find the position equivalent to half-of the circle
                    x_surf_lap = [x_surf(1:ind_x) 0];
                    z_surf_lap = [z_surf(1:ind_x) 0];
                    x_surf_lap = wrev(x_surf_lap);
                    z_surf_lap = wrev(z_surf_lap);
                    if needle_in
                        %Construction of the needle profile
                        needle_diameter_pxs = (needleDiam/10)*scale_cm;
                        pos_x_needle = round(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                        %if theta(j)>90
                        %indx_needle = find((abs(x_surf_lap(round(length(x_surf_lap)/2:end))-pos_x_needle)) == min(abs(x_surf_lap(round(length(x_surf_lap)/2:end))-pos_x_needle)),1,'first'); %Find the index of the instersection between the drop profile and the needle
                        %indx_needle = length(x_surf_lap)/2 + indx_needle;
                        %else
                        indx_needle = find((abs(x_surf_lap-pos_x_needle)) == min(abs(x_surf_lap-pos_x_needle)),1,'first'); %Find the index of the instersection between the drop profile and the needle
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
                    %figure
                    mesh(X,Y,Z)
                    axis equal
                    %Saving MATLAB figure
                    filename_tif = strcat('3D_Surface_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_Rsquare_',num2str(gof.rsquare,'%.2f'),'.tif');
                    fullfilename = fullfile(save_directory,filename_tif); %Object fullfilename
                    saveas(gcf,fullfilename); %Saving Matlab figure
                    %Export to .obj format
                    filename_obj = strcat('3D_Surface_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_Rsquare_',num2str(gof.rsquare,'%.2f'),'.obj');
                    fullfilename = fullfile(save_directory,filename_obj); %Object fullfilename
                    saveobjmesh(fullfilename,X,Y,Z); %Exporting 3D drop surface as .obj
                %--------------------------------------------------------------
                %------------ SIMULATION OF DIFFERENT EDGE WIDTHS -------------
                %--------------------------------------------------------------
                elseif option == 3 %SIMULATION OF BLUR IN LIQUID DROP AND SOLID SUBSTRATE PROFILES (EDGE WIDTH)
                    for q = 1:length(W_pxs_array)
                        %Creation of edge and normal profiles
                        if theta(j) <= 90 %theta<=90°
                            zc = res_v - hsubstrato + b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        else %theta>90°
                            zc = res_v - hsubstrato - b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        end
                        %zc = res_v-hsubstrato-b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        [x,z] = circle(xc,zc,r_pxs(j,i),theta(j),Ncircle); %Spherical cap coordinates in pixels
                        %Flip drop around x-axis that pass through zc
                        x_circ = x;
                        z_circ = 2*zc-z;

                        %Creation of spherical cap drop and needle edge normal profile
                        if needle_in
                            %Construction of the needle profile
                            needle_diameter_pxs = (needleDiam/10)*scale_cm;
                            pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                            pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                            if theta(j)>90
                                indx_l = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                                indx_l = length(x_circ)/4 + indx_l;
                                indx_r = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                                indx_r = length(x_circ)/4 + indx_r;
                            else
                                indx_l = find((abs(x_circ-pos_x_needle_l)) == min(abs(x_circ-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                                indx_r = find((abs(x_circ-pos_x_needle_r)) == min(abs(x_circ-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                            end
                            z_needle_max_l = z_circ(indx_l);
                            z_needle_max_r = z_circ(indx_r);
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
                            x_drop_l = x_circ(indx_l:end);
                            z_drop_l = z_circ(indx_l:end);
                            x_drop_r = x_circ(1:indx_r);
                            z_drop_r = z_circ(1:indx_r);
                            [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the left portion of the drop
                            [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the right portion of the drop
                            xn_drop = [xn_drop_l xn_drop_r];
                            zn_drop = [zn_drop_l zn_drop_r];
                            g_drop = [g_drop_l g_drop_r];
                        else
                            [xn_drop,zn_drop,g_drop] = grayscale_edge(x_circ,z_circ,n_normal,g1,g2,W_pxs_array(q)); %Calculation of the normal points to the drop profile
                        end
                        %Creation of solid subratre edge normal profile
                        x_sol_l = linspace(0,x_circ(end),ceil(x_circ(end)));
                        z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                        x_sol_r = linspace(x_circ(1),res_h,ceil(res_h-x_circ(1)));
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
                        for k=1:length(g)
                            if (xn(k)>=0) && (xn(k)<res_h) && (zn(k)>=0) && (zn(k)<res_v)
                                if double(I(floor(zn(k))+1,floor(xn(k))+1)) == 0
                                    I(floor(zn(k))+1,floor(xn(k))+1) = round(g(k));
                                else
                                    I(floor(zn(k))+1,floor(xn(k))+1) = round((double(I(floor(zn(k))+1,floor(xn(k))+1)) + g(k))/2);
                                end
                            end
                        end
                        %%{
                        if W_pxs_array(q) == 0
                            %- Filling background
                            for k = 1:(res_v-hsubstrato) %Search from left to right
                                for m = 1:res_h
                                    if I(k,m) == 0
                                        I(k,m) = g2;
                                    else
                                        break;
                                    end
                                end
                            end
                            for k = 1:(res_v-hsubstrato) %Search from right to left
                                for m = res_h:-1:1
                                    if I(k,m) == 0
                                        I(k,m) = g2;
                                    else
                                        break;
                                    end
                                end
                            end
                            %- Filling solid substrate,drop and needle
                            for k = 1:res_v %Search in all the image
                                for m = 1:res_h
                                    if I(k,m) == 0
                                        I(k,m) = g1;
                                    end
                                end
                            end
                        else
                            if needle_in == 1
                                %Filling needle
                                for k = 1:length(z_needle_r)
                                    if round(z_needle_r(k)) ~= 0
                                        for m = round(x_needle_r(k)):-1:1
                                            if I(round(z_needle_r(k)),m) == 0
                                                I(round(z_needle_r(k)),m) = g1;
                                            elseif I(round(z_needle_r(k)),m) == g2
                                                break;
                                            end
                                        end
                                    end
                                end
                                %Filling image inside drop
                                for k = 1:length(z_drop_r)
                                    for m = round(x_drop_r(k)):-1:1
                                        if I(round(z_drop_r(k)),m) == 0
                                            I(round(z_drop_r(k)),m) = g1;
                                        elseif I(round(z_drop_r(k)),m) == g2
                                            break;
                                        end
                                    end
                                end
                            else
                                %Find drop apex
                                indx_apex = find(abs(z_circ) == min(abs(z_circ)),1,'first');
                                %Divide drop profile and determine the right side
                                x_circ_r = x_circ(1:indx_apex);
                                z_circ_r = z_circ(1:indx_apex);
                                %Filling image inside drop
                                for k = 1:length(z_circ_r)
                                    for m = round(x_circ_r(k)):-1:1
                                        if I(round(z_circ_r(k)),m) == 0
                                            I(round(z_circ_r(k)),m) = g1;
                                        elseif I(round(z_circ_r(k)),m) == g2
                                            break;
                                        end
                                    end
                                end
                            end
                            %Filling solid substrate
                            for m = 1:res_h %Search from bottom to top - left to right
                                for k = res_v:-1:1
                                    if I(k,m) == g2
                                        break;
                                    elseif I(k,m) == 0
                                        I(k,m) = g1;
                                    end
                                end
                            end
                            %- Filling background
                            for k = 1:res_v %Search in all the image
                                for m = 1:res_h
                                    if I(k,m) == 0
                                        I(k,m) = g2;
                                    end
                                end
                            end
                        end
                        %%{
                        %Save images
                        imagename = strcat('Edge_Width_Grayscale_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs_array(q),'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
                        fullimagename = strcat(imagename,imFmt); %Image name + image format
                        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                        imwrite(I,fullfilename) %Save Image
                    end
                %--------------------------------------------------------------
                %------- SIMULATION OF CAMERA VERTICAL MISALIGNMENT -----------
                %--------------------------------------------------------------
                elseif option == 4 % SIMULATION OF CAMERA VERTICAL MISALIGNMENT
                    for q = 1:length(phi)
                        %Creation of edge and normal profiles
                        if theta(j) <= 90 %theta<=90°
                            zc = res_v - hsubstrato + b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        else %theta>90°
                            zc = res_v - hsubstrato - b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        end
                        %zc = res_v-hsubstrato-b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        [x,z] = circle(xc,zc,r_pxs(j,i),theta(j),Ncircle); %Spherical cap coordinates in pixels
                        %Flip drop around x-axis that pass through zc
                        x_circ = x;
                        z_circ = 2*zc-z;

                        %Creation of spherical cap drop and needle edge normal profile
                        if needle_in
                            %Construction of the needle profile
                            needle_diameter_pxs = (needleDiam/10)*scale_cm;
                            pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                            pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                            if theta(j)>90
                                indx_l = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                                indx_l = length(x_circ)/4 + indx_l;
                                indx_r = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                                indx_r = length(x_circ)/4 + indx_r;
                            else
                                indx_l = find((abs(x_circ-pos_x_needle_l)) == min(abs(x_circ-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                                indx_r = find((abs(x_circ-pos_x_needle_r)) == min(abs(x_circ-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                            end
                            %indx_l = find((abs(x_circ-pos_x_needle_l)) == min(abs(x_circ-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                            %indx_r = find((abs(x_circ-pos_x_needle_r)) == min(abs(x_circ-pos_x_needle_r)),1,'last'); %Find the index of the instersection between the right side of the needle and the drop
                            z_needle_max_l = z_circ(indx_l);
                            z_needle_max_r = z_circ(indx_r);
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
                            x_drop_l = x_circ(indx_l:end);
                            z_drop_l = z_circ(indx_l:end);
                            x_drop_r = x_circ(1:indx_r);
                            z_drop_r = z_circ(1:indx_r);
                            [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the drop
                            [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the drop
                            xn_drop = [xn_drop_l xn_drop_r];
                            zn_drop = [zn_drop_l zn_drop_r];
                            g_drop = [g_drop_l g_drop_r];
                        else
                            [xn_drop,zn_drop,g_drop] = grayscale_edge(x_circ,z_circ,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the drop profile
                        end
                        %Creation of solid subratre edge normal profile
                        x_sol_l = linspace(-100,x_circ(end),ceil(x_circ(end)+100));
                        z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                        x_sol_r = linspace(x_circ(1),res_h+100,ceil(res_h-x_circ(1)+100));
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
                            [x_circ_rot,z_circ_rot] = rotate_data(x_circ,z_circ,xc_rot,zc_rot,phi(q));
                        end
                        [xn_rot,zn_rot] = rotate_data(xn,zn,xc_rot,zc_rot,phi(q)); %all data

                        %Generation of synthetical rotated grayscale spherical drop images
                        M = zeros(res_v,res_h); %Array with image resolution size
                        I = uint8(M); %Image creation
                        %- Plotting solid substrate, drop profile and needle edge
                        for k=1:length(g)
                            if (xn_rot(k)>=0) && (xn_rot(k)<res_h) && (zn_rot(k)>=0) && (zn_rot(k)<res_v)
                                if double(I(floor(zn_rot(k))+1,floor(xn_rot(k))+1)) == 0
                                    I(floor(zn_rot(k))+1,floor(xn_rot(k))+1) = round(g(k));
                                else
                                    I(floor(zn_rot(k))+1,floor(xn_rot(k))+1) = round((double(I(floor(zn_rot(k))+1,floor(xn_rot(k))+1)) + g(k))/2);
                                end
                            end
                        end
                        %Spherical drop with needle
                        if needle_in
                            if phi(q)<0 %Counterclockwise rotation
                                %Filling image inside needle
                                for k = 1:length(z_needle_l_rot)
                                    if round(z_needle_l_rot(k)) > 0
                                        for m = round(x_needle_l_rot(k)):round(res_h)
                                            if I(round(z_needle_l_rot(k)),m) == 0
                                                I(round(z_needle_l_rot(k)),m) = g1;
                                            elseif I(round(z_needle_l_rot(k)),m) == g2
                                                break;
                                            end
                                        end
                                    end
                                end
                                %Filling image inside drop
                                for k = 1:length(z_drop_l_rot)
                                    for m = round(x_drop_l_rot(k)):round(res_h)
                                        if I(round(z_drop_l_rot(k)),m) == 0
                                            I(round(z_drop_l_rot(k)),m) = g1;
                                        elseif I(round(z_drop_l_rot(k)),m) == g2
                                            break;
                                        end
                                    end
                                end
                            else %Clockwise rotation
                                %Filling image inside needle
                                for k = 1:length(z_needle_r_rot)
                                    if round(z_needle_r_rot(k)) > 0
                                        for m = round(x_needle_r_rot(k)):-1:1
                                            if I(round(z_needle_r_rot(k)),m) == 0
                                                I(round(z_needle_r_rot(k)),m) = g1;
                                            elseif I(round(z_needle_r_rot(k)),m) == g2
                                                break;
                                            end
                                        end
                                    end
                                end
                                %Filling image inside drop
                                for k = 1:length(z_drop_r_rot)
                                    for m = round(x_drop_r_rot(k)):-1:1
                                        if I(round(z_drop_r_rot(k)),m) == 0
                                            I(round(z_drop_r_rot(k)),m) = g1;
                                        elseif I(round(z_drop_r_rot(k)),m) == g2
                                            break;
                                        end
                                    end
                                end
                            end
                            %Spherical drop without needle
                        else
                            %Finding drop apex and dividing drop profile
                            x_circ_rot = wrev(x_circ_rot);
                            z_circ_rot = wrev(z_circ_rot);
                            indx_z = find(z_circ_rot == min(z_circ_rot),1,'first'); %Find drop apex
                            x_circ_l_rot = x_circ_rot(1:indx_z);
                            z_circ_l_rot = z_circ_rot(1:indx_z);
                            x_circ_r_rot = x_circ_rot(indx_z:end);
                            z_circ_r_rot = z_circ_rot(indx_z:end);
                            if phi(q)<0 %Counterclockwise rotation
                                %Filling image inside drop
                                for k = 1:length(z_circ_l_rot)
                                    for m = round(x_circ_l_rot(k)):round(res_h)
                                        if I(round(z_circ_l_rot(k)),m) == 0
                                            I(round(z_circ_l_rot(k)),m) = g1;
                                        elseif I(round(z_circ_l_rot(k)),m) == g2
                                            break;
                                        end
                                    end
                                end
                            else %Clockwise rotation
                                %Filling image inside drop
                                for k = 1:length(z_circ_r_rot)
                                    for m = round(x_circ_r_rot(k)):-1:1
                                        if I(round(z_circ_r_rot(k)),m) == 0
                                            I(round(z_circ_r_rot(k)),m) = g1;
                                        elseif I(round(z_circ_r_rot(k)),m) == g2
                                            break;
                                        end
                                    end
                                end
                            end
                        end
                        %Filling solid substrate
                        for m=1:res_h %Search from bottom to top
                            for k =res_v:-1:1
                                if I(k,m) == 0
                                    I(k,m) = g1;
                                elseif I(k,m) == g2
                                    break;
                                end
                            end
                        end
                        %Filling background
                        for k = 1:res_v %Search in all the image
                            for m = 1:res_h
                                if I(k,m) == 0
                                    I(k,m) = g2;
                                end
                            end
                        end
                        %Save images
                        imagename = strcat('Vertical_Misalignment_Grayscale_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_phi(°)_',num2str(phi(q))); %Image filename
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
                            %Creation of edge and normal profiles
                            if theta(j) <= 90 %theta<=90°
                                zc = res_v - hsubstrato + b_pxs(j,i); %Coordenada y do centro da gota em pixels
                            else %theta>90°
                                zc = res_v - hsubstrato - b_pxs(j,i); %Coordenada y do centro da gota em pixels
                            end
                            %zc = res_v-hsubstrato-b_pxs(j,i); %Coordenada y do centro da gota em pixels
                            [x,z] = circle(xc,zc,r_pxs(j,i),theta(j),Ncircle); %Spherical cap coordinates in pixels
                            %Flip drop around x-axis that pass through zc
                            x_circ = x;
                            z_circ = 2*zc-z;

                            %Creation of spherical cap drop and needle edge normal profile
                            if needle_in
                                %Construction of the needle profile
                                needle_diameter_pxs = (needleDiam/10)*scale_cm;
                                pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                                pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                                if theta(j)>90
                                    indx_l = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                                    indx_l = length(x_circ)/4 + indx_l;
                                    indx_r = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                                    indx_r = length(x_circ)/4 + indx_r;
                                else
                                    indx_l = find((abs(x_circ-pos_x_needle_l)) == min(abs(x_circ-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                                    indx_r = find((abs(x_circ-pos_x_needle_r)) == min(abs(x_circ-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                                end
                                z_needle_max_l = z_circ(indx_l);
                                z_needle_max_r = z_circ(indx_r);
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
                                x_drop_l = x_circ(indx_l:end);
                                z_drop_l = z_circ(indx_l:end);
                                x_drop_r = x_circ(1:indx_r);
                                z_drop_r = z_circ(1:indx_r);
                                [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the drop
                                [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the drop
                                xn_drop = [xn_drop_l xn_drop_r];
                                zn_drop = [zn_drop_l zn_drop_r];
                                g_drop = [g_drop_l g_drop_r];
                            else
                                [xn_drop,zn_drop,g_drop] = grayscale_edge(x_circ,z_circ,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the drop profile
                            end

                            %Random perturbation of the drop profile
                            [xn_drop,zn_drop,xn_drop_edge,zn_drop_edge] = normal_perturbation_edge(xn_drop,zn_drop,pert(p),n_normal);

                            %Creation of solid subratre edge normal profile
                            x_sol_l = linspace(0,x_circ(end),ceil(x_circ(end)));
                            z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                            x_sol_r = linspace(x_circ(1),res_h,ceil(res_h-x_circ(1)));
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
                            for k=1:length(g)
                                if (xn(k)>=0) && (xn(k)<res_h) && (zn(k)>=0) && (zn(k)<res_v)
                                    if double(I(floor(zn(k))+1,floor(xn(k))+1)) == 0
                                        I(floor(zn(k))+1,floor(xn(k))+1) = round(g(k));
                                    else
                                        I(floor(zn(k))+1,floor(xn(k))+1) = round((double(I(floor(zn(k))+1,floor(xn(k))+1)) + g(k))/2);
                                    end
                                end
                            end
                            if needle_in == 1
                                for k = 1:length(z_needle_r)
                                    if round(z_needle_r(k)) ~= 0
                                        for m = round(x_needle_r(k)):-1:1
                                            if I(round(z_needle_r(k)),m) == 0
                                                I(round(z_needle_r(k)),m) = g1;
                                            elseif I(round(z_needle_r(k)),m) == g2
                                                break;
                                            end
                                        end
                                    end
                                end
                                %Find drop apex
                                indx_apex_sup = find(abs(zn_drop_edge) == min(abs(zn_drop_edge)),1,'last');
                                indx_apex_inf = find(abs(zn_drop_edge) == max(abs(zn_drop_edge)),1,'last');
                                %Divide drop profile and determine the right side
                                xn_drop_edge_r = xn_drop_edge(indx_apex_inf:indx_apex_sup);
                                zn_drop_edge_r = zn_drop_edge(indx_apex_inf:indx_apex_sup);
                                %Filling image inside drop
                                for k = 1:length(zn_drop_edge_r)
                                    for m = round(xn_drop_edge_r(k)):-1:1
                                        if I(round(zn_drop_edge_r(k)),m) == 0
                                            I(round(zn_drop_edge_r(k)),m) = g1;
                                        elseif I(round(zn_drop_edge_r(k)),m) == g2
                                            break;
                                        end
                                    end
                                end
                            else
                                %Find drop apex
                                indx_apex = find(abs(zn_drop_edge) == min(abs(zn_drop_edge)),1,'first');
                                %Divide drop profile and determine the right side
                                xn_drop_edge_r = xn_drop_edge(1:indx_apex);
                                zn_drop_edge_r = zn_drop_edge(1:indx_apex);
                                %Filling image inside drop
                                for k = 1:length(zn_drop_edge_r)
                                    for m = round(xn_drop_edge_r(k)):-1:1
                                        if I(round(zn_drop_edge_r(k)),m) == 0
                                            I(round(zn_drop_edge_r(k)),m) = g1;
                                        elseif I(round(zn_drop_edge_r(k)),m) == g2
                                            break;
                                        end
                                    end
                                end
                            end
                            %Filling solid substrate
                            for m = 1:res_h %Search from bottom to top - left to right
                                for k = res_v:-1:1
                                    if I(k,m) == g2
                                        break;
                                    elseif I(k,m) == 0
                                        I(k,m) = g1;
                                    end
                                end
                            end
                            %- Filling background
                            for k = 1:res_v %Search in all the image
                                for m = 1:res_h
                                    if I(k,m) == 0
                                        I(k,m) = g2;
                                    end
                                end
                            end
                            %Save images
                            imagename = strcat('Random_Edge_Perturbation_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_pert(pxs)_',num2str(pert(p),'%.2f'),'_image_',num2str(q)); %Image filename
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
                        g2 = g1 + diff_g(n); %Setting g2 value
                        %Creation of edge and normal profiles
                        if theta(j) <= 90 %theta<=90°
                            zc = res_v - hsubstrato + b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        else %theta>90°
                            zc = res_v - hsubstrato - b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        end
                        %zc = res_v-hsubstrato-b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        [x,z] = circle(xc,zc,r_pxs(j,i),theta(j),Ncircle); %Spherical cap coordinates in pixels
                        %Flip drop around x-axis that pass through zc
                        x_circ = x;
                        z_circ = 2*zc-z;

                        %Creation of spherical cap drop and needle edge normal profile
                        if needle_in
                            %Construction of the needle profile
                            needle_diameter_pxs = (needleDiam/10)*scale_cm;
                            pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                            pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                            if theta(j)>90
                                indx_l = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                                indx_l = length(x_circ)/4 + indx_l;
                                indx_r = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                                indx_r = length(x_circ)/4 + indx_r;
                            else
                                indx_l = find((abs(x_circ-pos_x_needle_l)) == min(abs(x_circ-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                                indx_r = find((abs(x_circ-pos_x_needle_r)) == min(abs(x_circ-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                            end
                            z_needle_max_l = z_circ(indx_l);
                            z_needle_max_r = z_circ(indx_r);
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
                            x_drop_l = x_circ(indx_l:end);
                            z_drop_l = z_circ(indx_l:end);
                            x_drop_r = x_circ(1:indx_r);
                            z_drop_r = z_circ(1:indx_r);
                            [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the left portion of the drop
                            %g_drop_l = wrev(g_drop_l);
                            [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the right portion of the drop
                            %g_drop_r = wrev(g_drop_r);
                            xn_drop = [xn_drop_l xn_drop_r];
                            zn_drop = [zn_drop_l zn_drop_r];
                            g_drop = [g_drop_l g_drop_r];
                        else
                            [xn_drop,zn_drop,g_drop] = grayscale_edge(x_circ,z_circ,n_normal,g1,g2,W_pxs); %Calculation of the normal points to the drop profile
                            %g_drop = wrev(g_drop);
                        end
                        %Creation of solid subratre edge normal profile
                        x_sol_l = linspace(0,x_circ(end),ceil(x_circ(end)));
                        z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                        x_sol_r = linspace(x_circ(1),res_h,ceil(res_h-x_circ(1)));
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
                        for k=1:length(g)
                            if (xn(k)>=0) && (xn(k)<res_h) && (zn(k)>=0) && (zn(k)<res_v)
                                if double(I(floor(zn(k))+1,floor(xn(k))+1)) == 0
                                    I(floor(zn(k))+1,floor(xn(k))+1) = round(g(k));
                                else
                                    I(floor(zn(k))+1,floor(xn(k))+1) = round((double(I(floor(zn(k))+1,floor(xn(k))+1)) + g(k))/2);
                                end
                            end
                        end
                        %Iorig = I;
                        %- Filling background
                        for k = 1:(res_v-hsubstrato) %Search from left to right
                            for m = 1:res_h
                                if I(k,m) == 0
                                    I(k,m) = g2;
                                else
                                    break;
                                end
                            end
                        end
                        for k = 1:(res_v-hsubstrato) %Search from right to left
                            for m = res_h:-1:1
                                if I(k,m) == 0
                                    I(k,m) = g2;
                                else
                                    break;
                                end
                            end
                        end
                        %- Filling solid substrate,drop and needle
                        for k = 1:res_v %Search in all the image
                            for m = 1:res_h
                                if I(k,m) == 0
                                    I(k,m) = g1;
                                end
                            end
                        end
                        %Save images
                        imagename = strcat('Grayscale_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
                        fullimagename = strcat(imagename,imFmt); %Image name + image format
                        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                        imwrite(I,fullfilename) %Save Image
                    end
                %--------------------------------------------------------------
                %--------- SIMULATION OF NON-UNIFORM ILLUMINATION -------------
                %--------------------------------------------------------------
                elseif option == 7 %Simulation of non-uniform illumination
                    %- Simulation of non-uniform illumination
                    iref = res_v/2; %Coordenada no eixo i do centro da iluminação não-uniforme
                    jref = res_h/2; %Coordenada no eixo j do centro da iluminação não-uniforme
                    for n = 1:length(deltag)
                        %Creation of edge and normal profiles
                        if theta(j) <= 90 %theta<=90°
                            zc = res_v - hsubstrato + b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        else %theta>90°
                            zc = res_v - hsubstrato - b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        end
                        %zc = res_v-hsubstrato-b_pxs(j,i); %Coordenada y do centro da gota em pixels
                        [x,z] = circle(xc,zc,r_pxs(j,i),theta(j),Ncircle); %Spherical cap coordinates in pixels
                        %Flip drop around x-axis that pass through zc
                        x_circ = x;
                        z_circ = 2*zc-z;
                        %Creation of spherical cap drop and needle edge normal profile
                        if needle_in
                            %Construction of the needle profile
                            needle_diameter_pxs = (needleDiam/10)*scale_cm;
                            pos_x_needle_l = xc-(needle_diameter_pxs/2); %x coordinate of the left side of the needle
                            pos_x_needle_r = xc+(needle_diameter_pxs/2); %x coordinate of the right side of the needle
                            if theta(j)>90
                                indx_l = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                                indx_l = length(x_circ)/4 + indx_l;
                                indx_r = find((abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)) == min(abs(x_circ(length(x_circ)/4:3*length(x_circ)/4)-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                                indx_r = length(x_circ)/4 + indx_r;
                            else
                                indx_l = find((abs(x_circ-pos_x_needle_l)) == min(abs(x_circ-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                                indx_r = find((abs(x_circ-pos_x_needle_r)) == min(abs(x_circ-pos_x_needle_r)),1,'first'); %Find the index of the instersection between the right side of the needle and the drop
                            end
                            %indx_l = find((abs(x_circ-pos_x_needle_l)) == min(abs(x_circ-pos_x_needle_l)),1,'first'); %Find the index of the instersection between the left side of the needle and the drop
                            %indx_r = find((abs(x_circ-pos_x_needle_r)) == min(abs(x_circ-pos_x_needle_r)),1,'last'); %Find the index of the instersection between the right side of the needle and the drop
                            z_needle_max_l = z_circ(indx_l);
                            z_needle_max_r = z_circ(indx_r);
                            z_needle_l = linspace(0,z_needle_max_l,ceil(z_needle_max_l)); %z coordinates of the needle profile
                            z_needle_r = linspace(0,z_needle_max_r,ceil(z_needle_max_r)); %z coordinates of the needle profile
                            x_needle_l = ones(1,length(z_needle_l))*pos_x_needle_l; %x coordinates of the left needle profile
                            x_needle_r = ones(1,length(z_needle_r))*pos_x_needle_r; %x coordinates of the right needle profile
                            [xn_needle_l,zn_needle_l,g_needle_l] = grayscale_edge_nonuniform_illumination(x_needle_l,z_needle_l,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the left needle profile
                            [xn_needle_r,zn_needle_r,g_needle_r] = grayscale_edge_nonuniform_illumination(x_needle_r,z_needle_r,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the right needle profile
                            %g_needle_r = wrev(g_needle_r); %Fix grayscale order in the right profile of the needle
                            for k = 1:(length(g_needle_r)/n_normal) %Correcting g_sol_r order
                                prov = g_needle_r(n_normal*(k-1)+1:n_normal*k);
                                for m = 1:floor(n_normal/2)
                                    g_needle_r((n_normal*(k-1))+m) = prov((n_normal)-(m-1));
                                    g_needle_r((n_normal*k)-(m-1)) = prov(m);
                                end
                            end
                            xn_needle = [xn_needle_l xn_needle_r]; %x coordinates of edge and normal points of the needle profile
                            zn_needle = [zn_needle_l zn_needle_r]; %z coordinates of edge and normal points of the needle profile
                            g_needle = [g_needle_l g_needle_r]; %grayscale value of edge and normal points of the needle profile

                            %Trim drop profile based on the needle profile
                            x_drop_l = x_circ(indx_l:end);
                            z_drop_l = z_circ(indx_l:end);
                            x_drop_r = x_circ(1:indx_r);
                            z_drop_r = z_circ(1:indx_r);
                            [xn_drop_l,zn_drop_l,g_drop_l] = grayscale_edge_nonuniform_illumination(x_drop_l,z_drop_l,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the left portion of the drop
                            %g_drop_l = wrev(g_drop_l);
                            [xn_drop_r,zn_drop_r,g_drop_r] = grayscale_edge_nonuniform_illumination(x_drop_r,z_drop_r,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the right portion of the drop
                            %g_drop_r = wrev(g_drop_r);
                            xn_drop = [xn_drop_l xn_drop_r];
                            zn_drop = [zn_drop_l zn_drop_r];
                            g_drop = [g_drop_l g_drop_r];
                        else
                            [xn_drop,zn_drop,g_drop] = grayscale_edge_nonuniform_illumination(x_circ,z_circ,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the drop profile
                            %g_drop = wrev(g_drop);
                        end
                        %Creation of solid subratre edge normal profile
                        x_sol_l = linspace(0,x_circ(end),ceil(x_circ(end)));
                        z_sol_l = ones(1,length(x_sol_l))*(res_v-hsubstrato);
                        x_sol_r = linspace(x_circ(1),res_h,ceil(res_h-x_circ(1)));
                        z_sol_r = ones(1,length(x_sol_r))*(res_v-hsubstrato);
                        x_sol_l = wrev(x_sol_l);
                        [xn_sol_l,zn_sol_l,g_sol_l] = grayscale_edge_nonuniform_illumination(x_sol_l,z_sol_l,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the left portion of the solid substrate surface
                        %g_sol_l = wrev(g_sol_l); %Fix grayscale order along the left solid substrate profile
                        [xn_sol_r,zn_sol_r,g_sol_r] = grayscale_edge_nonuniform_illumination(x_sol_r,z_sol_r,n_normal,g1,g2,W_pxs,resolution,deltag(n),jref,iref); %Calculation of the normal points to the right portion of the solid substrate surface
                        %g_sol_r = wrev(g_sol_r); %Fix grayscale order along the right solid substrate profile
                        for k = 1:(length(g_sol_r)/n_normal) %Correcting g_sol_r order
                            prov = g_sol_r(n_normal*(k-1)+1:n_normal*k);
                            for m = 1:floor(n_normal/2)
                                g_sol_r((n_normal*(k-1))+m) = prov((n_normal)-(m-1));
                                g_sol_r((n_normal*k)-(m-1)) = prov(m);
                            end
                        end
                        xn_sol = [wrev(xn_sol_l) xn_sol_r];
                        zn_sol = [wrev(zn_sol_l) zn_sol_r];
                        g_sol = [wrev(g_sol_l) g_sol_r];
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
                        for k=1:length(g)
                            if (xn(k)>=0) && (xn(k)<res_h) && (zn(k)>=0) && (zn(k)<res_v)
                                if double(I(floor(zn(k))+1,floor(xn(k))+1)) == 0
                                    I(floor(zn(k))+1,floor(xn(k))+1) = round(g(k));
                                else
                                    I(floor(zn(k))+1,floor(xn(k))+1) = round((double(I(floor(zn(k))+1,floor(xn(k))+1)) + g(k))/2);
                                end
                            end
                        end
                        %Iorig = I;
                        %- Filling background
                        for k = 1:(res_v-hsubstrato) %Search from left to right
                            for m = 1:res_h
                                if I(k,m) == 0
                                    r = sqrt((k-iref)^2+(m-jref)^2);
                                    I(k,m) = round(g2 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                    %I(k,m) = g2;
                                else
                                    break;
                                end
                            end
                        end
                        for k = 1:(res_v-hsubstrato) %Search from right to left
                            for m = res_h:-1:1
                                if I(k,m) == 0
                                    r = sqrt((k-iref)^2+(m-jref)^2);
                                    I(k,m) = round(g2 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                    %I(k,m) = g2;
                                else
                                    break;
                                end
                            end
                        end
                        %- Filling solid substrate,drop and needle
                        for k = 1:res_v %Search in all the image
                            for m = 1:res_h
                                if I(k,m) == 0
                                    r = sqrt((k-iref)^2+(m-jref)^2);
                                    I(k,m) = round(g1 + deltag(n)*(1-((12*r^2)/(res_h^2+res_v^2))));
                                    %I(k,m) = g1;
                                end
                            end
                        end
                        %Save images
                        imagename = strcat('Non_uniform_illumination_Grayscale_Spherical_cap_Needle_in_',num2str(needle_in),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(theta(j),'%.2f'),'_scale(pxs_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2),'_deltag_',num2str(round(deltag(n)))); %Image filename
                        fullimagename = strcat(imagename,imFmt); %Image name + image format
                        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
                        imwrite(I,fullfilename) %Save Image
                    end
                end
            end
            %Write excel file with important spherical drop profile properties
            fprintf('Exporting excel file containing important spherical cap properties... \n');
            filename_xls = strcat('Sim_spherical_cap_V&theta_Vmin(uL)_',num2str(Vmin,'%.2f'),'_Vmax(uL)_',num2str(Vmax,'%.2f'),'_thetamin(°)_',num2str(thetamin,'%.2f'),'_thetamax(°)_',num2str(thetamax,'%.2f'),'.xlsx'); %Excel file name
            fullfilename_xls = fullfile(save_directory,filename_xls);
            varnames = {'theta[°]','r[cm]','r[pxs]','a[pxs]','a[cm]','b[pxs]','b[cm]','h[pxs]','h[cm]','xCP_left[pxs]','yCP_left[pxs]','xCP_right[pxs]','yCP_right[pxs]'}; %Variable names in Excel file
            T = table(theta',r_cm(:,i),r_pxs(:,i),a_pxs(:,i),a_cm(:,i),b_pxs(:,i),b_cm(:,i),h_pxs(:,i),h_cm(:,i),xCPleft_pxs(:,i),yCPleft_pxs(:,i),xCPright_pxs(:,i),yCPright_pxs(:,i),'VariableNames',varnames); %Table creation
            sheet_name = strcat('V = ',num2str(V_uL(i)),' uL',' scale = ',num2str(scale),' pxs_mm'); %Sheet name of Excel file
            writetable(T,fullfilename_xls,'Sheet',sheet_name) %Write table to an excel file
        end
        fprintf('Simulation completed! \n');
        pause(2)
    end
    clc
    clear all
    close all
end