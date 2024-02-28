function simExpIncDrop()
%SIMEXPINCDROP Simulation of quasi-static inclined drop experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of quasi-static inclined drop experiment. Generation of
%   grayscale inclined drop images for a series of inclination angles 
%   (alpha) given liquid capillary constant (c), contact angle (CA) and 
%   drop volume (V). The routine requires the user to enter: liquid 
%   capillary constant [cm-2](c), minimum curvature at the apex possible
%   [cm-1](bmin), maximum curvature at the apex possible [cm-1](bmax), drop
%   volume [uL](V), contact angle [°](CA), initial inclination angle [°]
%   (incIni), final inclination angle [°](incMax), inclination rate [°/s], 
%   the possibility of simulating images with inclined or planed solid 
%   substrate (inclined), camera framerate [fps](framerate) and number of 
%   captured images (Nimages). The routine also allows the user to choose 
%   the directory where the images are going to be exported. An excel file 
%   (.xlsx) is created containing important properties of the inclined
%   drops. During the computation of inclined drop profiles errors may 
%   occur due to poor choices of input parameters. For certain combinations 
%   of curavture at the apex (bmin and bmax), contact angle (CA) and 
%   inclination angle (alpha), the inclined drop may not be simulated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ INCLINED DROP - SIMULATION OF QUASI-STATIC EXPERIMENT --------------------\n');
    resp = input("- Do you want to continue (Y/N) [Y]? ", "s");
    if isempty(resp) || (resp ~= 'Y' && resp ~= 'N')
        resp = 'Y';
    end
    if resp == 'N'
        break;
    end
    %----------------------------------------------------------------------
    %----------------------- INPUT PARAMETERS -----------------------------
    %----------------------------------------------------------------------
    %- INPUT PARAMETERS
    fprintf('Input parameters: \n');
    fprintf('1) Liquid propertes: \n');
    c = input("- Liquid capillary constant in cm-2 [e.g. water = 13.55]: "); %Cappilary constant in cm-2
    if isempty(c)
        c = 13.55; %Capillary constant of water in cm-2
    end
    fprintf('2) Drop configuration propertes: \n');
    CA = input("- Contact angle in degree [70]: "); %Contact angle in degree
    if isempty(CA)
        CA = 70; 
    end
    V = input("- Drop volume in uL [5]: "); %Drop volume in uL
    if isempty(V)
        V = 5; 
    end
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
    interval_b = 50; %Number of b values created between bmin and bmax for volume search
    tolV = 1e-3; %Tolerance for volume determination in uL
    itermax = 5; %Maximum number of iterations
    fprintf('3) Image propertes: \n');
    inclined = input("- Do you want to show the solid substrate inclined (Y/N)? [Y]: ", "s"); %Generate images with substrate inclined (inclined == 'Y') or planed (inclined == 'N')
    if isempty(inclined) || (inclined ~= 'Y' && inclined ~= 'N')
        inclined = 'Y';
    end
    if inclined == 'N'
        inclined = 0;
    else
        inclined = 1;
    end
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
    hsubstrato = 100; %Substrate height in pixels;
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
    fprintf('5) Quasi-static experiment propertes: \n');
    alphaIni = input("- Initial inclination angle in degree [0]: "); %Initial inclination angle in degree 
    if isempty(alphaIni)
        alphaIni = 0;
    end
    %}
    while 1
        alphaEnd = input("- Final inclination angle in degree [60]: "); %Final inclination angle in degree
        if isempty(alphaEnd)
            alphaEnd = 60;
        end
        if alphaEnd < alphaIni
            fprintf('Invalid final inclination angle. Please, enter a lower value. \n');
        else
            break;
        end
    end
    incRate = input("- Inclination rate in °/s [0.1]: ");
    if isempty(incRate)
        incRate = 0.1;
    end
    framerate = input("- Camera framerate in fps [30]: ");
    if isempty(framerate)
        framerate = 30;
    end
    totalframes = floor(framerate*(alphaEnd-alphaIni)/incRate); %Total number of frames acquired by the camera during the experiment
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
    alpha = linspace(alphaIni,alphaEnd,Nimages); %Inclination angle of each simulated drop image
    t_total = abs(alphaEnd - alphaIni)/incRate; %Total time of the experiment
    time = linspace(0,t_total,Nimages); %Time of each simulated drop image
    %Numerical integration options
    N = 13; %number of elements along theta/ Suggestion: give N odd numbers, so the middle element is given at theta=90°
    %- Save options
    fprintf('Selection of the directory to export files ... \n');
    save_directory = uigetdir('C:\'); %Save directory

    %- Variables initialization
    b_otm = zeros(length(alpha),1); %Curvature at the apex in cm-2
    r0_otm = zeros(length(alpha),1); %Radius at the apex in cm
    CA_sim_otm = zeros(length(alpha),1); %Intrinsic apparent CA for a W distance from the apex of the drop in degrees
    CA_adv_otm = zeros(length(alpha),1); %Advancing CA in degrees
    CA_rec_otm = zeros(length(alpha),1); %Receding CA in degrees
    Xcp_left = zeros(length(alpha),1); %X coordinate of the left triple contact point in pxs
    Ycp_left = zeros(length(alpha),1); %Y coordinate of the left triple contact point in pxs
    Xcp_right = zeros(length(alpha),1); %X coordinate of the right triple contact point in pxs
    Ycp_right = zeros(length(alpha),1); %Y coordinate of the right triple contact point in pxs
    d_w_cm_otm = zeros(length(alpha),1); %Distance between contact points in cm
    h_cm_otm = zeros(length(alpha),1); %Drop height in cm
    A_cm2_otm = zeros(length(alpha),1); %Drop surface area in cm²
    V_uL_otm = zeros(length(alpha),1); %Drop volume in mm³ or uL
        
    %----------------------------------------------------------------------
    %----------------- SIMULATION OF INCLINED DROP IMAGES -----------------
    %----------------------------------------------------------------------
    %- Generation of grayscale inclined drop images
    for u = 1:length(alpha)
        fprintf('-------------------------------------------------------------------------------------------------------- \n');
        fprintf('Simulating inclined drop profile for c = %.2f cm-2, alpha = %.2f°, CA = %.2f° and V = %.2f uL ... \n',c,alpha(u),CA,V);
        %------------------------------------------------------------------
        %--------- FIND B VALUE THAT CORRESPONDS TO THE DESIRED V ---------
        %------------------------------------------------------------------
        % Variable initialization
        x_ly_sim_cm = cell(1,interval_b);
        z_ly_sim_cm = cell(1,interval_b);
        CA_sim = zeros(1,interval_b); %Intrinsic apparent CA for a W distance from the apex of the drop in degrees
        CA_adv = zeros(1,interval_b); %Advancing CA in degrees
        CA_rec = zeros(1,interval_b); %Receding CA in degrees
        d_w_cm = zeros(1,interval_b); %Distance between contact points in cm
        h_cm = zeros(1,interval_b); %Drop height in cm
        A_cm2 = zeros(1,interval_b); %Drop surface area in cm²
        V_uL = zeros(1,interval_b); %Drop volume in mm³ or uL

        % Search for b value
        diffV = 1; %reset volume difference
        bmin = bmin_ini; %reset bmin
        bmax = bmax_ini; %reset bmax
        iter = 1;
        while diffV > tolV && iter <= itermax
            fprintf('Searching for V = %.2f uL, between bmin = %.4f and bmax = %.4f (iteration %d)... \n',V,bmin,bmax,iter);
            b = linspace(bmin,bmax,interval_b);
            r0 = 1./b;
            for i = 1:length(b)
                %Determination of the inclined drop profile by using numerical solution technique for solving Young-Laplace PDE
                [x_ly_sim_cm{i},z_ly_sim_cm{i},W,U,theta0,CA_sim(i),CA_adv(i),CA_rec(i),d_w_cm(i),h_cm(i),A_cm2(i),V_uL(i)] = profile_inclined_drop_ly(b(i),c,alpha(u),CA,N);
            end
            diffV = min(abs(V_uL-V));
            fprintf('diffV = %f\n',diffV);
            indx_V = find(abs(V_uL-V) == diffV,1,'last');
            if indx_V-1 <= 0
                bmin = b(indx_V);
                bmax = b(indx_V+1);
            elseif indx_V+1 >= length(b)
                bmin = b(indx_V-1);
                bmax = b(indx_V);
            else
                bmin = b(indx_V-1);
                bmax = b(indx_V+1);
            end
            iter = iter +1;
        end
        V_uL_otm(u) = V_uL(indx_V);
        x_ly_sim_cm_otm = x_ly_sim_cm{indx_V};
        z_ly_sim_cm_otm = z_ly_sim_cm{indx_V};
        b_otm(u) = b(indx_V);
        r0_otm(u) = r0(indx_V);
        CA_sim_otm(u) = CA_sim(indx_V);
        CA_adv_otm(u) = CA_adv(indx_V);
        CA_rec_otm(u) = CA_rec(indx_V);
        d_w_cm_otm(u) = d_w_cm(indx_V);
        h_cm_otm(u) = h_cm(indx_V);
        A_cm2_otm(u) = A_cm2(indx_V);

        %----------------------------------------------------------------------
        %------------------ GENERATION OF GRAYSCALE IMAGES --------------------
        %----------------------------------------------------------------------
        fprintf('Generating grayscale drop image... \n');
        %for k = 1:length(b)
        %fprintf('----------------------------------------------------------------------------------------- \n');
        %{
        fprintf('Simulating inclined drop profile for c = %.2f cm-2, alpha = %.2f°, CA = %.2f° and V = %.2f uL ... \n',c,alpha_deg,mean(CA),V(k));
        %Determination of the inclined drop profile by using numerical solution technique for solving Young-Laplace PDE
        [x_ly_sim_cm,z_ly_sim_cm,W,U,theta0,CA_sim(k),CA_adv(k),CA_rec(k),d_w_cm(k),h_cm(k),A_cm2(k),V_uL(k)] = profile_inclined_drop_ly(b(k),c,alpha_deg,CA(k),N);
        %}
        %Conversion to pixels
        x_ly_sim_pxs = x_ly_sim_cm_otm*scale_cm;   %Drop width in pixels
        z_ly_sim_pxs = z_ly_sim_cm_otm*scale_cm;   %Drop height in pixels
        %Construction of synthetical image
        res_v = resolution(2);
        res_h = resolution(1);
        xc = round(res_h/2); %Coordenada x do centro da imagem em pixels
        zc = round(res_v/2); %Coordenada z do centro da imagem em pixels

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
            [xn,zn] = rotate_data(xn,zn,xc,zc,alpha(u));
        end

        %Determination of triple contact points
        if inclined
            [x_plot_rot,z_plot_rot] = rotate_data(x_plot,z_plot,xc,zc,alpha(u));
            Xcp_left(u) = x_plot_rot(1);
            Ycp_left(u) = z_plot_rot(1);
            Xcp_right(u) = x_plot_rot(end);
            Ycp_right(u) = z_plot_rot(end);
        else
            Xcp_left(u) = x_plot(1);
            Ycp_left(u) = z_plot(1);
            Xcp_right(u) = x_plot(end);
            Ycp_right(u) = z_plot(end);
        end

        %Generation of synthetical illuminated (grayscale) inclined drop images
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
            [x_plot_r_rot,z_plot_r_rot] = rotate_data(x_plot_r,z_plot_r,xc,zc,alpha(u));
            x_plot_l = x_plot(1:indx_apex);
            z_plot_l = z_plot(1:indx_apex);
            [x_plot_l_rot,z_plot_l_rot] = rotate_data(x_plot_l,z_plot_l,xc,zc,alpha(u));
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
        %Save image
        imagename = strcat('Grayscale_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_CA(°)_',num2str(CA_sim_otm(u),'%.2f'),'_alpha(°)_',num2str(alpha(u),'%.2f'),'_b(cm-1)_',num2str(b_otm(u),'%.3f'),'_V(uL)_',num2str(V_uL_otm(u),'%.2f'),'_CA_max(°)_',num2str(CA_adv_otm(u),'%.2f'),'_CA_min(°)_',num2str(CA_rec_otm(u),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm),'_W(pxs)_',num2str(W_pxs,'%.4f'),'_g1(pxs)_',num2str(g1),'_g2(pxs)_',num2str(g2)); %Image filename
        fullimagename = strcat(imagename,imFmt); %Image name + image format
        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
        imwrite(I,fullfilename) %Save Image
    end
    %Write excel file with important inclined drop profile properties
    fprintf('Exporting inclined drop properties... \n');
    varnames = {'time[s]','alpha[°]','b[cm-1]','r0[cm]','CA[°]','CA_max[°]','CA_min[°]','Xcp_left[pxs]','Ycp_left[pxs]','Xcp_right[pxs]','Ycp_right[pxs]','d_w[cm]','h[cm]','A[cm2]','V[uL]'}; %Variables name in the excel file
    T = table(time',alpha',b_otm,r0_otm,CA_sim_otm,CA_adv_otm,CA_rec_otm,Xcp_left,Ycp_left,Xcp_right,Ycp_right,d_w_cm_otm,h_cm_otm,A_cm2_otm,V_uL_otm,'VariableNames',varnames); %Build a table with important drop properties
    sheet_name = strcat('c = ',num2str(c),' cm^-2'); %Excel sheet name
    savename_xls = strcat('Dimensions_InclinedDrop_',num2str(inclined),'(c(cm-2)_',num2str(c),'_CA(°)_',num2str(mean(CA_sim_otm(u)),'%.2f'),'_alphaIni(°)_',num2str(alphaIni),'_alphaEnd(°)_',num2str(alphaEnd),'_incRate(uL_s)_',num2str(incRate),').xlsx'); %Excel file name
    fullsavename_xls = fullfile(save_directory,savename_xls);
    writetable(T,fullsavename_xls,'Sheet',sheet_name) %Write table to an excel file

    fprintf('Simulation completed! \n');
    pause(2)
end
end