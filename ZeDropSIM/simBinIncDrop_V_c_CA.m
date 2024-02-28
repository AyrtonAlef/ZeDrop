function simBinIncDrop_V_c_CA()
%SIMBININCDROP_V_C_CA Simulation of inclined sessile drop images for 
% different drop volume (V) given liquid capillary constant (c), contact
% angle (CA) and inclination angle (alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of binary inclined drop images for different drop volumes 
%   (V) given liquid capillary constant (c), inclination angle (alpha) and 
%   contact angle (CA). It requires the user to enter: the liquid capillary 
%   constant [cm-2](c), the inclination angle [°](CA), the possibility of 
%   simulating images inclined or planed solid substrate (inclined), the 
%   image resolution [pixels](res_h and res_v) and the image scale 
%   [pixels/mm](scale). In order to make the execution of the algorithm 
%   faster, it allows the user to load an excel file containing essential 
%   drop properties: curvature at the apex (b), drop volume (V) and contact
%   angle (CA). The height of the solid substrate in the image is fixed. 
%   The routine also allows the user to choose the directory where the 
%   images are going to be exported and the image format. An excel file 
%   (.xlsx) is created containing important properties of the simulated 
%   sessile drops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ INCLINED DROP - SIMULATION OF BINARY IMAGES --------------------\n');
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
    % Input parameters
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
    scale_cm = scale*10; %Escala em pixels/cm
    imFmt = input("- Image format [.tif]: ", "s");
    if isempty(imFmt)
        imFmt = ".tif";
    end
    %Numerical integration options
    N = 13; %number of elements along theta/ Suggestion: give n odd numbers, so the middle element is given at theta=90°

    %- Save options
    fprintf('Selection of the directory to export files ... \n');
    save_directory = uigetdir('C:\'); %Save directory
   
    %- Variables initialization
    CA_sim = zeros(length(b),1); %Intrinsic apparent CA for a W distance from the apex of the drop in degrees
    CA_adv = zeros(length(b),1); %Advancing (downhill) CA in degrees
    CA_rec = zeros(length(b),1); %Receding (uphill) CA in degrees
    d_w_cm = zeros(length(b),1); %Distance between contact points in cm
    h_cm = zeros(length(b),1); %Drop height in cm
    A_cm2 = zeros(length(b),1); %Drop surface area in cm²
    V_uL = zeros(length(b),1); %Drop volume in mm³ or uL
    
    %----------------------------------------------------------------------
    %----------------- SIMULATION OF INCLINED DROP IMAGES -----------------
    %----------------------------------------------------------------------
    %- Generation of binary inclined drop images
    for i = 1:length(b)
        fprintf('----------------------------------------------------------------------------------------- \n');     
        fprintf('Simulating binary drop profile for c = %.2f cm-2, alpha = %.2f°, CA = %.2f° and V = %.2f uL ... \n',c,alpha_deg,mean(CA),V(i));
        if loadExc == 'Y'   
            %Determination of the inclined drop profile by using numerical solution technique for solving Young-Laplace PDE
            [x_ly_sim_cm,z_ly_sim_cm,W,U,theta0,CA_sim(i),CA_adv(i),CA_rec(i),d_w_cm(i),h_cm(i),A_cm2(i),V_uL(i)] = profile_inclined_drop_ly(b(i),c,alpha_deg,CA(i),N);
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
                fprintf('Searching for V = %.2f uL, between bmin = %.4f and bmax = %.4f (iteration %d)... \n',V(i),bmin,bmax,iter);
                b_prov = linspace(bmin,bmax,interval_b);
                r0_prov = 1./b_prov;
                for u = 1:length(b_prov)
                    %Determination of the inclined drop profile by using numerical solution technique for solving Young-Laplace PDE
                    [x_ly_sim_cm_prov{u},z_ly_sim_cm_prov{u},W,U,theta0,CA_sim_prov(u),CA_adv_prov(u),CA_rec_prov(u),d_w_cm_prov(u),h_cm_prov(u),A_cm2_prov(u),V_uL_prov(u)] = profile_inclined_drop_ly(b_prov(u),c,alpha_deg,CA,N);
                end
                diffV = min(abs(V_uL_prov-V(i)));
                fprintf('diffV = %f\n',diffV);
                indx_V = find(abs(V_uL_prov-V(i)) == diffV,1,'last');
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
            V_uL(i) = V_uL_prov(indx_V);
            x_ly_sim_cm = x_ly_sim_cm_prov{indx_V};
            z_ly_sim_cm = z_ly_sim_cm_prov{indx_V};
            b(i) = b_prov(indx_V);
            r0(i) = r0_prov(indx_V);
            CA_sim(i) = CA_sim_prov(indx_V);
            CA_adv(i) = CA_adv_prov(indx_V);
            CA_rec(i) = CA_rec_prov(indx_V);
            d_w_cm(i) = d_w_cm_prov(indx_V);
            h_cm(i) = h_cm_prov(indx_V);
            A_cm2(i) = A_cm2_prov(indx_V);
        end

        %Creation of synthetical binary image
        fprintf('Generation of binary sessile drop image ... \n');
        [Ibw] = binary_inclined_drop_ly(x_ly_sim_cm,z_ly_sim_cm,hsubstrato,resolution,scale_cm,inclined,alpha_deg);

        %Save image
        fprintf('Saving image... \n');
        %imagename = strcat('Binary_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(i),'%.3f'),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(mean(theta),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm)); %Image filename
        imagename = strcat('Binary_InclinedDrop_Inclined_',num2str(inclined),'_c(cm^-2)_',num2str(c),'_CA(°)_',num2str(CA_sim(i),'%.2f'),'_alpha(°)_',num2str(alpha_deg),'_b(cm-1)_',num2str(b(i),'%.3f'),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_CA_max(°)_',num2str(CA_adv(i),'%.2f'),'_CA_min(°)_',num2str(CA_rec(i),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm)); %Image filename
        fullimagename = strcat(imagename,imFmt); %Image name + image format
        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
        imwrite(Ibw,fullfilename) %Save Image
    end
    %Write excel file with important inclined drop profile properties
    fprintf('Exporting excel file containing important inclined drop properties... \n');
    %varnames = {'b[cm-1]','r0[cm]','theta[°]','xCP_left[pxs]','yCP_left[pxs]','xCP_right[pxs]','yCP_right[pxs]','r_w[cm]','r_eq[cm]','h[cm]','A [cm-2]','V[uL]'}; %Vairables name in the excel file
    %T = table(b,r0,theta_sim,xCPleft_pxs,yCPleft_pxs,xCPright_pxs,yCPright_pxs,r_w_cm,r_eq_cm,h_cm,A_cm2,V_uL,'VariableNames',varnames); %Build a table with important drop properties
    varnames = {'b[cm-1]','r0[cm]','CA[°]','CA_max[°]','CA_min[°]','d_w[cm]','h[cm]','A[cm2]','V[uL]'}; %Variables names in the excel file
    T = table(b,r0,CA_sim,CA_adv,CA_rec,d_w_cm,h_cm,A_cm2,V_uL,'VariableNames',varnames); %Build a table with important drop properties
    sheet_name = strcat('c = ',num2str(c),' cm^-2',' alpha = ',num2str(alpha_deg),' °'); %Excel sheet name
    savename_xls = strcat('Binary_Inclined_drop(c(cm-2)_',num2str(c),'_CA(°)_',num2str(mean(CA_sim(i)),'%.2f'),'_alpha(°)_',num2str(alpha_deg),').xlsx'); %Excel file name
    fullsavename_xls = fullfile(save_directory,savename_xls);
    writetable(T,fullsavename_xls,'Sheet',sheet_name) %Write table to an excel file

    fprintf('Simulation completed! \n');
    pause(2)
end
end