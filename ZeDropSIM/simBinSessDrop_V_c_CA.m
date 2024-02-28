function simBinSessDrop_V_c_CA()
%SIMBINSPHCAPDROP_V_C_CA Simulation of binary sessile drop images for 
% different drop volume (V) given liquid capillary constant (c) and contact
% angle (CA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of binary sessile drop images for different drop volumes (V)
%   given liquid capillary constant (c) and contact angle (CA). 
%   It requires the user to set: the liquid capillary constant [cm-2](c), 
%   the minimum drop volume [uL](Vmin), the maximum drop volume [uL](Vmax),
%   the number of drop volume levels that is going to be simulated between 
%   Vmin and Vmax (nVol), the contact angle [°](theta), the needle diameter
%   [mm](needleDiam), the possibility of simulating images with and without 
%   needle (needleon), the image resolution [pixels](res_h and res_v) and 
%   the image scale [pixels/mm](scale). In order to make the execution of 
%   the algorithm faster, it allows the user to load an excel file 
%   containing important drop properties such as curvature at the apex (b),
%   drop volume (V) and contact angle (CA). The height of the solid 
%   substrate in the image is fixed. The routine also allows the user to 
%   choose the directory where the images are going to be exported and the 
%   image format. An excel file (.xlsx) is created containing important 
%   properties of the simulated sessile drops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ SESSILE DROP - SIMULATION OF BINARY IMAGES --------------------\n');
    resp = input("- Do you want to continue (Y/N) [Y]? ", "s");
    if isempty(resp) || (resp ~= 'Y' && resp ~= 'N')
        resp = 'Y';
    end
    if resp == 'N'
        break;
    end
    % Input parameters
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
    scale_cm = scale*10; %Escala em pixels/cm
    imFmt = input("- Image format [.tif]: ", "s");
    if isempty(imFmt)
        imFmt = ".tif";
    end
    %Numerical integration options
    S_span = (0:.0001:8); % S_span is the step variable for ode45 solver
    initialCo = [0 1e-100 0 0 0]; %Initial conditions for integration
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
        
    %- Generation of binary sessile drop images
    %fprintf('Generating binary sessile drop images... \n');
    for i = 1:nSimImages
        fprintf('----------------------------------------------------------------------------------------- \n');     
        fprintf('Simulating sessile drop profile for c = %.2f cm-2, CA = %.2f° and V = %.2f uL ... \n',c,mean(theta),V(i));
        if loadExc == 'Y'   
            %Determination of sessile drop profile by numerical integration of the Young-Laplace equation - given c, b and CA values
            [x_ly_sim_cm,z_ly_sim_cm,theta_sim(i),r_w_cm(i),r_eq_cm(i),h_cm(i),A_cm2(i),V_uL(i)] = profile_sessile_drop_ly(b(i),c,theta(i),S_span,initialCo);
        else
            %Determination of sessile drop profile by numerical integration of the Young-Laplace equation - given c, V and CA values
            [x_ly_sim_cm,z_ly_sim_cm,theta_sim(i),b(i),r_w_cm(i),r_eq_cm(i),h_cm(i),A_cm2(i),V_uL(i)] = profile_sessile_drop_ly_V(V(i),c,theta,S_span,initialCo,bmin_ini,bmax_ini,interval_b,tolV,itermax);
        end
        %Creation of synthetical binary image
        fprintf('Generation of binary sessile drop image ... \n');
        [Ibw] = binary_sessile_drop_ly(x_ly_sim_cm,z_ly_sim_cm,hsubstrato,resolution,scale_cm,needleDiam,needle_in);

        %Determination of the coordinates of triple contact points
        x_ly_sim_pxs = x_ly_sim_cm*scale_cm;   % Drop width in pixels
        z_ly_sim_pxs = z_ly_sim_cm*scale_cm;   % Drop height in pixels
        xc = round(res_h/2); %x coordinate of drop center in pixels
        x_plot_lap_full = [wrev(-x_ly_sim_pxs+xc); x_ly_sim_pxs+xc]; %Creation of drop profile (x coordinate)
        z_plot = z_ly_sim_pxs - (hsubstrato + z_ly_sim_pxs(end) - res_v);
        z_plot_lap_full = [wrev(z_plot);z_plot]; %Creation of drop profile (z coordinate)
        xCPleft_pxs(i) = x_plot_lap_full(1);
        xCPright_pxs(i) = x_plot_lap_full(end);
        yCPleft_pxs(i) = z_plot_lap_full(1);
        yCPright_pxs(i) = z_plot_lap_full(end);

        %Save image
        fprintf('Saving image... \n');
        imagename = strcat('Binary_Sessile_drop_Needle_in_',num2str(needle_in),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(i),'%.3f'),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_theta(°)_',num2str(mean(theta),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm)); %Image filename
        fullimagename = strcat(imagename,imFmt); %Image name + image format
        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
        imwrite(Ibw,fullfilename) %Save Image
    end
    r0 = 1./b; %radius at the apex in cm
    %Write excel file with important spherical drop profile properties
    fprintf('Exporting excel file containing important sessile drop properties... \n');
    varnames = {'b[cm-1]','r0[cm]','theta[°]','xCP_left[pxs]','yCP_left[pxs]','xCP_right[pxs]','yCP_right[pxs]','r_w[cm]','r_eq[cm]','h[cm]','A [cm2]','V[uL]'}; %Variables name in the excel file
    T = table(b,r0,theta_sim,xCPleft_pxs,yCPleft_pxs,xCPright_pxs,yCPright_pxs,r_w_cm,r_eq_cm,h_cm,A_cm2,V_uL,'VariableNames',varnames); %Build a table with important drop properties
    savename_xls = strcat('Binary_Sessiledrop_(needle_in_',num2str(needle_in),'_c(cm-2)_',num2str(c),'_theta(°)_',num2str(mean(theta),'%.2f'),').xlsx'); %Excel file name
    fullsavename_xls = fullfile(save_directory,savename_xls);
    writetable(T,fullsavename_xls) %Write table to an excel file

    fprintf('Simulation completed! \n');
    pause(2)
end
end