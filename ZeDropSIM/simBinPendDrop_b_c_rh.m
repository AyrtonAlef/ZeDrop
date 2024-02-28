function simBinPendDrop_b_c_rh()
%SIMBINPENDDROP_B_C_RH Simulation of binary images of pendant drop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of binary images of pendant drop. It requires the user to 
%   set: the liquid capillary constant [cm-2](c), the minimum curvature at 
%   the apex [cm-1](bmin), the maximum curvature at the apex [cm-1](bmax), 
%   the number of pendant drop images that is going to be simulated between
%   bmin and bmax (nProfiles) and the needle diameter [mm](needleDiam). The
%   user may generate pendant drop images with and without needle. The 
%   pendant drop is centered in the image. The routine also allows the user
%   to choose the directory where the images are going to be exported and 
%   the image format. An excel file (.xlsx) is created containing important
%   properties of the pendant drops. During the computation of pendant drop
%   profiles errors may occur due to poor choices of input parameters. For 
%   certain combinations of curavture at the apex (bmin and bmax) and 
%   needle diameter, the pendant drop may not be simulated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ PENDANT DROP - SIMULATION OF BINARY IMAGES --------------------\n');
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
    bmin = input("- Minimum curvature at the apex in cm-1 [7.2]: "); %Minimum curvature at the apex in cm-1
    if isempty(bmin)
        bmin = 7.2;
    end
    while 1
        bmax = input("- Maximum curvature at the apex in cm-1 [9.7]: "); %Maximum curvature at the apex in cm-1;
        if isempty(bmax)
            bmax = 9.7;
        end
        if bmax < bmin
            fprintf('Invalid maximum curavture at the apex. Please, set a lower value. \n');
        else
            break;
        end
    end
    nb = input("- Number of curvature at the apex level simulated [3]: "); %Number of simulated images between minimum and maximum curvature at the apex
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
    r_h_cm = (needleDiam/2)/10; %Holder outer radius in cm
    needleON = input("- Needle ON (Y/N)? [Y]: ", "s"); %Generate images with (needle == 'Y') or without needle (needle == 'N')
    if isempty(needleON) || (needleON ~= 'Y' && needleON ~= 'N')
        needleON = 'Y';
    end
    if needleON == 'N'
        needle = 0;
    else
        needle = 1;
    end
    tol = 10^-5; %Tolerance to find the intersection points of the drop with the holder radius
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
    fprintf('Selection of the directory to export files ... \n');
    save_directory = uigetdir('C:\'); %Save directory

    % Variables calculation and initialization
    S_span=(0:.0001:8); % S_span is the step variable for ode45 solver
    initialCo = [0 1e-100 0 0 0]; %Initial conditions for integration

    d_eq_cm = zeros(length(b),length(c)); %Maximum drop diameter in cm
    h_eq_cm = zeros(length(b),length(c)); %Distance from the apex of the drop to its maximum diameter (dop height) in cm
    h_cm = zeros(length(b),length(c)); %Distance from the apex of the drop to the needle tip in cm
    d_s_cm = zeros(length(b),length(c)); %Drop diameter from a distance equivalent to the maximum drop diameter from the apex of the drop in cm
    A_cm2 = zeros(length(b),length(c)); %Drop surface area in cm²
    V_uL = zeros(length(b),length(c)); %Drop volume in uL calculated concomitantly with the Laplace curve generation

    % Simulation of binary pendant drop images
    fprintf("Simulation of binary pendant drop images... \n");
    for i=1:length(b)
        %Determination of the pendant drop profile
        [x_ly_cm,z_ly_cm,x_ly_sim_cm,z_ly_sim_cm,indx_last,indx_second,indx_deq,indx_ds,d_eq_cm(i),h_eq_cm(i),h_cm(i),d_s_cm(i),A_cm2(i),V_uL(i)] = profile_pendant_drop_ly(b(i),c,r_h_cm,S_span,initialCo,tol);

        %Generation of synthetical binary pendant drop image
        [Ibw] = binary_pendant_drop_ly(x_ly_sim_cm,z_ly_sim_cm,r_h_cm,resolution,scale_cm,needle);

        %Save image
        %imagename = strcat('Pendant_drop_needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(i),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_V(uL)_',num2str(V_uL(i),'%.2f'),'.tif'); %Image filename
        imagename = strcat('Binary_Pendant_drop_Needle_',num2str(needle),'_c(cm^-2)_',num2str(c),'_b(cm-1)_',num2str(b(i),'%.3f'),'_r_h(cm)_',num2str(r_h_cm,'%.4f'),'_V(uL)_',num2str(V_uL(i),'%.2f'),'_scale(pixels_cm)_',num2str(scale_cm)); %Image filename
        fullimagename = strcat(imagename,imFmt); %Image name + image format
        fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
        imwrite(Ibw,fullfilename) %Save Image
        %}
    end
    %Write excel file with important pendant drop profile properties
    varnames = {'b[cm-1]','r0[cm]','De[cm]','h_De[cm]','h[cm]','Ds[cm]','A[cm²]','V[uL]'}; %Variable names
    T = table(b',r0',d_eq_cm,h_eq_cm,h_cm,d_s_cm,A_cm2,V_uL,'VariableNames',varnames);
    sheet_name = strcat('c = ',num2str(c,'%.2f'),' cm-2',' r_h = ',num2str(r_h_cm,'%.4f'),' cm'); %Sheet name
    savename_xls = strcat('Binary_Pendant_drop_(needle_',num2str(needle),'_c_cm-2_',num2str(c),'_r_h_cm_',num2str(r_h_cm,'%.4f'),').xlsx'); %Excel file name
    fullsavename_xls = fullfile(save_directory,savename_xls);
    writetable(T,fullsavename_xls,'Sheet',sheet_name) %Write table to an excel file

    fprintf('Simulation completed! \n');
    pause(2)
end
end

