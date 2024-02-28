function simBinSphCapDrop_V_theta_rh()
%SIMBINSPHCAPDROP_V_THETA_RH Simulation of binary images of spherical cap
%drops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of binary images of spherical cap. It requires the user to
%   set: the minimum drop volume [uL](Vmin), the maximum drop volume
%   [uL](Vmax), the number of drop volume levels that is going to be 
%   simulated between Vmin and Vmax (nVol), the minimum contact angle [°]
%   (thetamin), the maximum contact angle [°](thetamax), the number of 
%   contact angle levels that is going to be simulated between thetamin and
%   thetamax (ntheta), the needle diameter [mm](needleDiam), the
%   possibility of simulating images with and without needle (needleon),
%   the image resolution [pixels](res_h and res_v) and the image scale
%   [pixels/mm](scale). The height of the solid substrate in the image is
%   fixed. The routine also allows the user to choose the directory where 
%   the images are going to be exported and the image format. An excel file
%   (.xlsx) is created containing important properties of the simulated 
%   spherical drops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ SPHERICAL CAP DROP - SIMULATION OF BINARY IMAGES --------------------\n');
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
    %- Save options
    fprintf('Selection of the directory to export files ... \n');
    save_directory = uigetdir('C:\'); %Save directory
    filename_xls = strcat('Binary_SphericalCap_(needle_',num2str(needle),'_rh_cm_',num2str(r_h_cm,'%.4f'),'_Vmin_uL_',num2str(Vmin),'_Vmax_uL_',num2str(Vmax),'_thetamin_',num2str(thetamin),'_thetamax_',num2str(thetamax),').xlsx'); %Name of Excel file
    %filename_xls = 'Sim_sessile_spherical_V&theta.xlsx'; %Nome da planilha excel criada
    fullfilename_xls = fullfile(save_directory,filename_xls);

    %- Variables initialization
    r_pxs = zeros(length(theta),length(V_uL)); %radius of the droplet's spawning sphere in pixels
    a_pxs = zeros(length(theta),length(V_uL)); %drop wetted radius in pixels
    b_pxs = zeros(length(theta),length(V_uL)); %distance from the drop center to the cutting level in pixels
    h_pxs = zeros(length(theta),length(V_uL)); %drop height in pixels
    r_cm = zeros(length(theta),length(V_uL));  %radius of the droplet's spawning sphere in cm
    a_cm = zeros(length(theta),length(V_uL)); %drop wetted radius in cm
    b_cm = zeros(length(theta),length(V_uL)); %distance from the drop center to the cutting level in cm
    h_cm = zeros(length(theta),length(V_uL)); %drop height in cm
    xCPleft_pxs = zeros(length(theta),length(V_uL)); %x coordinate of the left triple contact point
    yCPleft_pxs = zeros(length(theta),length(V_uL)); %y coordinate of the left triple contact point
    xCPright_pxs = zeros(length(theta),length(V_uL)); %x coordinate of the right triple contact point
    yCPright_pxs = zeros(length(theta),length(V_uL)); %y coordinate of the right triple contact point

    % Simulation of binary spherical cap images
    fprintf('Generating binary spherical cap images... \n');
    xc = round(res_h/2); %x coordinate of the spherical drop center in pxs
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
            r_pxs(j,i) = r_cm(j,i)*scale_cm; %in pixels
            a_pxs(j,i) = a_cm(j,i)*scale_cm; %in pixels
            b_pxs(j,i) = b_cm(j,i)*scale_cm; %in pixels
            h_pxs(j,i) = h_cm(j,i)*scale_cm; %in pixels

            xCPleft_pxs(j,i) = xc - a_pxs(j,i);
            xCPright_pxs(j,i) = xc + a_pxs(j,i);

            %Generate of synthetical image
            if theta(j) <= 90 %theta<=90°
                yc = res_v - hsubstrato + b_pxs(j,i); %y coordinate of the spherical drop center in pxs
                yCPleft_pxs(j,i) = yc - b_pxs(j,i);
                yCPright_pxs(j,i) = yc - b_pxs(j,i);
            else %theta>90°
                yc = res_v - hsubstrato - b_pxs(j,i); %y coordinate of the spherical drop center in pxs
                yCPleft_pxs(j,i) = yc + b_pxs(j,i);
                yCPright_pxs(j,i) = yc + b_pxs(j,i);
            end
            M = ones(res_v,res_h)*255; %Matrix with image resolution
            I = uint8(M); %Image generation
            I = insertShape(I,'FilledCircle',[xc yc r_pxs(j,i)],'Color','black','Opacity',1); %Insert circle that represents the drop
            I = insertShape(I,'FilledRectangle',[1 (res_v-hsubstrato) res_h hsubstrato],'Color','black','Opacity',1); %Insert rectangle that represents the solid substrate
            if needle
                needle_diameter_pxs = needleDiam*scale; %Needle outer diameter in pixels
                I = insertShape(I,'FilledRectangle',[(xc-needle_diameter_pxs/2) 1 needle_diameter_pxs yc],'Color','black','Opacity',1);%Insert rectangle that represents the needle
            end
            %Conversion to binary image
            Igray = rgb2gray(I); %Convert image to grayscale
            Ibw = imbinarize(Igray); %Binarize image

            %Save Image
            filename = strcat('Binary_SphericalCap_Needle_',num2str(needle),'_V(uL)_',num2str(V_uL(i),"%.2f"),'_theta(°)_',num2str(theta(j),"%.2f"),'_scale(pxs_mm)_',num2str(scale)); %Image filename
            fullimagename = strcat(filename,imFmt); %Image name + image format
            fullfilename = fullfile(save_directory,fullimagename); %Image fullfilename
            imwrite(Ibw,fullfilename) %Save Image

        end
        %Write excel file with important spherical cap properties
        fprintf('Exporting excel file containing important spherical cap properties... \n');
        varnames = {'theta[°]','r[cm]','r[pxs]','a[pxs]','a[cm]','b[pxs]','b[cm]','h[pxs]','h[cm]','xCP_left[pxs]','yCP_left[pxs]','xCP_right[pxs]','yCP_right[pxs]'}; %Variable names
        T = table(theta',r_cm(:,i),r_pxs(:,i),a_pxs(:,i),a_cm(:,i),b_pxs(:,i),b_cm(:,i),h_pxs(:,i),h_cm(:,i),xCPleft_pxs(:,i),yCPleft_pxs(:,i),xCPright_pxs(:,i),yCPright_pxs(:,i),'VariableNames',varnames); %Table creation
        sheet_name = strcat('V = ',num2str(V_uL(i)),' uL',' scale = ',num2str(scale),' pxs_mm'); %Excel file name
        writetable(T,fullfilename_xls,'Sheet',sheet_name) %Write table to an excel file
    end
    fprintf('Simulation completed! \n');
    pause(2)
end
end
