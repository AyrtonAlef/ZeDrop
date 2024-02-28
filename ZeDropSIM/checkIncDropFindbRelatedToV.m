function checkIncDropFindbRelatedToV()
%CHECKINCDROPFINDBRELATEDTOV Find curvatures at the apex (b) that 
% corresponds to desired drop volumes (V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Find curvatures at the apex (b) that corresponds to desired drop
%   volumes (V). It requires the uer to enter: the liquid capillary 
%   constant [cm-2](c), the minimum curvature at the apex level that is 
%   going to be searched [cm-1](bmin_ini), the maximum curvature at the 
%   apex level that is going to be searched [cm-1](bmax_ini), the number of
%   curvature at the apex levels created for volume search (interval_b),
%   the minimum drop volume desired in uL (Vmin), the maximum drop volume
%   desired in uL (Vmax), the number of drop volumes levels desires (nVol),
%   the tolerance for volume determination (tVol), the inclination angle in
%   degree (alpha), the maximum number of iterations for volume search 
%   (itermax), the contact angle in degree (CA). The routine allows the 
%   user to choose the directory where the results are going to be 
%   exported. An excel file (.xlsx) is created containing important 
%   properties of the inclined drop for the different volumes. During the 
%   computation of inclined drop profiles errors may occur due to poor 
%   choices of input parameters. For certain combinations of curvature at 
%   the apex (bmin_ini and bmax_ini), inclination angle, contact angle and 
%   liquid capillary constant, the volume may not be found.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('-- INCLINED DROP - FIND CURVATURES AT THE APEX (b) THAT CORRESPONDS TO DESIRED DROP VOLUMES (V) --\n');
    resp = input("- Do you want to continue (Y/N) [Y]? ", "s");
    if isempty(resp) || (resp ~= 'Y' && resp ~= 'N')
        resp = 'Y';
    end
    if resp == 'N'
        break;
    end
    % Input parameters
    fprintf('Input parameters: \n');
    c = input("- Liquid capillary constant in cm-2 [e.g. water = 13.55]: "); %Capillary constant in cm-2
    if isempty(c)
        c = 13.55; %Capillary constant of water in cm-2
    end
    alpha_deg = input("- Inclination angle in degree [15]: ");
    if isempty(alpha_deg)
        alpha_deg = 15;
    end
    CA = input("- Contact angle in degree [70]: "); %Contact angle in degrees
    if isempty(CA)
        CA = 70;
    end
    bmin_ini = input("- Minimum curvature at the apex in cm-1 [4]: "); %Minimum curvature at the apex that is going to be searched in cm-1
    if isempty(bmin_ini)
        bmin_ini = 4; %For water
    end
    while 1
        bmax_ini = input("- Maximum curvature at the apex in cm-1 [9]: "); %Maximum curvature at the apex that is going to be searched in cm-1
        if isempty(bmax_ini)
            bmax_ini = 9; %For water
        end
        if bmax_ini < bmin_ini
            fprintf('Invalid maximum curavture at the apex. Please, set a lower value. \n');
        else
            break;
        end
    end
    interval_b = input("- Number of curvature at the apex levels created for volume search [50]: ");  %Number of b values created between bmin and bmax for volume search
    if isempty(interval_b)
        interval_b = 50;
    end
    Vmin = input("- Minimum drop volume in uL [2]: "); %Minimum drop volume in uL
    if isempty(Vmin)
        Vmin = 2;
    end
    while 1
        Vmax = input("- Maximum drop volume in uL [15]: "); %Maximum drop volume in uL
        if isempty(Vmax)
            Vmax = 15;
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
    tolV = input("- Tolerance for volume determination in uL [0.001]: ");  %Tolerance for volume determination in uL 
    if isempty(tolV)
        tolV = 1e-3;
    end
    itermax = input("- Maximum number of iterations for volume search [5]: "); %Maximum number of iterations for volume search
    if isempty(itermax)
        itermax = 5;
    end
    
    fprintf('Selection of the directory to export files ... \n');
    save_directory = uigetdir('C:\'); %Save directory
    %filename_xls = strcat('Dimensions_sessiledrop(c_',num2str(c),'_theta_',num2str(theta),').xlsx'); %Nome da planilha excel criada

    %- Numerical integration options
    N = 13; %number of elements along theta/ Suggestion: give N odd numbers, so the middle element is given at theta = 90°

    %VARIABLES INITIALIZATION
    %{
    theta_sim = zeros(1,interval_b); %Contact angle (theta) effectively used for drop simulation 
    wetted_radius_cm = zeros(1,interval_b);  %Wetting radius in cm
    r_eq_cm = zeros(1,interval_b); %Equatorial (maximum) radius of the drop in cm
    drop_height_cm = zeros(1,interval_b); %Drop height in cm
    A_cm2 = zeros(1,interval_b); %Drop surface area in cm²
    V_uL = zeros(1,interval_b); %Drop volume in uL calculated concomitantly with the Laplace curve generation
    V_otm = length(V);
    theta_otm = length(V);
    b_otm = length(V);
    r0_otm = length(V);
    wetted_radius_otm = length(V);
    r_eq_otm = length(V);
    drop_height_otm =  length(V);
    A_otm = length(V);
    %}
    CA_sim = zeros(1,interval_b); %Intrinsic apparent CA for a W distance from the apex of the drop in degrees
    CA_adv = zeros(1,interval_b); %Advancing CA in degrees
    CA_rec = zeros(1,interval_b); %Receding CA in degrees
    d_w_cm = zeros(1,interval_b); %Distance between contact points in cm
    h_cm = zeros(1,interval_b); %Drop height in cm
    A_cm2 = zeros(1,interval_b); %Drop surface area in cm²
    V_uL = zeros(1,interval_b); %Drop volume in mm³ or uL
    b_otm = length(V);
    r0_otm = length(V);
    CA_otm = length(V);
    CA_adv_otm = length(V);
    CA_rec_otm = length(V);
    d_w_otm = length(V);
    h_otm = length(V);
    A_otm = length(V);
    V_otm = length(V);

    % Find b values related to V
    for j = 1:length(V)
        diffV = 1; %reset volume difference
        bmin = bmin_ini; %reset bmin
        bmax = bmax_ini; %reset bmax
        iter = 1;
        while diffV > tolV && iter <= itermax
            fprintf('Searching for V = %.2f uL, between bmin = %.4f and bmax = %.4f (iteration %d)... \n',V(j),bmin,bmax,iter);
            b = linspace(bmin,bmax,interval_b);
            r0 = 1./b;
            for i=1:length(b)
                %Determination of the inclined drop profile by using numerical solution technique for solving Young-Laplace PDE
                [x_ly_sim_cm,z_ly_sim_cm,W,U,theta0,CA_sim(i),CA_adv(i),CA_rec(i),d_w_cm(i),h_cm(i),A_cm2(i),V_uL(i)] = profile_inclined_drop_ly(b(i),c,alpha_deg,CA,N);
            end
            diffV = min(abs(V_uL-V(j)));
            fprintf('diffV = %f\n',diffV);
            indx_V = find(abs(V_uL-V(j)) == diffV,1,'last');
            if (indx_V-1) <= 0
                bmin = b(indx_V);
                bmax = b(indx_V+1);
            elseif (indx_V+1) >= length(b)
                bmin = b(indx_V-1);
                bmax = b(indx_V);
            else
                bmin = b(indx_V-1);
                bmax = b(indx_V+1);
            end
            iter = iter +1;
        end
        V_otm(j) = V_uL(indx_V);
        b_otm(j) = b(indx_V);
        r0_otm(j) = r0(indx_V);
        CA_otm(j) = CA_sim(indx_V);
        CA_adv_otm(j) = CA_adv(indx_V);
        CA_rec_otm(j) = CA_rec(indx_V);
        d_w_otm(j) = d_w_cm(indx_V);
        h_otm(j) = h_cm(indx_V);
        A_otm(j) = A_cm2(indx_V);
    end
    %EXPORT RESULTS
    fprintf('Exporting results... \n');
    %- Write excel file with important drop properties
    varnames = {'b[cm-1]','r0[cm]','CA[°]','CA_max[°]','CA_min[°]','d_w[cm]','h[cm]','A[cm2]','V[uL]'}; %Variables names in excel file
    T = table(b_otm',r0_otm',CA_otm',CA_adv_otm',CA_rec_otm',d_w_otm',h_otm',A_otm',V_otm','VariableNames',varnames); %Build a table with important drop properties
    sheet_name = strcat('c = ',num2str(c),' cm^-2',' alpha = ',num2str(alpha_deg),' °'); %Excel sheet name
    savename_xls = strcat('Dimensions_Inclined_drop(c(cm-2)_',num2str(c),'_CA(°)_',num2str(mean(CA_sim),'%.0f'),'_alpha(°)_',num2str(alpha_deg),').xlsx'); %Excel file name
    fullsavename_xls = fullfile(save_directory,savename_xls);
    writetable(T,fullsavename_xls,'Sheet',sheet_name) %Write table to an excel file
    
    fprintf('Routine completed! \n');
    pause(2)
end
