function checkDimSessDrop_b_c_CA()
%CHECKDIMSESSDROP Check sessile drop dimensions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check sessile drop dimensions. It requires the user to set: the liquid 
%   capillary constant [cm-2](c), the minimum curvature at the apex [cm-1]
%   (bmin), the maximum curvature at the apex [cm-1](bmax), the number of 
%   curvatures at the apex that is going to be checked between bmin and 
%   bmax (nb), the minimum contact angle [°](thetamin), the maximum contact
%   angle [°](thetamax), the number of contact angle levels that is going 
%   to be simulated between thetamin and thetamax (ntheta). The routine 
%   allows the user to choose the directory where the results are going to 
%   be exported. Figures and an excel file (.xlsx) are created containing 
%   important properties of the sessile drop for different curvatures at 
%   the apex. During the computation of sessile drop profiles errors may 
%   occur due to poor choices of input parameters. For certain combinations
%   of curvature at the apex (bmin and bmax) the sessile drop may not be 
%   simulated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ SESSILE DROP - CHECK DROP DIMENSIONS --------------------\n');
    resp = input("- Do you want to continue (Y/N) [Y]? ", "s");
    if isempty(resp) || (resp ~= 'Y' && resp ~= 'N')
        resp = 'Y';
    end
    if resp == 'N'
        break;
    end
    % Input parameters
    fprintf('Input parameters: \n');
    c = input("- Liquid capillary constant in cm-2 [e.g. water = 13.55]: "); %Cappilary constant in cm-2
    if isempty(c)
        c = 13.55; %Capillary constant of water in cm-2
    end
    bmin = input("- Minimum curvature at the apex in cm-1 [0.05]: "); %Minimum curvature at the apex in cm-1
    if isempty(bmin)
        bmin = 0.05; %For water
     end
    while 1
        bmax = input("- Maximum curvature at the apex in cm-1 [10]: "); %Maximum curvature at the apex in cm-1;
        if isempty(bmax)
            bmax = 10; %For water
        end
        if bmax < bmin
            fprintf('Invalid maximum curavture at the apex. Please, set a lower value. \n');
        else
            break;
        end
    end
    nb = input("- Number of curvature at the apex levels simulated [100]: "); %Number of simulated curvature at the apex between bmin and bmax
    if isempty(nb)
            nb = 100;
    end
    b = linspace(bmin,bmax,nb);
    r0 = 1./b; %radius at the apex in cm
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
    imFmt = input("- Image format [.png]: ", "s");
    if isempty(imFmt)
        imFmt = ".png";
    end
    fprintf('Selection of the directory to export files ... \n');
    save_directory = uigetdir('C:\'); %Save directory
    filename_xls = strcat('Dimensiosn_sessiledrop(c_',num2str(c),'_bmin_',num2str(bmin),'_bmax_',num2str(bmax),'_thetamin_',num2str(thetamin),'_thetamax_',num2str(thetamax),').xlsx'); %Excel file name

    % Variables calculation and initialization
    S_span = (0:.0001:8); % S_span is the step variable for ode45 solver
    initialCo = [0 1e-100 0 0 0]; %Initial conditions for integration

    %VARIABLES INITIALIZATION
    theta_sim = zeros(length(b),length(theta)); %Contact angle (theta) effectively used for drop simulation 
    wetted_radius_cm = zeros(length(b),length(theta)); %Wetting radius in cm
    r_eq_cm = zeros(length(b),length(theta)); %Equatorial (maximum) radius of the drop in cm
    drop_height_cm = zeros(length(b),length(theta)); %Drop height in cm
    A_cm2 = zeros(length(b),length(theta)); %Drop surface area in cm²
    V_uL = zeros(length(b),length(theta)); %Drop volume in uL calculated concomitantly with the Laplace curve generation

    % Simulation of sessile drop profile
    for j=1:length(theta)
        for i=1:length(b)
            fprintf('Running theta = %.2f°/ b = %.2f cm-1\n',theta(j),b(i));
            %Determination of the sessile drop profile by numerical integration of the Young-Laplace equation
            [x_ly_sim_cm,z_ly_sim_cm,theta_sim(i,j),wetted_radius_cm(i,j),r_eq_cm(i,j),drop_height_cm(i,j),A_cm2(i,j),V_uL(i,j)] = profile_sessile_drop_ly(b(i),c,theta(j),S_span,initialCo);
        end
        %Plot of drop dimensions for differente b values
        figure
        yyaxis left
        plot(b',V_uL(:,j),'-o','DisplayName','V (uL)')
        graphname = strcat('Drop dimensions for c = ',num2str(c),' cm-2/ theta = ',num2str(theta(j),'%.2f'), ' °');
        title(graphname)
        xlabel('b (cm-1)')
        ylabel('V (uL)')
        yyaxis right
        plot(b',wetted_radius_cm(:,j),'--*','DisplayName','r (cm)')
        hold on
        ylabel('r or h (cm)')
        yyaxis right
        plot(b',drop_height_cm(:,j),':x','DisplayName','h (cm)')
        legend
        
        %EXPORT RESULTS
        fprintf('Exporting results... \n');
        %- Save graphics b x V(uL)/r(cm)/h(cm)
        filename_imFmt = strcat('dim_sessiledrop_c_',num2str(c),'_theta_',num2str(theta(j),'%.2f'));
        filename_fig = strcat('dim_sessiledrop_c_',num2str(c),'_theta_',num2str(theta(j),'%.2f'),'.fig');
        fullfilename = strcat(filename_imFmt,imFmt); %Image name + image format
        fullfilename_imFmt = fullfile(save_directory,fullfilename); %Image fullfilename
        fullfilename_fig = fullfile(save_directory,filename_fig);
        saveas(gcf,fullfilename_imFmt)
        saveas(gcf,fullfilename_fig)
        %- Write excel file with important drop properties
        varnames = {'b[cm-1]','r0[cm]','theta[°]','r_w[cm]','r_eq[cm]','h[cm]','V[uL]'}; %Excel file name
        T = table(b',r0',theta_sim(:,j),wetted_radius_cm(:,j),r_eq_cm(:,j),drop_height_cm(:,j),V_uL(:,j),'VariableNames',varnames);
        sheet_name = strcat('c = ',num2str(c),' cm^-2',' theta = ',num2str(theta(j),'%.2f'),' °'); %Excel sheet name
        fullfilename_xls = fullfile(save_directory,filename_xls);
        writetable(T,fullfilename_xls,'Sheet',sheet_name) %Write table to an excel file
    end
    fprintf('Routine completed! \n');
    pause(2)
end
end
