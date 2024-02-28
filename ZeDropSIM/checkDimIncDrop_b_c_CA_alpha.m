function checkDimIncDrop_b_c_CA_alpha()
%CHECKDIMSINCDROP_B_C_CA_ALPHA Check inclined drop dimensions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check inclined drop dimensions. It requires the user to enter: the 
%   liquid capillary constant [cm-2](c), the apparent contact angle
%   [°](theta), the minimum inclination angle [°](alphamin), the maximum
%   inclination angle [°](alphamax), the number of inclination angle levels
%   that is going to be simulated between alphamin and alphamax (nalpha),
%   the minimum curvature at the apex for each inclination angle
%   [cm-1](bmin), the maximum curvature at the apex for each inclination
%   angle [cm-1](bmax) and the number of levels of curvature at the apex 
%   that is going to be simulated between bmin and bmax (nb). The routine 
%   allows the user to choose the directory where the results are going to 
%   be exported. Graphs and an excel file (.xlsx) are created containing 
%   important properties of the inclined drop for different curvatures at 
%   the apex. During the computation of inclined drop profiles errors may 
%   occur due to poor choices of input parameters. For certain combinations
%   of curvature at the apex (bmin and bmax) the inclined drop may not be 
%   simulated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ INCLINED DROP - CHECK DROP DIMENSIONS --------------------\n');
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
    theta = input("- Apparent contact angle in degree [70]: "); %Minimum contact angle in degree
    if isempty(theta)
        theta = 70;
    end
    alphamin = input("- Minimum inclination angle in degree [15]: "); %Minimum inclination angle in degree
    if isempty(alphamin)
        alphamin = 15;
    end
    while 1
        alphamax = input("- Maximum inclination angle in degree [60]: "); %Maximum inclination angle in degree
        if isempty(alphamax)
            alphamax = 60;
        end
        if alphamax < alphamin
            fprintf('Invalid maximum inclination angle. Please, set a lower value. \n');
        else
            break;
        end
    end
    nalpha = input("- Number of inclination angle levels simulated [3]: "); %Number of simulated inclination angles between alphamin and alphamax
    if isempty(nalpha)
        nalpha = 3;
    end
    alpha_deg = linspace(alphamin,alphamax,nalpha); %inclination angle in degrees
    %bmin = zeros(1,nalpha);
    %bmax = zeros(1,nalpha);
    nb = input("- Number of curvature at the apex levels simulated [100]: "); %Number of simulated curvature at the apex between bmin and bmax
    if isempty(nb)
            nb = 100;
    end
    b = zeros(nb,nalpha);
    for i = 1:nalpha
        while 1
            bminStr = strcat("- Minimum curvature at the apex in cm-1 for an inclination angle of ",num2str(alpha_deg(i),'%.2f'),"°: ");
            bmin = input(bminStr); %Minimum curvature at the apex in cm-1
            if ~isempty(bmin)
                break;
            end
        end
        while 1
            bmaxStr = strcat("- Maximum curvature at the apex in cm-1 for an inclination angle of ",num2str(alpha_deg(i),'%.2f'),"°: ");
            bmax = input(bmaxStr); %Maximum curvature at the apex in cm-1;
            if isempty(bmax)

            elseif bmax < bmin
                fprintf('Invalid maximum curavture at the apex. Please, set a lower value. \n');
            else
                break;
            end
        end
        b(:,i) = linspace(bmin,bmax,nb);
    end
    r0 = 1./b; %radius at the apex in cm  
    imFmt = input("- Image format [.png]: ", "s");
    if isempty(imFmt)
        imFmt = ".png";
    end
    %- Numerical integration options
    N = 13; %number of elements along theta/ Suggestion: give N odd numbers, so the middle element is given at theta = 90°

    %- Save options
    fprintf('Selection of the directory to export files ... \n');
    save_directory = uigetdir('C:\'); %Save directory
    filename_xls = strcat('Dimensions_Inclined_drop(c(cm-2)_',num2str(c),'_CA(°)_',num2str(theta),'_alphamin(°)_',num2str(alphamin,'%.2f'),'_alphamax(°)_',num2str(alphamax,'%.2f'),').xlsx'); %Excel file name
    graphname1 = strcat("Drop dimensions for c = ", num2str(c)," cm-2/ CA = ",num2str(theta),"°");
    graphname2 = strcat("Drop contact angles for c = ", num2str(c)," cm-2/ CA = ",num2str(theta),"°");
    all_marks = {'o','s','+','*','.','x','d','^','v','>','<','p','h'};
    all_colors ={'r','g','b','c','m','y','k'};

    %VARIABLES INITIALIZATION 
    CA_sim = zeros(nb,nalpha); %Intrinsic apparent CA for a W distance from the apex of the drop in degrees
    CA_adv = zeros(nb,nalpha); %Advancing CA in degrees
    CA_rec = zeros(nb,nalpha); %Receding CA in degrees
    d_w_cm = zeros(nb,nalpha); %Distance between contact points in cm
    h_cm = zeros(nb,nalpha); %Drop height in cm
    A_cm2 = zeros(nb,nalpha); %Drop surface area in cm²
    V_uL = zeros(nb,nalpha); %Drop volume in mm³ or uL
    
    % Simulation of inclined drop profile
    for j = 1:nalpha
        %fprintf('Running for alpha = %.2f°, bmin = %.3f cm-1 and bmax = %.3f cm-1 ... \n',alpha_deg(j),b(1,j),b(end,j));
        fprintf('-------------------------------------------------------------------------- \n');
        for i = 1:nb
            fprintf('Running for alpha = %.2f°, b = %.3f cm-1 ...\n',alpha_deg(j),b(i,j));
            %fprintf('b = %.3f cm-1 ... \n',b(i,j));
            %Determination of the inclined drop profile by using numerical solution technique for solving Young-Laplace PDE
            [x_ly_sim_cm,z_ly_sim_cm,W,U,theta0,CA_sim(i,j),CA_adv(i,j),CA_rec(i,j),d_w_cm(i,j),h_cm(i,j),A_cm2(i,j),V_uL(i,j)] = profile_inclined_drop_ly(b(i,j),c,alpha_deg(j),theta,N);
        end
        %Plot drop volume, weting diameter and drop height for different values of b and export results
        fig1 = figure(1);
        yyaxis left
        dispNameV = strcat("V(uL)/ alpha = ",num2str(alpha_deg(j),'%.2f'),"°");
        plot(b(:,j),V_uL(:,j),'LineStyle','--','Marker','o','Color',all_colors{mod(j,7)},'DisplayName',dispNameV)
        title(graphname1)
        xlabel('b (cm-1)')
        ylabel('V (uL)')
        yyaxis right
        dispNamedw = strcat("dw(cm)/ alpha = ",num2str(alpha_deg(j),'%.2f'),"°");
        plot(b(:,j),d_w_cm(:,j),'LineStyle','--','Marker','s','Color',all_colors{mod(j,7)},'DisplayName',dispNamedw)
        hold on
        ylabel('dw or h (cm)')
        yyaxis right
        dispNameh = strcat("h(cm)/ alpha = ",num2str(alpha_deg(j),'%.2f'),"°");
        plot(b(:,j),h_cm(:,j),'LineStyle','--','Marker','^','Color',all_colors{mod(j,7)},'DisplayName',dispNameh)
        lgd1 = legend;
        lgd1.FontSize = 8;
        hold on

        %Plot CA for different values of b and export results
        fig2 = figure(2);
        %yyaxis left
        dispNameCAmax = strcat("CAmax(°)/ alpha = ",num2str(alpha_deg(j),'%.2f'),"°");
        plot(b(:,j),CA_adv(:,j),'LineStyle','--','Marker','o','Color',all_colors{mod(j,7)},'DisplayName',dispNameCAmax)
        title(graphname2)
        xlabel('b (cm-1)')
        ylabel('CA (°)')
        hold on
        dispNameCAmin = strcat("CAmin(°)/ alpha = ",num2str(alpha_deg(j),'%.2f'),"°");
        plot(b(:,j),CA_rec(:,j),'LineStyle','--','Marker','s','Color',all_colors{mod(j,7)},'DisplayName',dispNameCAmin)
        lgd2 = legend;
        lgd2.FontSize = 8;
        hold on

        %Write excel file with important inclined drop profile properties
        fprintf('Exporting excel file containing important inclined drop properties... \n');
        varnames = {'b[cm-1]','r0[cm]','CA[°]','CA_max[°]','CA_min[°]','d_w[cm]','h[cm]','A[cm2]','V[uL]'}; %Variable names in the excel file
        T = table(b(:,j),r0(:,j),CA_sim(:,j),CA_adv(:,j),CA_rec(:,j),d_w_cm(:,j),h_cm(:,j),A_cm2(:,j),V_uL(:,j),'VariableNames',varnames);
        sheet_name = strcat('c = ',num2str(c),' cm^-2',' alpha = ',num2str(alpha_deg(j),'%.2f'),' °'); %Excel sheet name
        fullsavename_xls = fullfile(save_directory,filename_xls);
        writetable(T,fullsavename_xls,'Sheet',sheet_name) %Write table to an excel file
    end
    %Saving graphic b(cm-1) x V(uL)
    fprintf('Saving graphics... \n');
    filenamebV = strcat("DropDimensionIncDrop_bxV_c(cm-2)_",num2str(c),"_CA_",num2str(theta,'%.2f'));
    filenamebV_imFmt = strcat(filenamebV,imFmt);
    filenamebV_fig = strcat(filenamebV,'.fig');
    fullfilenamebV_imFmt = fullfile(save_directory,filenamebV_imFmt);
    fullfilenamebV_fig = fullfile(save_directory,filenamebV_fig);
    saveas(fig1,fullfilenamebV_imFmt)
    saveas(fig1,fullfilenamebV_fig)
    filenamebCA = strcat("DropCAsIncDrop_bxCA_c(cm-2)_",num2str(c),"_CA_",num2str(theta,'%.2f'));
    filenamebCA_imFmt = strcat(filenamebCA,imFmt);
    filenamebCA_fig = strcat(filenamebCA,'.fig');
    fullfilenamebCA_imFmt = fullfile(save_directory,filenamebCA_imFmt);
    fullfilenamebCA_fig = fullfile(save_directory,filenamebCA_fig);
    saveas(fig2,fullfilenamebCA_imFmt)
    saveas(fig2,fullfilenamebCA_fig)
    %close all    
    fprintf('Routine completed! \n');
    pause(2)
end      
end
