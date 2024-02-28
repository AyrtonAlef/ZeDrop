function checkDimPendDrop_b_c_rh()
%CHECKDIMPENDDROP Check pendant drop profile and dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check pendant drop profile and dimensions. It requires the user to set:
%   the liquid capillary constant [cm-2](c), the minimum curvature at the 
%   apex [cm-1](bmin), the maximum curvature at the apex [cm-1](bmax), the 
%   number of profiles that is going to be checked between bmin and bmax
%   (nProfiles) and the needle diameter [mm](needleDiam). In addition to
%   the simulation of the pendant drop profiles for different values of 
%   curvature at the apex, the routine offers additional options such as: 
%   (i) highlight of important profile points, like the maximum drop 
%   diameter (indexDeq), the possbile intersections between drop profile 
%   and the needle diameter (indexLast and indexSecond), and the drop 
%   diameter at a distant equivalent to the maximum droplet diameter from 
%   the apex of the drop (indexDs), and (ii) presentation of pendant drop 
%   profiles trimmed by the diameter of the needle. The routine allows the 
%   user to choose the directory where the results are going to be 
%   exported. An excel file (.xlsx) is created containing  important 
%   properties of the pendant drops. During the computation of pendant 
%   drop profiles errors may occur due to poor choices of input 
%   parameters. For certain combinations of curvature at the apex 
%   (bmin and bmax) and needle diameter, the pendant drop may not be 
%   simulated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    clc
    clear all
    close all
    fprintf('------------------ PENDANT DROP - CHECK DROP PROFILE AND DIMENSIONS --------------------\n');
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
            fprintf('Invalid maximum curavture at the apex. Please, set a lower value. \n');
        else
            break;
        end
    end
    nProfiles = input("- Number of simulated profiles [3]: "); %Number of simulated profiles between minimum and maximum curvature at the apex
    if isempty(nProfiles)
            nProfiles = 3;
    end
    b = linspace(bmin,bmax,nProfiles);
    r0 = 1./b; %radius at the apex in cm
    needleDiam = input("- Needle diameter in mm [0.9088]: "); %Needle diameter in mm
    if isempty(needleDiam)
        needleDiam = 0.9088;
    end
    r_h_cm = (needleDiam/2)/10; %Holder outer radius in cm
    tol = 10^-5; %Tolerance to find the intersection points of the drop with the holder radius
    highlightPoints = input("- Highlight relevant drop profile points (indexDeq, indexDs, indexLast and indexSecond)(Y/N) [Y]? ", "s");
    if isempty(highlightPoints) || (highlightPoints ~= 'Y' && highlightPoints ~= 'N')
        highlightPoints = 'Y';
    end
    trimDrop = input("- Trim drop profile (Y/N) [Y]? ", "s");
    if isempty(trimDrop) || (trimDrop ~= 'Y' && trimDrop ~= 'N')
        trimDrop = 'Y';
    end
    imFmt = input("- Image format [.png]: ", "s");
    if isempty(imFmt)
        imFmt = ".png";
    end
    fprintf('Selection of the directory to export resulting files ... \n');
    save_directory = uigetdir('C:\'); %Save directory

    % Variables calculation and initialization
    S_span=(0:.0001:8); % S_span is the step variable for ode45 solver
    initialCo = [0 1e-100 0 0 0]; %Initial conditions for integration

    d_eq_cm = zeros(length(b),length(c)); %Maximum drop diameter in cm
    h_eq_cm = zeros(length(b),length(c)); %Distance from the apex of the drop to its maximum diameter in cm
    h_cm = zeros(length(b),length(c)); %Distance from the apex of the drop to the needle tip in cm
    d_s_cm = zeros(length(b),length(c)); %Drop diameter from a distance equivalent to the maximum drop diameter from the apex of the drop in cm
    A_cm2 = zeros(length(b),length(c)); %Drop surface area in cm²
    V_uL = zeros(length(b),length(c)); %Drop volume in uL calculated concomitantly with the Laplace curve generation
    fig = []; %Cell array containing epndant drop profile figures

    % Simulation of pendant drop profile
    for j = 1:length(c)
        for i = 1:length(b)
            fprintf('Running c = %.2f cm-2/ b = %.3f cm-1\n',c(j),b(i));
            % Determination of pendant drop profile
            [x_ly_cm,z_ly_cm,x_ly_sim_cm,z_ly_sim_cm,indx_last,indx_second,indx_deq,indx_ds,d_eq_cm(i,j),h_eq_cm(i,j),h_cm(i,j),d_s_cm(i,j),A_cm2(i,j),V_uL(i,j)] = profile_pendant_drop_ly(b(i),c(j),r_h_cm,S_span,initialCo,tol);

            % Plot pendant drop profile
            % - Plot pendant drop profile for different b values
            x_plot_lap_full = [wrev(-x_ly_cm); x_ly_cm];
            z_plot_lap_full = [wrev(z_ly_cm);z_ly_cm];
            displayname = strcat('c = ',num2str(c(j),'%.2f'),' cm-2/ ','b = ',num2str(b(i),'%.3f'), 'cm-1');
            figure(1)
            plot(x_plot_lap_full,z_plot_lap_full,'DisplayName',displayname);
            hold on
            % - Plot pendant drop profile highlighting important points
            if highlightPoints == 'Y'
                figure(2)
                plot(x_plot_lap_full,z_plot_lap_full,'DisplayName',displayname);
                hold on
                plot(x_ly_cm(indx_last),z_ly_cm(indx_last),'or','DisplayName','indxLast');
                hold on
                plot(x_ly_cm(indx_second),z_ly_cm(indx_second),'ob','DisplayName','indxSecond');
                hold on
                plot(x_ly_cm(indx_deq),z_ly_cm(indx_deq),'og','DisplayName','indxDeq');
                hold on
                if d_s_cm(i,j)~= 0
                    plot(x_ly_cm(indx_ds),z_ly_cm(indx_ds),'om','DisplayName','indxDs');
                    hold on
                end
                xlabel('[cm]');
                ylabel('[cm]');
            end
            % - Plot trimmed pendant drop profile according to holder radius
            if trimDrop == 'Y'
                x_plot_lap_sim_full = [wrev(-x_ly_sim_cm); x_ly_sim_cm];
                z_plot_lap_sim_full = [wrev(z_ly_sim_cm);z_ly_sim_cm];
                figure(3)
                plot(x_plot_lap_sim_full,z_plot_lap_sim_full,'DisplayName',displayname);
                hold on
                %plot(x_ly_cm(indx_last),z_ly_cm(indx_last),'or','DisplayName','indxLast');
                %hold on
                plot(x_ly_sim_cm(indx_second),z_ly_sim_cm(indx_second),'ob','DisplayName','indxSecond');
                hold on
                plot(x_ly_sim_cm(indx_deq),z_ly_sim_cm(indx_deq),'og','DisplayName','indxDeq');
                hold on
                if d_s_cm(i,j)~= 0
                    plot(x_ly_sim_cm(indx_ds),z_ly_sim_cm(indx_ds),'om','DisplayName','indxDs');
                end
                if i == 1
                    zaxis_max = z_ly_sim_cm(end);
                end
            end
        end
        % - Plot legend and needle diameter lines
        figure(1);
        lgd = legend;
        lgd.Location = "northwest";
        lgd.FontSize = 6;
        axis equal
        hold on
        x_rh = ones(length(z_plot_lap_full),1)*r_h_cm;
        plot([-x_rh;x_rh],[z_plot_lap_full;z_plot_lap_full],':k','DisplayName','holder simulation');
        xlabel('[cm]');
        ylabel('[cm]');
        fig{1} = figure(1);
        if highlightPoints == 'Y'
            figure(2)
            lgd = legend;
            lgd.Location = "northwest";
            lgd.FontSize = 6;
            axis equal
            hold on
            plot([-x_rh;x_rh],[z_plot_lap_full;z_plot_lap_full],':k','DisplayName','holder simulation');
            xlabel('[cm]');
            ylabel('[cm]');
            fig{2} = figure(2);
        end
        if trimDrop == 'Y'
            figure(3)
            lgd = legend;
            lgd.Location = "northwest";
            lgd.FontSize = 6;
            axis equal
            hold on
            x_rh = ones(2,1)*r_h_cm;
            zneedle_plot = [0; zaxis_max];
            %plot([-x_rh;x_rh],[zneedle_plot;zneedle_plot],':k','DisplayName','holder simulation');
            plot(-x_rh,zneedle_plot,':k','DisplayName','holder simulation');
            plot(x_rh,zneedle_plot,':k','HandleVisibility','off');
            xlabel('[cm]');
            ylabel('[cm]');
            fig{3} = figure(3);
        end
        % Export results
        fprintf('Exporting result... \n');
        %- Save graphics b x V(uL)/r(cm)/h(cm)
        for k = 1:3
            filename = strcat('dim_pendantdrop_',num2str(k),'_c(cm-2)_',num2str(c),'_rh(cm)_',num2str(r_h_cm,'%.4f'),'_bmin(cm-1)_',num2str(bmin,'%.2f'),'_bmax(cm-1)_',num2str(bmax,'%.2f'));
            filename_imFmt = strcat(filename,imFmt);
            filename_fig = strcat(filename,'.fig');
            fullfilename_imFmt = fullfile(save_directory,filename_imFmt);
            fullfilename_fig = fullfile(save_directory,filename_fig);
            saveas(fig{k},fullfilename_imFmt)
            saveas(fig{k},fullfilename_fig)
        end
        % Write excel file with important pendant drop profile properties
        filename_xls = strcat('Drop_dimensions(c_',num2str(c(j),'%.2f'),').xlsx'); %Excel file name
        fullfilename_xls = fullfile(save_directory,filename_xls);
        varnames = {'b[cm-1]','r0[cm]','De[cm]','h_De[cm]','h[cm]','Ds[cm]','A[cm2]','V[uL]'}; %Variable names
        T = table(b',r0',d_eq_cm,h_eq_cm,h_cm,d_s_cm,A_cm2,V_uL,'VariableNames',varnames);
        sheet_name = strcat('c = ',num2str(c(j),'%.2f'),' cm-2',' r_h = ',num2str(r_h_cm,'%.4f'),' cm'); %Excel sheet name
        writetable(T,fullfilename_xls,'Sheet',sheet_name) %Write table to an excel file
        fprintf('Routine completed! \n');
        pause(2)
    end
end
