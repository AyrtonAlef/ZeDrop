%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Software developed for checking and simulating different drop 
% configurations like pendant drop, sessile drop and inclined drop.
% Two-dimensional (2D) and three-dimensional (3D) drop profiles may be 
% generated. Binary and grayscale images may be simulated. The sofwtare
% also makes possible the creation of a set of consecutive images
% simulating quasi-static experiments. Some source of errors, such as edge
% width, camera vertical misalignment, random perturbation in the drop 
% profile, lack of contrast, non-uniform illumination, random noise, 
% Gaussian noise, impulsive salt-pepper noise and satellite droplets, may 
% also be simulated. Graphs and excel files are created containing important
% drop properties.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

% Main
while 1
    % Menu
    fprintf('--------------------------------------- ZEDROP_SIM SOFTWARE -----------------------------------------\n');
    fprintf('Definitions: \n');
    fprintf('c - liquid capillary constant \n');
    fprintf('b - curvature at the apex \n');
    fprintf('dh - holder (needle) diameter \n');
    fprintf('V - drop volume \n');
    fprintf('CA - contact angle \n');
    fprintf('alpha - inclination angle \n\n');

    fprintf('MAIN MENU: \n');
    fprintf('1. Pendant drop - check drop profile and dimensions; \n');
    fprintf('2. Pendant drop - simulation of binary images for different b - given c and dh; \n');
    fprintf('3. Pendant drop - simulation of grayscale images for different b - given c and dh; \n');
    fprintf('4. Pendant drop - simulation of quasi-static experiment; \n');
    fprintf('5. Sessile spherical cap drop - simulation of binary images for different V - given CA; \n');
    fprintf('6. Sessile spherical cap drop - simulation of grayscale images for different V - given CA; \n');
    fprintf('7. Sessile drop - check drop dimensions for different b and CA - given c; \n');
    fprintf('8. Sessile drop - create excel file containing b, V and CA; \n');
    fprintf('9. Sessile drop - simulation of binary images for different V - given c and CA; \n');
    fprintf('10. Sessile drop - simulation of grayscale images for different V - given c and CA; \n');
    %fprintf('11. Sessile drop - simulation of advancing/receding quasi-static experiments; \n');
    fprintf('11. Inclined drop - check drop dimensions for different alpha and b - given c and CA; \n');
    fprintf('12. Inclined drop - create excel file containing b, V, CA and alpha; \n');
    fprintf('13. Inclined drop - simulation of binary images for different V - given c, alpha and CA; \n');
    fprintf('14. Inclined drop - simulation of grayscale images for different V - given c, alpha and CA; \n');
    fprintf('15. Inclined drop - simulation of quasi-static experiment; \n');

    option = input('Enter one of the options above (0 to exit) [0]: ');
    
    if option == 0
        break;
    elseif option == 1 %Pendant drop - check drop profile and dimensions
        clc
        checkDimPendDrop_b_c_rh(); %Routine for check pendant drop profile and dimensions
    elseif option == 2 %Pendant drop - simulation of binary images
        clc
        simBinPendDrop_b_c_rh();
    elseif option == 3 %Pendant drop - simulation of grayscale images
        clc
        simGrayPendDrop_b_c_rh();
    elseif option == 4 %Pendant drop - simulation of quasi-static experiment
        clc
        simExpPendDrop();
     elseif option == 5 %Spherical cap - simulation of binary images
        clc
        simBinSphCapDrop_V_theta_rh();
    elseif option == 6 %Spherical cap - simulation of grayscale images
        clc
        simGraySphCapDrop_V_theta_rh();
    elseif option == 7 %Sessile drop - check drop profile and dimensions
        clc
        checkDimSessDrop_b_c_CA();
    elseif option == 8 %Sessile drop - create excel file containg b, V and CA
        clc
        checkSessDropFindbRelatedToV;
    elseif option == 9 %Sessile drop - simulation of binary images
        clc
        simBinSessDrop_V_c_CA();
    elseif option == 10 %Sessile drop - simulation of grayscale images
        clc
        simGraySessDrop_V_c_CA();
    %elseif option == 11 %Sessile drop - simulation of advancing/receding quasi-static experiments
        %clc
        %simExpAdvRecSessDrop();
    elseif option == 11 %Inclined drop - check drop dimensions for different alpha and b (given c and CA)
        clc
        checkDimIncDrop_b_c_CA_alpha();
    elseif option == 12 %Inclined drop - create excel file containg b, V, CA and alpha
        clc
        checkIncDropFindbRelatedToV;
    elseif option == 13 %Inclined drop - simulation of binary images for different V - given c, alpha and CA
        clc
        simBinIncDrop_V_c_CA();
    elseif option == 14 %Inclined drop - simulation of grayscale images for different V - given c, alpha and CA
        clc
        simGrayIncDrop_V_c_CA();
    elseif option == 15 %Inclined drop - simulation of quasi-static drop experiment
        clc
        simExpIncDrop();
    end

    clc
    close all
end