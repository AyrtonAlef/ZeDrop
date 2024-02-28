%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Software developed for the evaluation of intermittent and continuous drop
% quasi-static experiments. Pendant, sessile and inclined drop 
% configurations may be evaluated. For pendant drops, the software analyses
% .gam results from Giada software and organize them in tables that are
% exported to Excel files allowing a detailed analysis of drop liquid 
% surface tension. For sessile and inclined drops, the software analyses 
% directly the capture images and determines important angles and drop 
% properties like contact angle, position, displacement and velocity of 
% contact points, and drop volume. Various pre-processing, edge detection, 
% post-processing, baseline determination and fitting options are
% available. During pre-processing of images, the sofwtare enables the user
% to enhance/correct contrast and to use smoothing filters. For edge
% detection, the software offers the use of gradient and non-gradient
% methods with pixel and sub-pixel precision. Morphological operations may
% be used during post-processing of images. For determining CA, different
% fitting methods may be used like circle fitting, polynomial fitting, 
% ellipse fitting and mask method. For baseline determination, the software
% offers options from the manual selection of baseline coordinates to the
% automatic identification of contact points for reflective and
% non-reflective drops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

% Main
while 1
    % Menu
    fprintf('--------------------------------------- ZEDROP_EVAL SOFTWARE -----------------------------------------\n');
    fprintf('MAIN MENU: \n');
    fprintf('1. Evaluation of intermittent PENDANT DROP quasi-static experiments (results from Giada); \n');
    fprintf('2. Evaluation of continuous PENDANT DROP quasi-static experiments (results from Giada); \n');
    fprintf('3. Evaluation of intermittent (needle out) SESSILE DROP (advancing and receding) quasi-static experiments; \n');
    fprintf('4. Evaluation of continuous (needle in) SESSILE DROP (advancing and receding) quasi-static experiments; \n');
    fprintf('5. Evaluation of continuous INCLINED DROP quasi-static experiments; \n');

    option = input('Enter one of the options above (0 to exit) [0]: ');
    
    if option == 0
        break;
    elseif option == 1 %Pendant drop (results from Giada)(intermittent)
        clc
        EvalIntPenDrop();
    elseif option == 2 %Pendant drop (results from Giada)(continuous)
        clc
        EvalContPenDrop();
    elseif option == 3 %Sessile drop (intermittent)
        clc
        EvalIntSesDropNeedleOut()
    elseif option == 4 %Sessile drop (continuous)
        clc
        EvalContSesDropNeedIn()
    elseif option == 5 %Inclined drop (continuous)
        clc
        EvalContIncDrop();
    end

    clc
    close all
end