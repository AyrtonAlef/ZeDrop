function [BaselineDetOpt,yDistWeight,marginBaselineDet,yAverageLR,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize] = SettingsBaselineDet(BaselineDetOpt,yDistWeight,marginBaselineDet,yAverageLR,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize)
%SETTINGSBASELINEDET Enable the user to edit baseline determination settings
%   INPUT/OUTPUT:
% BaselineDetOpt - Baselie determination option
% (0 - Manual/
% 1 - Automatic for reflective drop (Andersen and Taboryski 2017)/ 
% 2 - Automatic for non-reflective drop (Liu et al 2019)
% 3- Automatic for reflective and non-reflective drop (Akbari and Antonini 2021))
% yDistWeight = 40; %Weigth given to the vertical distance between edge points and drop baseline center (BaselineDetOpt == 0)
% marginBaselineDet = 10; %Margin from baseline selected for the search of contact points (BaselineDetOpt == 0)
% yAverageLR = 1; %Average y coordinate between left and right baseline selections (BaselineDetOpt == 0)
% nBaselineDet - Integer number of points that are fitted. Option for Andersen and Taboryski 2017 code (BaselineDetOpt == 1)
% cutBaselineDet - Number of pixels from the top where the script should not look for reflections. Option for Andersen and Taboryski 2017 code (BaselineDetOpt == 1)
% percentBaselineDet - Percent used to calculate cutoff value. Option for Liu et al 2019 code (BaselineDetOpt == 2)
% fillBlankHorLines = 1; %Option for filling blank horizontal lines (1 = enable and 0 = disable). Option for Akbari and Antonini 2021 code (BaselineDetOpt == 3).
% marginMaskBaselineDet = 10; %Minimum pixel distance (in y coordinate) between points along the detected edge to be considered part of the drop profile (Suggestion = jum_dist). Option for Akbari and Antonini 2021 code (BaselineDetOpt == 3). 
% maskSize = 20; %Mask size (neighbour index)

%Baseline determination settings
%- Current settings
fprintf('-------------------------------------------------------------- \n');
fprintf('- Current baseline determination settings: \n');
fprintf(strcat("BaselineDetOpt = ",num2str(BaselineDetOpt)," \n"));
fprintf(strcat("yDistWeight = ",num2str(yDistWeight)," \n"));
fprintf(strcat("marginBaselineDet = ",num2str(marginBaselineDet)," \n"));
fprintf(strcat("yAverageLR = ",num2str(yAverageLR)," \n"));
fprintf(strcat("nBaselineDet = ",num2str(nBaselineDet)," \n"));
fprintf(strcat("cutBaselineDet = ",num2str(cutBaselineDet)," \n"));
fprintf(strcat("percentBaselineDet = ",num2str(percentBaselineDet)," \n"));
fprintf(strcat("fillBlankHorLines = ",num2str(fillBlankHorLines)," \n"));
fprintf(strcat("marginMaskBaselineDet = ",num2str(marginMaskBaselineDet)," \n"));
fprintf(strcat("maskSize = ",num2str(maskSize)," \n"));
fprintf('-------------------------------------------------------------- \n');

%- Edit settings
fprintf('- Edit baseline determination settings... \n');
fprintf('Baseline determination methods available: \n');
fprintf('0. Manual baseline determination; \n');
fprintf('1. Automatic baseline detection for reflective drops; \n'); %from Andersen and Taboryski 2017 code
fprintf('2. Automatic baseline detection for non-reflective drops; \n'); %from Liu et. al 2019 code
fprintf('3. Automatic baseline detection for reflective and non-relfective drops; \n'); %from Akbari and Antonini 2021 code
BaselineDetOpt = input('Select an option [3]: ');
if isempty(BaselineDetOpt)
    BaselineDetOpt = 3;
end
if BaselineDetOpt == 0 %Manual baseline determination
    fprintf('Manual baseline determination selected. \n');
    yDistWeight = input('Weigth given to vertical distances between edge points and baseline center [40]: ');
    if isempty(yDistWeight)
        yDistWeight = 40;
    end
    marginBaselineDet = input('Margin for the identification of contact points [10]: ');
    if isempty(marginBaselineDet)
        marginBaselineDet = 10;
    end
    yAverageLR = input('Average y coordinate between left and right baseline selections (Y/N) [Y]? ','s');
    if isempty(yAverageLR)
        yAverageLR = 'Y';
    end
    if yAverageLR == 'N'
        yAverageLR = 0;
    else
        yAverageLR = 1;
    end
elseif BaselineDetOpt == 1 %Andersen and Taboryski 2017 code
    fprintf('Automatic baseline detection for reflective drops selected. \n');
    fprintf('Automatic baseline detection for reflective drops settings... \n');
    nBaselineDet = input('Integer number of points that are fitted [10]: ');
    if isempty(nBaselineDet)
        nBaselineDet = 10;
    end
    cutBaselineDet = input('Number of pixels from the top where the script should not look for reflections [260]: ');
    if isempty(cutBaselineDet)
        cutBaselineDet = 260;
    end
elseif BaselineDetOpt == 2 %Liu et al 2019 code
    fprintf('Automatic baseline detection for non-reflective drops selected. \n');
    fprintf('Automatic baseline detection for non-reflective drops settings... \n');
    while 1
        percentBaselineDet = input('Percent used to calculate cutoff value (min. = 0/ max. = 1) [0.99]: ');
        if isempty(percentBaselineDet) || percentBaselineDet < 0
            percentBaselineDet = 0.99;
        end
        if percentBaselineDet > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
elseif BaselineDetOpt == 3 %Akbari and Antonini 2021 code
    fprintf('Automatic baseline detection for reflective and non-reflective drops selected. \n');
    ansfillBlankHorLines = input('Fill blank horizontal lines (Y/N) [Y]? ','s');
    if isempty(ansfillBlankHorLines)
        ansfillBlankHorLines = 'Y';
    end
    if ansfillBlankHorLines == 'N'
        fillBlankHorLines = 0;
    else
        fillBlankHorLines = 1;
    end
    marginMaskBaselineDet = input('Minimum pixel distance (y coordinate) along the detected edge to be considered part of the drop profile [10]: ');
    if isempty(marginMaskBaselineDet )
        marginMaskBaselineDet  = 10;
    end
    maskSize = input('Mask size (neighbour index) [20]: ');
    if isempty(maskSize)
        maskSize  = 20;
    end
end
end

