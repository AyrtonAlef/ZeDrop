function [EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter] = SettingsPreprocessing(EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter)
%SETTINGSPREPROCESSING Enable the user to edit pre-processing image settings
%   INPUT/OUTPUT:
% EnhanceContrast - Option to enhance contrast
% low_in - Low contrast limit for input image
% high_in - High contrast limit for input image
% low_out - Low contrast limit for output image
% high_out - High contrast limit for output image
% SmoothingFilter - Option for applying Gaussian smoothing filter
% sigmaSmoothFilter - Standard deviation of the 2-D Gaussian smoothing kernel
% sizeSmoothFilter - Filter size

%Pre-processing settings
%Current settings
fprintf('-------------------------------------------------------------- \n');
fprintf('- Current pre-processing settings: \n');
fprintf(strcat("EnhanceContrast = ",num2str(EnhanceContrast)," \n"));
fprintf(strcat("low_in = ",num2str(low_in)," \n"));
fprintf(strcat("high_in = ",num2str(high_in)," \n"));
fprintf(strcat("low_out = ",num2str(low_out)," \n"));
fprintf(strcat("high_out = ",num2str(high_out)," \n"));
fprintf(strcat("SmoothingFilter = ",num2str(SmoothingFilter)," \n"));
fprintf(strcat("sigmaSmoothFilter = ",num2str(sigmaSmoothFilter)," \n"));
fprintf(strcat("sizeSmoothFilter = ",num2str(sizeSmoothFilter)," \n"));
fprintf('-------------------------------------------------------------- \n');
%Edit settings
fprintf('- Edit pre-processing image settings... \n');
%Contrast settings
fprintf('Enhance contrast settings... \n');
ansContrast = input('Do you like to enhance contrast (Y/N) [Y]? ','s');
if isempty(ansContrast)
    ansContrast = 'Y';
end
if ansContrast == 'Y'
    EnhanceContrast = 1; %Enable enhance contrast option
    while 1
        low_in = input('Low contrast limit for input image (min. = 0/ max. = 1) [0]: ');
        if isempty(low_in) || low_in < 0
            low_in = 0;
        end
        if low_in > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
    while 1
        high_in = input('High contrast limit for input image (min. = 0/ max. = 1) [1]: ');
        if isempty(high_in) || high_in < 0
            high_in = 1;
        end
        if high_in < low_in
            fprintf('High contrast limit lower than the low contrast limit. Please, enter a new value. \n')
        elseif high_in > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
    while 1
        low_out = input('Low contrast limit for output image (min. = 0/ max. = 1) [0]: ');
        if isempty(low_out) || low_out < 0
            low_out = 0;
        end
        if low_out > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
    while 1
        high_out = input('High contrast limit for output image (min. = 0/ max. = 1) [1]: ');
        if isempty(high_out) || high_out < 0
            high_out = 1;
        end
        if high_out < low_out
            fprintf('High contrast limit lower than the low contrast limit. Please, enter a new value. \n')
        elseif high_out > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
else
    EnhanceContrast = 0; %Disable enhance contrast option
end
%Smoothing settings
fprintf('Smoothing filter settings... \n');
ansSmoothing = input('Do you like to enable smoothing filter (Gaussian filter) (Y/N) [Y]? ','s');
if isempty(ansSmoothing)
    ansSmoothing = 'Y';
end
if ansSmoothing == 'Y'
    SmoothingFilter = 1;
    sigmaSmoothFilter = input('Standard deviation of Gaussian distribution [0.5]: ');
    if isempty(sigmaSmoothFilter)
        sigmaSmoothFilter = 0.5;
    end
    defaultFilterSize = 2*ceil(2*sigmaSmoothFilter)+1; %Default value for filter size;
    strSmoothFilter = strcat('Filter size [',num2str(defaultFilterSize),']: ');
    sizeSmoothFilter = input(strSmoothFilter);
    if isempty(sizeSmoothFilter)
        sizeSmoothFilter = defaultFilterSize;
    end
else
    SmoothingFilter = 0;
end
end

