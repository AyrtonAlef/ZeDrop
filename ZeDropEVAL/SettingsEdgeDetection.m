function [EdgeDetMethodOpt,edgeDetMetName,edgeDetDirection,edgeDetSigma,edgeDetThreshold,edgeDetFilter,edgeDetFilterTypeOpt,edgeDetFilterTypeName,edgeDetFilterSize,edgeDetFilterSigma,edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter] = SettingsEdgeDetection(EdgeDetMethodOpt,edgeDetMetName,edgeDetDirection,edgeDetSigma,edgeDetThreshold,edgeDetFilter,edgeDetFilterTypeOpt,edgeDetFilterTypeName,edgeDetFilterSize,edgeDetFilterSigma,edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter)
%SETTINGSEDGEDETECTION Enable the user to edit edge detection settings
%   INPUT/OUTPUT:
% EdgeDetMethodOpt - Edge detection method option (1. Sobel/ 2. Prewitt/ 3. Roberts/ 4. LoG/ 5. Canny/ 6. Zerocross/ 7. Partial area effect)
% edgeDetMetName - Name of the edge detection method
% edgeDetDirection - Direction of edges to detect (both/ horizontal/ vertical). Valid for Sobel, Prewitt and Roberts methods
% edgeDetSigma - Standard deviation of filter. Valid for LoG and Canny gradient methods
% edgeDetThreshold - Sensitivity threshold (min. = 0/ max. = 1). Valid for Canny gradient methods
% edgeDetFilter - Option to enable the definition of filter to zerocross method
% edgeDetFilterTypeOpt - Edge detection filter type option. Option for zerocross function
% edgeDetFilterTypeName - Type of filter. Option for zerocross function
% edgeDetFilterSize - Filter size. Option for zerocross function and filter types: average, disk, gaussian and laplacian
% edgeDetFilterSigma - Standard deviation. Option for zerocross function and filter type: gaussian
% edgeDetFilterShape - Shape of the Laplacian (min. = 0/ max. = 1). Option for zerocross function and filter type: laplacian
% edgeDetDifThreshold - Minimum difference of intensity at both sides of a pixel to be considered an edge
% edgeDetOrder - Order of the edge to find(1. First order edges (straight lines)/ 2. Second order edges)
% edgeDetSmoothingIter - Smoothing iterations needed to find final edges (0. Oriented to noise free images/ 1. Oriented to low-noise images/ >1. Oriented to high-noise images) 

%Edge detection method
    %- Current settings
    fprintf('-------------------------------------------------------------- \n');
    fprintf('- Current edge detection settings: \n');
    fprintf(strcat("EdgeDetMethodOpt = ",num2str(EdgeDetMethodOpt)," \n"));
    fprintf(strcat("edgeDetMetName = ",edgeDetMetName," \n"));
    fprintf(strcat("edgeDetDirection = ",edgeDetDirection," \n"));
    fprintf(strcat("edgeDetSigma = ",num2str(edgeDetSigma)," \n"));
    fprintf(strcat("edgeDetThreshold = ",num2str(edgeDetThreshold)," \n"));
    fprintf(strcat("edgeDetFilter = ",num2str(edgeDetFilter)," \n"));
    fprintf(strcat("edgeDetFilterTypeOpt = ",edgeDetFilterTypeOpt," \n"));
    fprintf(strcat("edgeDetFilterTypeName = ",edgeDetFilterTypeName," \n"));
    fprintf(strcat("edgeDetFilterSize = ",num2str(edgeDetFilterSize)," \n"));
    fprintf(strcat("edgeDetFilterSigma = ",num2str(edgeDetFilterSigma)," \n"));
    fprintf(strcat("edgeDetFilterShape = ",num2str(edgeDetFilterShape)," \n"));
    fprintf(strcat("edgeDetDifThreshold = ",num2str(edgeDetDifThreshold)," \n"));
    fprintf(strcat("edgeDetOrder = ",num2str(edgeDetOrder)," \n"));
    fprintf(strcat("edgeDetSmoothingIter = ",num2str(edgeDetSmoothingIter)," \n"));
    fprintf('-------------------------------------------------------------- \n');

    %- Edit settings
while 1
    fprintf('- Edit edge detection settings... \n');
    fprintf('Edge detection methods available: \n');
    fprintf('1. Sobel \n');
    fprintf('2. Prewitt \n');
    fprintf('3. Roberts \n');
    fprintf('4. LoG \n');
    fprintf('5. Canny \n');
    fprintf('6. Zerocross \n');
    fprintf('7. Partial area effect \n');
    EdgeDetMethodOpt = input('Select a method [7]: ');
    if isempty(EdgeDetMethodOpt)
        EdgeDetMethodOpt = 7;
    end
    if EdgeDetMethodOpt > 0 || EdgeDetMethodOpt <= 7
        break;
    end
end
if EdgeDetMethodOpt == 1 %Sobel method
    fprintf('Sobel edge detection method selected. \n');
    fprintf('Sobel edge detection method settings... \n');
    edgeDetMetName = 'Sobel';
    while 1
        edgeDetThreshold = input('Sensitivity threshold (min. = 0/ max. = 1) []: ');
        if isempty(edgeDetThreshold) || edgeDetThreshold < 0
            edgeDetThreshold = [];
        end
        if edgeDetThreshold > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
    edgeDetDirection = input('Direction of edges to detect (both, horizontal or vertical) [both]: ','s');
    if isempty(edgeDetDirection) || (~strcmp(edgeDetDirection,'horizontal') && ~strcmp(edgeDetDirection,'vertical'))
        edgeDetDirection = 'both';
    end
elseif EdgeDetMethodOpt == 2 %Prewitt method
    fprintf('Prewitt edge detection method selected. \n');
    fprintf('Prewitt edge detection method settings... \n');
    edgeDetMetName = 'Prewitt';
    while 1
        edgeDetThreshold = input('Sensitivity threshold (min. = 0/ max. = 1) []: ');
        if isempty(edgeDetThreshold) || edgeDetThreshold < 0
            edgeDetThreshold = [];
        end
        if edgeDetThreshold > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
    edgeDetDirection = input('Direction of edges to detect (both, horizontal or vertical) [both]: ','s');
    if isempty(edgeDetDirection) || (~strcmp(edgeDetDirection,'horizontal') && ~strcmp(edgeDetDirection,'vertical'))
        edgeDetDirection = 'both';
    end
elseif EdgeDetMethodOpt == 3 %Roberts method
    fprintf('Roberts edge detection method selected. \n');
    fprintf('Roberts edge detection method settings... \n');
    edgeDetMetName = 'Roberts';
    while 1
        edgeDetThreshold = input('Sensitivity threshold (min. = 0/ max. = 1) []: ');
        if isempty(edgeDetThreshold) || edgeDetThreshold < 0
            edgeDetThreshold = [];
        end
        if edgeDetThreshold > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
    edgeDetDirection = input('Direction of edges to detect (both, horizontal or vertical) [both]: ','s');
    if isempty(edgeDetDirection) || (~strcmp(edgeDetDirection,'horizontal') && ~strcmp(edgeDetDirection,'vertical'))
        edgeDetDirection = 'both';
    end
elseif EdgeDetMethodOpt == 4 %Laplace of Gaussian (LoG) method
    fprintf('Laplace of Gaussian (LoG) edge detection method selected. \n');
    fprintf('Laplace of Gaussian (LoG) edge detection method settings... \n');
    edgeDetMetName = 'log';
    while 1
        edgeDetThreshold = input('Sensitivity threshold (min. = 0/ max. = 1) []: ');
        if isempty(edgeDetThreshold) || edgeDetThreshold < 0
            edgeDetThreshold = [];
        end
        if edgeDetThreshold > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
    edgeDetSigma = input('Standard deviation of filter [2]: ');
    if isempty(edgeDetSigma)
        edgeDetSigma = 2;
    end
elseif EdgeDetMethodOpt == 5 %Canny method
    fprintf('Canny edge detection method selected. \n');
    fprintf('Canny edge detection method settings... \n');
    edgeDetMetName = 'Canny';
    while 1
        edgeDetThreshold = input('Sensitivity threshold (min. = 0/ max. = 1) []: ');
        if isempty(edgeDetThreshold) || edgeDetThreshold < 0
            edgeDetThreshold = [];
        end
        if edgeDetThreshold > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
    edgeDetSigma = input('Standard deviation of filter [2]: ');
    if isempty(edgeDetSigma)
        edgeDetSigma = 2;
    end
elseif EdgeDetMethodOpt == 6 %Zero-crossing method
    fprintf('Zero-crossing edge detection method selected. \n');
    fprintf('Zero-crossing detection method settings... \n');
    edgeDetMetName = 'zerocross';
    while 1
        edgeDetThreshold = input('Sensitivity threshold (min. = 0/ max. = 1) []: ');
        if isempty(edgeDetThreshold) || edgeDetThreshold < 0
            edgeDetThreshold = [];
        end
        if edgeDetThreshold > 1
            fprintf('Maximum limit reached. Please, enter a new value. \n')
        else
            break;
        end
    end
    ansGradFilter = input('Do you want to determine the filter of zero-crossing method (Y/N) [N]?','s');
    if isempty(ansGradFilter) || (~strcmp(ansGradFilter,'Y') && ~strcmp(ansGradFilter,'N'))
        ansGradFilter = 'N';
        edgeDetFilter = 0;
    end
    if ansGradFilter == 'Y'
        edgeDetFilter = 1;
        edgeDetFilterTypeOpt = input('Filter type: average (A), disk (D), Gaussian (G) or Laplacian (L) [G]?','s');
        if isempty(edgeDetFilterTypeOpt) || (~strcmp(edgeDetFilterTypeOpt,'A') && ~strcmp(edgeDetFilterTypeOpt,'D') && ~strcmp(edgeDetFilterTypeOpt,'L'))
            edgeDetFilterTypeOpt = 'G';
        end
        if edgeDetFilterTypeOpt == 'A'
            edgeDetFilterTypeName = 'average';
            edgeDetFilterSize = input('Size of the filter [3]: ');
            if isempty(edgeDetFilterSize)
                edgeDetFilterSize = 3;
            end
        elseif edgeDetFilterTypeOpt == 'D'
            edgeDetFilterTypeName = 'disk';
            edgeDetFilterSize = input('Radius of a disk-shaped filter [5]: ');
            if isempty(edgeDetFilterSize)
                edgeDetFilterSize = 5;
            end
        elseif edgeDetFilterTypeOpt == 'G'
            edgeDetFilterTypeName = 'gaussian';
            edgeDetFilterSize = input('Size of the filter [5]: ');
            if isempty(edgeDetFilterSize)
                edgeDetFilterSize = 5;
            end
            edgeDetFilterSigma = input('Standard deviation [0.5]: ');
            if isempty(edgeDetFilterSigma)
                edgeDetFilterSigma = 0.5;
            end
        elseif edgeDetFilterTypeOpt == 'L'
            edgeDetFilterTypeName = 'laplacian';
            edgeDetFilterShape = input('Shape of the Laplacian (min. = 0/ max. = 1) [0.2]: ');
            if isempty(edgeDetFilterShape) || (edgeDetFilterShape < 0 && edgeDetFilterShape > 1)
                edgeDetFilterShape = 0.2;
            end
        end
    end
elseif EdgeDetMethodOpt == 7 %Partial area effect method
    fprintf('Partial area effect edge detection method (subpixel detection technique) selected. \n');
    fprintf('Partial area effect edge detection method settings... \n');
    edgeDetDifThreshold = input('Threshold value (min. difference of intensity at both sides of a pixel to be considered an edge) [5]: ');
    if isempty(edgeDetDifThreshold)
        edgeDetDifThreshold = 5;
    end
    fprintf('Order of the edge to find: \n');
    fprintf('1. First order edges (straight lines) \n');
    fprintf('2. Second order edges \n');
    edgeDetOrder = input('Select an order [2]: ');
    if isempty(edgeDetOrder)
        edgeDetOrder = 2;
    end
    fprintf('Smoothing iterations needed to find final edges: \n');
    fprintf('0. Oriented to noise free images \n');
    fprintf('1. Oriented to low-noise images \n');
    fprintf('>1. Oriented to high-noise images \n');
    edgeDetSmoothingIter = input('Enter the number of smoothing iterations [1]: ');
    if isempty(edgeDetSmoothingIter)
        edgeDetSmoothingIter = 1;
    end
end
end

