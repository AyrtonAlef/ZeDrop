function [EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter,...
EdgeDetMethodOpt,edgeDetMetName,edgeDetDirection,edgeDetSigma,edgeDetThreshold,edgeDetFilter,edgeDetFilterTypeOpt,...
edgeDetFilterTypeName,edgeDetFilterSize,edgeDetFilterSigma,edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter,...
MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode,...
BaselineDetOpt,yDistWeight,marginBaselineDet,yAverageLR,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize] = SettingsGeneralEdge(EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter,...
EdgeDetMethodOpt,edgeDetMetName,edgeDetDirection,edgeDetSigma,edgeDetThreshold,edgeDetFilter,edgeDetFilterTypeOpt,...
edgeDetFilterTypeName,edgeDetFilterSize,edgeDetFilterSigma,edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter,...
MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode,...
BaselineDetOpt,yDistWeight,marginBaselineDet,yAverageLR,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize)
%SETTINGSGENERALEDGE Enable the user to edit edge detection and baseline
%determination settings
%   INPUT/OUTPUT:
%1) Pre-processing
% EnhanceContrast - Option to enhance contrast
% low_in - Low contrast limit for input image
% high_in - High contrast limit for input image
% low_out - Low contrast limit for output image
% high_out - High contrast limit for output image
% SmoothingFilter - Option for applying Gaussian smoothing filter
% sigmaSmoothFilter - Standard deviation of the 2-D Gaussian smoothing kernel
% sizeSmoothFilter - Filter size
%2) EdgeDetection
% EdgeDetMethodOpt - Edge detection method option (1. Sobel/ 2. Prewitt/ 3. Roberts/ 4. LoG/ 5. Canny/ 6. Zerocross/ 7. Partial area effect)
% edgeDetMetName- Name of the edge detection method
% edgeDetDirection - Direction of edges to detect (both/ horizontal/ vertical). Valid for Sobel, Prewitt and Roberts methods
% edgeDetSigma - Standard deviation of filter. Valid for LoG and Canny gradient methods
% edgeDetThreshold - Sensitivity threshold (min. = 0/ max. = 1). Valid for Canny gradient methods
% edgeDetFilter - Option to enable the definition of filter to zerocross method
% edgeDetFilterTypeOpt - Edge detection filter type option. Option for zerocross function
% edgeDetFilterType - Type of filter. Option for zerocross function
% edgeDetFilterSize - Filter size. Option for zerocross function and filter types: average, disk, gaussian and laplacian
% edgeDetFilterSigma - Standard deviation. Option for zerocross function and filter type: gaussian
% edgeDetFilterShape - Shape of the Laplacian (min. = 0/ max. = 1). Option for zerocross function and filter type: laplacian
% edgeDetDifThreshold - Minimum difference of intensity at both sides of a pixel to be considered an edge
% edgeDetOrder - Order of the edge to find(1. First order edges (straight lines)/ 2. Second order edges)
% edgeDetSmoothingIter - Smoothing iterations needed to find final edges (0. Oriented to noise free images/ 1. Oriented to low-noise images/ >1. Oriented to high-noise images) 
%3) Post-processing
% MorphOpen - Option to enable open morphological operation (remove small objects)
% openP - Maximum number of pixels in objects
% openConn - Pixel connectivity
% MorphBridge - Option to enable bridge morphological operation (bridge unconnected pixels)
% bridgeN - Number of times to perform bridge operation
% MorphClose - Option to enable close morphological operation (dilatation followed by erosion)
% closeN - Number of times to perform close operation [Inf]
% MorphThin - Option to enable thin morphological operation (thin objects to line by removing pixels from the boundary of objects
% thinN - Number of times to perfrom thin operation
% margin - Pixel distance between edges that should be considered part of same contour
% jump_dist - Distance value in pixels that will stop the search for new boundaries
% filterMode - Filter mode of drop edge profile (H - horizontal and V - Vertical)
%4) Baseline Determination
% BaselineDetOpt - Baselie determination option
% (0 - Manual/
% 1 - Automatic for reflective drop (Andersen and Taboryski 2017)/ 
% 2 - Automatic for non-reflective drop (Liu et al 2019)
% 3- Automatic for reflective and non-reflective drop (Akbari and Antonini 2021))
% yDistWeigth = 40; %Weigth given to the vertical distance between edge points and drop baseline center (BaselineDetOpt == 0)
% marginBaselineDet = 10; %Margin from baseline selected for the search of contact points (BaselineDetOpt == 0)
% yAverageLR = 1; %Average y coordinate between left and right baseline selections (BaselineDetOpt == 0)
% nBaselineDet - Integer number of points that are fitted. Option for Andersen and Taboryski 2017 code (BaselineDetOpt == 1)
% cutBaselineDet - Number of pixels from the top where the script should not look for reflections. Option for Andersen and Taboryski 2017 code (BaselineDetOpt == 1)
% percentBaselineDet - Percent used to calculate cutoff value. Option for Liu et al 2019 code (BaselineDetOpt == 2)
% fillBlankHorLines = 1; %Option for filling blank horizontal lines (1 = enable and 0 = disable). Option for Akbari and Antonini 2021 code (BaselineDetOpt == 3).
% marginMaskBaselineDet = 10; %Minimum pixel distance (in y coordinate) between points along the detected edge to be considered part of the drop profile (Suggestion = jum_dist). Option for Akbari and Antonini 2021 code (BaselineDetOpt == 3). 
% maskSize = 20; %Mask size (neighbour index)

while 1
    fprintf('- General settings of edge detection and baseline identification: \n');
    fprintf('1. Pre-processing (contrast and smoothing) \n');
    fprintf('2. Edge detection \n');
    fprintf('3. Post-processing (morphological operations, longest edge, sort left and right edges) \n');
    fprintf('4. Baseline determination \n');
    optEdgeConfig = input ('Enter an option [0 to exit]: ');
    if isempty(optEdgeConfig)
        optEdgeConfig = 0;
    end
    if optEdgeConfig == 0
        break;
    elseif optEdgeConfig == 1 %Pre-processing (contrast and smoothing)
        [EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter] = SettingsPreprocessing(EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter);
    elseif optEdgeConfig == 2 %Edge detection
        [EdgeDetMethodOpt,edgeDetMetName,edgeDetDirection,edgeDetSigma,edgeDetThreshold,edgeDetFilter,edgeDetFilterTypeOpt,edgeDetFilterTypeName,edgeDetFilterSize,edgeDetFilterSigma,edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter] = SettingsEdgeDetection(EdgeDetMethodOpt,edgeDetMetName,edgeDetDirection,edgeDetSigma,edgeDetThreshold,edgeDetFilter,edgeDetFilterTypeOpt,edgeDetFilterTypeName,edgeDetFilterSize,edgeDetFilterSigma,edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter);
    elseif optEdgeConfig == 3 %Post-processing (morphological operations)
        [MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode] = SettingsPostprocessing(MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode);
    elseif optEdgeConfig == 4 %Baseline determination
        [BaselineDetOpt,yDistWeight,marginBaselineDet,yAverageLR,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize] = SettingsBaselineDet(BaselineDetOpt,yDistWeight,marginBaselineDet,yAverageLR,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize);
    end
end
end

