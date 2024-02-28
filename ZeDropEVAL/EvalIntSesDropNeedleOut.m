function EvalIntSesDropNeedleOut()
%EVALINTSESDROPNEEDLEOUT Evaluation of intermittent (needle out) sessile 
% drop quasi-static experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used for the analysis of intermittent sessile drop (needle out) 
% quasi-static experiment. Based on Example1.m and Example2.m codes from 
% Andersen and Taboryski 2017.
% Created by Ayrton Pereira at 02/18/23.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%
%%{
%---------------------------------------------------------------------------
%                         0. Reset settings
%---------------------------------------------------------------------------
if ~isfile('VarIntSessDrop.mat') %Initialize default settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %4) Edge detection and baseline determination
    %- Pre-processing settings
    EnhanceContrast = 1; %Option to enhance contrast
    low_in = 0; %Low contrast limit for input image
    high_in = 1; %High contrast limit for input image
    low_out = 0; %Low contrast limit for output image
    high_out = 1; %High contrast limit for output image
    SmoothingFilter = 1; %Option for applying Gaussian smoothing filter
    %sigma = 2; %Standard deviation of the 2-D Gaussian smoothing kernel
    sigmaSmoothFilter = 2; %Standard deviation of the 2-D Gaussian smoothing kernel
    sizeSmoothFilter = 2*ceil(2*sigmaSmoothFilter)+1; %Filter size

    %- Edge detection settings
    %EdgeDetMethodOpt = 2; %Edge detection method option (1. Gradient Methods/ 2. Non-Gradient Methods)
    EdgeDetMethodOpt = 7; %Edge detection method option (1. Sobel/ 2. Prewitt/ 3. Roberts/ 4. LoG/ 5. Canny/ 6. Zerocross/ 7. Partial area effect)
    edgeDetMetName = 'zerocross'; %Name of the edge detection method
    edgeDetDirection = 'both'; %Direction of edges to detect (both/ horizontal/ vertical). Valid for Sober, Prewitt and Roberts methods
    edgeDetSigma = 2; %Standard deviation of filter. Valid for LoG and Canny methods
    edgeDetThreshold = []; %Sensitivity threshold (min. = 0/ max. = 1). Valid for Canny methods
    edgeDetFilter = 0; %Option to enable the definition of filter to zerocross method
    edgeDetFilterTypeOpt = 'G'; %Edge detection filter type option. Option for zerocross function
    edgeDetFilterTypeName = 'gaussian'; %Type of filter. Option for zerocross function
    edgeDetFilterSize = 5; %Filter size. Option for zerocross function and filter types: average, disk, gaussian and laplacian
    edgeDetFilterSigma = 0.5; %Standard deviation. Option for zerocross function and filter type: gaussian
    edgeDetFilterShape = 0.2; %Shape of the Laplacian (min. = 0/ max. = 1). Option for zerocross function and filter type: laplacian
    h = fspecial(edgeDetFilterTypeName,edgeDetFilterSize); %Filter. Option for zerocross function
    %edgeThreshold = 5; %Threshold value (for subpixelEdge function)
    edgeDetDifThreshold = 5; %Minimum difference of intensity at both sides of a pixel to be considered an edge
    %order = 2; %Specifies the order of the edge to find
    edgeDetOrder = 2; %Order of the edge to find(1. First order edges (straight lines)/ 2. Second order edges)
    %smoothingIter = 1; %Specifies how many smoothing iterations are needed to find final edges
    edgeDetSmoothingIter = 1; %Smoothing iterations needed to find final edges (0. Oriented to noise free images/ 1. Oriented to low-noise images/ >1. Oriented to high-noise images)

    %- Pos-processing settings
    %-- Identification of the longest edge
    margin = 10; % pixel distance between edges that should be considered part of same contour (for findlongestedge function)
    %-- Sort left and right edges
    jump_dist = margin; %distance value in pixels that will stop the search for new branches during sorting left and right edges
    filterMode = 'V'; %Filter mode of drop edge profile (H - horizontal and V - Vertical)
    %-- Morphological operations
    MorphOpen = 0; %Option to enable open morphological operation (remove small objects)
    openP = 100; %Maximum number of pixels in objects
    defaultopenP = 1; %Change if the openP value is modified
    openConn = 8; %Pixel connectivity
    MorphBridge = 0; %Option to enable bridge morphological operation (bridge unconnected pixels)
    bridgeN = Inf; %Number of times to perform bridge operation
    MorphClose = 0; %Option to enable close morphological operation (dilatation followed by erosion)
    closeN = Inf; %Number of times to perform close operation [Inf]
    MorphThin = 0; %Option to enable thin morphological operation (thin objects to line by removing pixels from the boundary of objects
    thinN = Inf; %Number of times to perfrom thin operation

    %- Baseline determination settings
    %BaselineOption = 1; %Option for the method usaed to determine baseline
    BaselineDetOpt = 3; %Baselie determination option
    % (0 - Manual determination/
    % 1 - Automatic for reflective drop (Andersen and Taboryski 2017)/
    % 2 - Automatic for non-reflective drop (Liu et al 2019)
    % 3- Automatic for reflective and non-reflective drop (Akbari and Antonini 2021))
    yDistWeight = 40; %Weigth given to the vertical distance between edge points and drop baseline center
    marginBaselineDet = 10; %Margin from baseline selected for the search of contact points (BaselineDetOpt == 0)
    yAverageLR = 1; %Average y coordinate between left and right baseline selections (BaselineDetOpt == 0)
    nBaselineDet = 10; %Integer number of points that are fitted. Option for Andersen and Taboryski 2017 code (BaselineDetOpt == 1)
    cutBaselineDet = 260; %Number of pixels from the top where the script should not look for reflections. Option for Andersen and Taboryski 2017 code (BaselineDetOpt == 1)
    percentBaselineDet = 0.98; %Percent used to calculate cutoff value. Option for Liu et al 2019 code (BaselineDetOpt == 2)
    fillBlankHorLines = 1; %Option for filling blank horizontal lines (1 = enable and 0 = disable). Option for Akbari and Antonini 2021 code (BaselineDetOpt == 3).
    marginMaskBaselineDet = 10; %Minimum pixel distance (in y coordinate) between points along the detected edge to be considered part of the drop profile (Suggestion = jum_dist). Option for Akbari and Antonini 2021 code (BaselineDetOpt == 3).
    maskSize = 20; %Mask size (neighbour index)

    %- Multiple image analysis settings
    limBaselineTilt = 1; %Baseline tilt limit in degrees allowed to consider an image suitable for subsequent fitting and determination of
    plotIntEdgeDetImage = 1; %Plot interval of resulted edge detection and baseline determination images

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %5) Average baselines found in all frames
    UseAverageBaseVec = 0; %Option to use average baseline using all baseline coordinates

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %6) Fitting data and calculating CA
    UseCircFit = 1; %Option for using circle fitting
    n_circfit = 100; %number of fitting points for circle fitting (Suggested: 100)
    PlotCircleImages = 1; %Option to plot circle fitting (0 = false, 1 = true)
    UsePolyFit = 1; %Option for using polynomial fitting
    n_polyfit = 10; %number of fitting points for polynomial fitting (Suggested: 10)
    poly_ord = 4; %polynomial order for polynomial fitting (Suggested: 4)
    PlotPolynomialImages = 1; %Option to plot polynomial fitting (0 = false, 1 = true)
    UsePolyFitDropen = 1; %Option for using polynomial fitting from Dropen
    n_polyfitDropen = 100; %Number of fitting points for polynomial fitting from Dropen (Suggested: 100)
    poly_ordDropen = 3; %polynomial order for polynomial fitting from Dropen (Suggested: (i) for CA < 60° -> poly_ord = 2 and (ii) for CA > 60° -> poly_ord = 3)
    PlotPolynomialDropenImages = 1; %Option to plot polynomial fitting (Dropen) (0 = false, 1 = true)
    UseEllipseFit = 1; %Option for using double sided elliptical fitting
    PlotEllipticImages = 1; %Option to plot ellipse fitting (0 = false, 1 = true)
    UseMaskMethod = 1; %Option for using mask method
    PlotMaskImages = 1; %Option to plot mask method (0 = false, 1 = true)

    plotIntFittingImages = 5; %Plot interval of fitting images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           SAVE VARIABLES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    routVarNames = {'EnhanceContrast','low_in','high_in','low_out','high_out','SmoothingFilter','sigmaSmoothFilter','sizeSmoothFilter','EdgeDetMethodOpt','edgeDetMetName','edgeDetDirection','edgeDetSigma','edgeDetThreshold','edgeDetFilter','edgeDetFilterTypeOpt','edgeDetFilterTypeName','edgeDetFilterSize','edgeDetFilterSigma','edgeDetFilterShape','edgeDetDifThreshold','edgeDetOrder','edgeDetSmoothingIter','margin','jump_dist','filterMode','MorphOpen','openP','defaultopenP','openConn','MorphBridge','bridgeN','MorphClose','closeN','MorphThin','thinN','BaselineDetOpt','yDistWeight','marginBaselineDet','yAverageLR','nBaselineDet','cutBaselineDet','percentBaselineDet','fillBlankHorLines','marginMaskBaselineDet','maskSize','limBaselineTilt','plotIntEdgeDetImage','UseAverageBaseVec','UseCircFit','n_circfit','PlotCircleImages','UsePolyFit','n_polyfit','poly_ord','PlotPolynomialImages','UsePolyFitDropen','n_polyfitDropen','poly_ordDropen','PlotPolynomialDropenImages','UseEllipseFit','PlotEllipticImages','UseMaskMethod','PlotMaskImages','plotIntFittingImages'};
    save('VarIntSessDrop.mat',routVarNames{:});
end
%}
%%
routVarNames = {'EnhanceContrast','low_in','high_in','low_out','high_out','SmoothingFilter','sigmaSmoothFilter','sizeSmoothFilter','EdgeDetMethodOpt','edgeDetMetName','edgeDetDirection','edgeDetSigma','edgeDetThreshold','edgeDetFilter','edgeDetFilterTypeOpt','edgeDetFilterTypeName','edgeDetFilterSize','edgeDetFilterSigma','edgeDetFilterShape','edgeDetDifThreshold','edgeDetOrder','edgeDetSmoothingIter','margin','jump_dist','filterMode','MorphOpen','openP','defaultopenP','openConn','MorphBridge','bridgeN','MorphClose','closeN','MorphThin','thinN','BaselineDetOpt','yDistWeight','marginBaselineDet','yAverageLR','nBaselineDet','cutBaselineDet','percentBaselineDet','fillBlankHorLines','marginMaskBaselineDet','maskSize','limBaselineTilt','plotIntEdgeDetImage','UseAverageBaseVec','UseCircFit','n_circfit','PlotCircleImages','UsePolyFit','n_polyfit','poly_ord','PlotPolynomialImages','UsePolyFitDropen','n_polyfitDropen','poly_ordDropen','PlotPolynomialDropenImages','UseEllipseFit','PlotEllipticImages','UseMaskMethod','PlotMaskImages','plotIntFittingImages'};
fprintf('---- EVALUATION OF INTERMITTENT (NEEDLE OUT) SESSILE DROP (ADVANCING AND RECEDING) EXPERIMENT ----\n');
%---------------------------------------------------------------------------
%                1.       Loading of data
%---------------------------------------------------------------------------
%- Load variables
load VarIntSessDrop.mat
%- Read files
fprintf('--------------------------------- 1. DATA LOADING ------------------------------------ \n');
fprintf('Select images... \n');
[fileNames, inputPath] = uigetfile({'C:\*.*'},'MultiSelect','on'); %Selecting images
fprintf('Select a directory for exporting results ... \n');
outputPath = uigetdir(inputPath);
outputName = input('Enter a prefix name for the exported data files: ','s');
outputFmt = '.txt'; %Data output format
nStops = input('Enter the no. of stops performed during the test (advancing and receding) [20]: ');
if isempty(nStops)
    nStops = 20;
end
tCapt = input('Enter the capture time at each stop in s [10]: '); %Considering a capture of 1 image per second
if isempty(tCapt)
    tCapt = 10;
end
%- Creation of variables
if ~iscell(fileNames) %Selection of a single image
    a = fileNames;
    fileNames = cell(1);
    fileNames{1} = a;
end
nImages = size(fileNames,2); %Number of images
nImagesStops = nImages/nStops; %Number of images captured per stop
im = cell(1,length(fileNames));
volExp = zeros(nImages,1); %Drop volume in uL

%- Loading images and drop volume
ansVol = input('Load volume information from the image name (I) or from a file (F) [I]? ','s');
if isempty(ansVol) || (ansVol ~= 'I' && ansVol ~= 'F')
    ansVol = 'I';
end
for i = 1:nImages
    im{i} = imread([inputPath,char(fileNames(i))]); %Reading images
    if ansVol == 'I' %Extract volume information from imagename
        %extractVol = extractBetween(fileNames{i},"V(uL)_",".tif"); %Extract volume from file name
        extractVol = extractBetween(fileNames{i},"V(uL)_","_"); %Extract volume from file name
        volExp(i) = str2num(extractVol{1}); %Convert extracted string to volume value in uL
    end
end
if ansVol == 'F'
    fprintf('Select volume information file... \n')
    [fileNameVol, inputPathVol] = uigetfile({'*.dat;*.csv;*.txt;*.LOG','Data File(*.dat,*.csv,*.txt,*.LOG)';'*.*','All Files (*.*)'},'Select a file',inputPath);
    fullfileNameVol = fullfile(inputPathVol,fileNameVol);
    fprintf('Reading volume information file... \n')
    DataVol = readtable(fullfileNameVol,'VariableNamingRule','preserve'); %Import data from file
    col = input('Enter the column which contain the volume data [3]: ');
    if isempty(col)
        col = 3;
    end
    volExp = DataVol{:,col}; %Extract volume information
end
fprintf('Data loading completed! \n')
%%
%---------------------------------------------------------------------------
%                2.       Image scale
%---------------------------------------------------------------------------
%scale = 137.2; %Scale in pixels/mm
% Determination of image scale
fprintf('------------------------- 2. SCALE FACTOR DETERMINATION ------------------------------- \n');
while 1
    %- Check the existence of the scale factor file (ScaleFactor.mat)
    if ~isfile('ScaleFactor.mat') %First time running the code
        fprintf('An scale factor file was not found in the current folder. \n');
        fprintf('To continue it is required to enter a scale factor... \n');
    else
        load ScaleFactor.mat scale  %Load scale factor from root folder
        fprintf('An scale factor file was found in the root folder (scale = %.2f pxs/mm). \n',scale);
        ansScale = input('Do you like to continue with it (Y/N) [Y]? ', 's'); 
        if isempty(ansScale)
            ansScale = 'Y';
        end
        if ansScale == 'Y'
            break;
        end
    end
    ansScalMode= input('Do you like to type the scale factor manually (M) or define it from an image (I) [I]? ','s');
    if isempty(ansScalMode)
        ansScalMode = 'I';
    end
    if ansScalMode == 'M'
        scale = input('Enter the scale factor [pxs/mm]: ');
        
    else
        fprintf('Select an image for scale factor determination [pxs/mm]... \n');
        [fileScaleName, scalePath] = uigetfile({'*.*'},'Select a file',inputPath); %Select image for determining scale factor
        imScale = imread([scalePath,fileScaleName]);
        scale = setScaleFactor(imScale);
    end
    save ScaleFactor.mat scale %Saving scale factor
    fprintf('Scale factor saved. \n');
    break;
end
fprintf('Scale factor determination completed! \n')
%%
%---------------------------------------------------------------------------
%                3.       Crop images
%---------------------------------------------------------------------------
% Crop images
fprintf('---------------------------------- 3. IMAGE CROP ------------------------------------ \n');
ansCrop = input('Do you like to crop images (Y/N) [N]? ','s');
if isempty(ansCrop)
    ansCrop = 'N';
end
if ansCrop == 'Y'
    % Selecting reference image for determining ROI
    fprintf('Select an image for the determination of the region of interest... \n');
    [fileCropName, cropPath] = uigetfile({'*.*'},'Select a file',inputPath); %Select image for determining crop region
    imCrop = imread([cropPath,fileCropName]);
    [imCropRef, rect] = imcrop(imCrop);
    close(gcf)

    %Cropping images
    fprintf('Please wait, cropping images... \n')
    for i = 1:nImages
        im{i} = imcrop(im{i},rect); %crop image
    end
    fprintf('Image crop completed! \n')
end
%%
%---------------------------------------------------------------------------
%                4.       Edge detection
%---------------------------------------------------------------------------
fprintf('----------------- 4. EDGE DETECTION AND BASELINE DETERMINATION -------------------- \n');
while 1
    % Edge detection and baseline determination settings
    ansEdgeConfig = input('Do you like to change edge detection and baseline determination settings (Y/N) [N]? ','s');
    if isempty(ansEdgeConfig) || (ansEdgeConfig ~= 'Y' && ansEdgeConfig ~= 'N')
        ansEdgeConfig = 'N';
    end
    if ansEdgeConfig == 'Y'
        [EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter,...
            EdgeDetMethodOpt,edgeDetMetName,edgeDetDirection,edgeDetSigma,edgeDetThreshold,edgeDetFilter,edgeDetFilterTypeOpt,...
            edgeDetFilterTypeName,edgeDetFilterSize,edgeDetFilterSigma,edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter,...
            MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode,...
            BaselineDetOpt,yDistWeight,marginBaselineDet,yAverageLR,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize] = SettingsGeneralEdge(EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter,...
            EdgeDetMethodOpt,edgeDetMetName,edgeDetDirection,edgeDetSigma,edgeDetThreshold,edgeDetFilter,edgeDetFilterTypeOpt,...
            edgeDetFilterTypeName,edgeDetFilterSize,edgeDetFilterSigma,edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter,...
            MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode,...
            BaselineDetOpt,yDistWeight,marginBaselineDet,yAverageLR,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize);
    end
    %----------------------------------------------------------------------
    %- Check edge detection and baseline determination routine
    %----------------------------------------------------------------------
    while 1
        close all
        ansTestEdgeDetection = input('Do you like to check edge detection and baseline determination routine (Y/N) [N]? ','s');
        if isempty(ansTestEdgeDetection) || (ansTestEdgeDetection ~= 'Y' && ansTestEdgeDetection ~= 'N')
            ansTestEdgeDetection = 'N';
        end
        if ansTestEdgeDetection == 'Y'
            RefImageEdgeDetection = input('Enter an image number [1]: ');
            if isempty(RefImageEdgeDetection)
                RefImageEdgeDetection = 1;
            end
            %1) Load reference image
            fprintf('Loading reference image... \n');
            imRef = im{RefImageEdgeDetection};
            %- Show reference image
            refImage = figure;
            imshow(imRef)
            set(refImage,'Name','Reference Image','Units','normalized','Position',[0,0.6,0.3,0.3]);

            %2) Pre-processing of image
            while 1
                fprintf('Pre-processing of image... \n');
                imPreProc = imRef;
                if EnhanceContrast
                    if low_in == 0 && high_in == 1 && low_out == 0 && high_out == 1
                        imPreProc = imadjust(imPreProc); %Adjust image intensity values or colormap. Saturates the bottom 1% and the top 1% of all pixel values. Increases the contrast of the output image
                    else
                        imPreProc = imadjust(imPreProc,[low_in high_in],[low_out high_out]);
                    end
                end
                if SmoothingFilter
                    imPreProc = imgaussfilt(imPreProc,sigmaSmoothFilter,"FilterSize",sizeSmoothFilter); %Apply Gaussian smoothing filter (Gaussian smoothing kernel). Sigma is the standard deviation of the Gaussian distribution (increasing sigma, increases the smoothing/blur effect).
                end
                %- Show pre-processed image
                preProcImage = figure;
                imshow(imPreProc)
                set(preProcImage,'Name','Pre-processed Image','Units','normalized','Position',[0.31,0.6,0.3,0.3]);
                %- Ask for changes of pre-processing settings
                ansPreProcessing = input('Do you like to change pre-processing settings and redo pre-processing step (Y/N) [N]? ','s');
                if isempty(ansPreProcessing)
                    ansPreProcessing = 'N';
                end
                if ansPreProcessing == 'Y' %Open pre-processing image settings
                    [EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter] = SettingsPreprocessing(EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter);
                    close(preProcImage)
                else
                    break;
                end
            end

            %3) Edge Detection
            fprintf('Edge detection... \n');
            while 1
                checkEdgeDetImage = [];
                ansCheckEdgeDetection = input('Do you like to run all edge detection methods available (Y/N) [N]? ','s');
                if isempty(ansCheckEdgeDetection)
                    ansCheckEdgeDetection = 'N';
                end
                if ansCheckEdgeDetection == 'Y' %Run all edge detection methods available
                    checkEdgeDetImage = figure;
                    set(checkEdgeDetImage,'Name','Edge detection methods','Units','normalized','Position',[0.1,0.1,0.8,0.8]);
                    titleColor = 'black';
                    %Sobel
                    fprintf('Sobel edge detection method... \n');
                    edgeBinImageSobel = edge(imPreProc,'Sobel',edgeDetThreshold,edgeDetDirection); %Edge binary image
                    subplot(3,3,1)
                    imshow(edgeBinImageSobel);
                    if EdgeDetMethodOpt == 1 %Change title color if this method was selected
                        titleColor = 'red';
                        edgeBinImage = edgeBinImageSobel;
                    end
                    title(["(1) Sobel",strcat("Threshold: ",num2str(edgeDetThreshold),"/ Direction: ",edgeDetDirection)],'Color',titleColor)
                    titleColor = 'black';
                    %Prewitt
                    fprintf('Prewitt edge detection method... \n');
                    edgeBinImagePrewitt = edge(imPreProc,'Prewitt',edgeDetThreshold,edgeDetDirection); %Edge binary image
                    subplot(3,3,2)
                    if EdgeDetMethodOpt == 2 %Change title color if this method was selected
                        titleColor = 'red';
                        edgeBinImage = edgeBinImagePrewitt;
                    end
                    imshow(edgeBinImagePrewitt);
                    title(["(2) Prewitt",strcat("Threshold: ",num2str(edgeDetThreshold),"/ Direction: ",edgeDetDirection)],'Color',titleColor)
                    titleColor = 'black';
                    %Roberts
                    fprintf('Roberts edge detection method... \n');
                    edgeBinImageRoberts = edge(imPreProc,'Roberts',edgeDetThreshold,edgeDetDirection); %Edge binary image
                    subplot(3,3,3)
                    imshow(edgeBinImageRoberts);
                    if EdgeDetMethodOpt == 3 %Change title color if this method was selected
                        titleColor = 'red';
                        edgeBinImage = edgeBinImageRoberts;
                    end
                    title(["(3) Roberts",strcat("Threshold: ",num2str(edgeDetThreshold),"/ Direction: ",edgeDetDirection)],'Color',titleColor)
                    titleColor = 'black';
                    %Laplace of Gaussian (LoG)
                    fprintf('Laplacian of Gaussian edge detection method... \n');
                    edgeBinImageLoG = edge(imPreProc,'log',edgeDetThreshold,edgeDetSigma); %Edge binary image
                    subplot(3,3,4)
                    imshow(edgeBinImageLoG);
                    if EdgeDetMethodOpt == 4 %Change title color if this method was selected
                        titleColor = 'red';
                        edgeBinImage = edgeBinImageLoG;
                    end
                    title(["(4) LoG",strcat("Threshold: ",num2str(edgeDetThreshold),"/ Sigma: ",num2str(edgeDetSigma))],'Color',titleColor)
                    titleColor = 'black';
                    %Canny
                    fprintf('Canny edge detection method... \n');
                    edgeBinImageCanny = edge(imPreProc,'Canny',edgeDetThreshold,edgeDetSigma); %Edge binary image
                    subplot(3,3,5)
                    imshow(edgeBinImageCanny);
                    if EdgeDetMethodOpt == 5 %Change title color if this method was selected
                        titleColor = 'red';
                        edgeBinImage = edgeBinImageCanny;
                    end
                    title(["(5) Canny",strcat("Threshold: ",num2str(edgeDetThreshold),"/ Sigma: ",num2str(edgeDetSigma))],'Color',titleColor)
                    titleColor = 'black';
                    %Zero-crossing
                    fprintf('Zero-crossing edge detection method... \n');
                    if edgeDetFilter %FIlter modification and creation
                        edgeBinImageZeroCross = edge(imPreProc,'zerocross',edgeDetThreshold,h); %Edge binary image
                        subplot(3,3,6)
                        imshow(edgeBinImageZeroCross);
                        if EdgeDetMethodOpt == 6 %Change title color if this method was selected
                            titleColor = 'red';
                            edgeBinImage = edgeBinImageZeroCross;
                        end
                        title(['(6) Zero-Crossing',strcat("Threshold: ",num2str(edgeDetThreshold),"/ FilterName: ",edgeDetFilterTypeName)],'Color',titleColor);
                    else
                        edgeBinImageZeroCross = edge(imPreProc,'zerocross',edgeDetThreshold);
                        subplot(3,3,6)
                        imshow(edgeBinImageZeroCross);
                        if EdgeDetMethodOpt == 6 %Change title color if this method was selected
                            titleColor = 'red';
                            edgeBinImage = edgeBinImageZeroCross;
                        end
                        title(['(6) Zero-Crossing',strcat("Threshold: ",num2str(edgeDetThreshold))],'Color',titleColor);
                    end
                    titleColor = 'black';
                    %Partial area effect
                    fprintf('Partial area effect edge detection method... \n');
                    [edges, RI] = subpixelEdges(imPreProc, edgeDetDifThreshold,'Order',edgeDetOrder,'SmoothingIter',edgeDetSmoothingIter); %find edge
                    %- Compute the edge binary image
                    %-- Fast but less precision process
                    %edgeBinImagePAE = zeros(size(imPreProc));
                    %edgeBinImagePAE(edges.position) = 1;
                    %-- Computes a high resolution edge binary image
                    edgeBinImagePAE = subpixelImage(edges,size(imPreProc),1);
                    subplot(3,3,8)
                    imshow(edgeBinImagePAE);
                    if EdgeDetMethodOpt == 7 %Change title color if this method was selected
                        titleColor = 'red';
                        edgeBinImage = edgeBinImagePAE;
                    end
                    title(['(7) PAE',strcat("ThresholdDif: ",num2str(edgeDetDifThreshold)),strcat("Order: ",num2str(edgeDetOrder),"/ SmoothingIter: ",num2str(edgeDetSmoothingIter))],'Color',titleColor)
                    titleColor = 'black';
                else %Run only the selected edge detection method
                    if EdgeDetMethodOpt == 1 || EdgeDetMethodOpt == 2 || EdgeDetMethodOpt == 3 %Sobel, Prewitt or Roberts
                        if EdgeDetMethodOpt == 1
                            fprintf('Sobel edge detection method... \n');
                        elseif EdgeDetMethodOpt == 2
                            fprintf('Prewitt edge detection method... \n');
                        elseif EdgeDetMethodOpt == 3
                            fprintf('Roberts edge detection method... \n');
                        end
                        edgeBinImage = edge(imPreProc,edgeDetMetName,edgeDetThreshold,edgeDetDirection); %Edge binary image
                    elseif EdgeDetMethodOpt == 4 %LoG
                        fprintf('Laplacian of Gaussian edge detection method... \n');
                        edgeBinImage = edge(imPreProc,edgeDetMetName,edgeDetThreshold,edgeDetSigma); %Edge binary image
                    elseif EdgeDetMethodOpt == 5 %Canny
                        fprintf('Canny edge detection method... \n');
                        edgeBinImage = edge(imPreProc,edgeDetMetName,edgeDetThreshold,edgeDetSigma); %Edge binary image
                    elseif EdgeDetMethodOpt == 6 %Zerocross
                        fprintf('Zero-crossing edge detection method... \n');
                        if edgeDetFilter %FIlter modification and creation
                            edgeBinImage = edge(imPreProc,edgeDetMetName,edgeDetThreshold,h); %Edge binary image
                        else
                            edgeBinImage = edge(imPreProc,edgeDetMetName,edgeDetThreshold); %Edge binary image
                        end
                    elseif EdgeDetMethodOpt == 7 %Partial Area Effect
                        fprintf('Partial area effect edge detection method... \n');
                        [edges, RI] = subpixelEdges(imPreProc, edgeDetDifThreshold,'Order',edgeDetOrder,'SmoothingIter',edgeDetSmoothingIter); %find edge
                        %- Compute the edge binary image
                        %-- Fast but less precision process
                        %edgeBinImage = zeros(size(imPreProc));
                        %edgeBinImage(edges.position) = 1;
                        %-- Computes a high resolution edge binary image
                        edgeBinImage = subpixelImage(edges,size(imPreProc),1);
                    end
                end
                %- Show edge detection image
                edgeDetImage = figure;
                imshow(edgeBinImage)
                set(edgeDetImage,'Name',strcat('Edge detected binary image (EdgeDetMetOpt = ',num2str(EdgeDetMethodOpt),')'),'Units','normalized','Position',[0.62,0.6,0.3,0.3]);

                %- Ask for changes of edge detection settings
                ansEdgeDetection = input('Do you like to change edge detection settings and redo edge detection step (Y/N) [N]? ','s');
                if isempty(ansEdgeDetection)
                    ansEdgeDetection = 'N';
                end
                if ansEdgeDetection == 'Y' %Open edge detection settings
                    [EdgeDetMethodOpt,edgeDetMetName,edgeDetDirection,edgeDetSigma,edgeDetThreshold,edgeDetFilter,edgeDetFilterTypeOpt,edgeDetFilterTypeName,...
                        edgeDetFilterSize,edgeDetFilterSigma,edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter] = SettingsEdgeDetection(EdgeDetMethodOpt,edgeDetMetName,...
                        edgeDetDirection,edgeDetSigma,edgeDetThreshold,edgeDetFilter,edgeDetFilterTypeOpt,edgeDetFilterTypeName,edgeDetFilterSize,edgeDetFilterSigma,...
                        edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter);
                    if ~isempty(checkEdgeDetImage)
                        close(checkEdgeDetImage)
                        checkEdgeDetImage = [];
                    end
                    close(edgeDetImage)
                else
                    if ~isempty(checkEdgeDetImage)
                        close(checkEdgeDetImage)
                        checkEdgeDetImage = [];
                    end
                    break;
                end
            end

            %4) Post-processing
            fprintf('Post-processing of image... \n');
            %(i) Morphological operations
            while 1
                fprintf('Running morphological operations... \n');
                edgeBinPostProcImage = edgeBinImage;
                if MorphOpen
                    fprintf('Running open morphological operations... \n');
                    if defaultopenP
                        openP = round(size(edgeBinPostProcImage,1)/2);
                    end
                    edgeBinPostProcImage = bwareaopen(edgeBinPostProcImage,openP,openConn);
                elseif MorphBridge
                    fprintf('Running bridge morphological operations... \n');
                    edgeBinPostProcImage = bwmorph(edgeBinPostProcImage,'bridge',bridgeN);
                elseif MorphClose
                    fprintf('Running close morphological operations... \n');
                    edgeBinPostProcImage = bwmorph(edgeBinPostProcImage,'close',closeN);
                elseif MorphThin
                    fprintf('Running thin morphological operations... \n');
                    edgeBinPostProcImage = bwmorph(edgeBinPostProcImage,'thin',thinN);
                end
                %- Show edge detection image
                postProcMorphImage = figure;
                imshow(edgeBinPostProcImage)
                set(postProcMorphImage,'Name','Edge detected post-processed binary image','Units','normalized','Position',[0,0.3,0.3,0.3]);
                %- Ask for changes of morphological operations settings
                ansPostProcessing = input('Do you like to change morphological operation settings and redo morphological operation step (Y/N) [N]? ','s');
                if isempty(ansPostProcessing)
                    ansPostProcessing = 'N';
                end
                if ansPostProcessing == 'Y' %Open post-processing image settings
                    [MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode] = SettingsPostprocessing(MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode);
                    close(postProcMorphImage)
                else
                    break;
                end
            end

            %(ii) Identification of the longest edge
            while 1
                fprintf('Identification of the longest edge... \n');
                if EdgeDetMethodOpt ~= 7 || MorphOpen || MorphBridge || MorphClose || MorphThin
                    B = bwboundaries(edgeBinPostProcImage);
                    edgeCoords = B{1}; %Initialize edgeCoords array
                    for j = 2:length(B)
                        edgeCoords = cat(1,edgeCoords,B{j});
                    end
                    edges = EdgePixel;
                    edges.x = edgeCoords(:,2);
                    edges.y = edgeCoords(:,1);
                    edges.position = sub2ind(size(imPreProc),edges.y,edges.x); %Convert subscripts to linear indices
                    edges.curv = LineCurvature2D([edges.x,edges.y]);
                    N = LineNormals2D([edges.x,edges.y]);
                    edges.nx = N(:,1);
                    edges.ny = N(:,2);
                    edges.i0 = zeros(length(edges.position),1);
                    edges.i1 = zeros(length(edges.position),1);
                end
                longestedge = findlongestedge(edges,size(edgeBinPostProcImage),margin); % Selection of the longest edge
                %- Show longest edge image
                postProcLongEdgeImage = figure;
                imshow(imRef)
                hold on
                plot(longestedge.x,longestedge.y,'.')
                set(postProcLongEdgeImage,'Name','Identification of the longest edge','Units','normalized','Position',[0.31,0.3,0.3,0.3]);
                %- Ask for changes of identification of the longest edge settings
                ansPostProcessing = input('Do you like to change the identification of the longest edge settings and redo this step (Y/N) [N]? ','s');
                if isempty(ansPostProcessing)
                    ansPostProcessing = 'N';
                end
                if ansPostProcessing == 'Y' %Open post-processing image settings
                    [MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode] = SettingsPostprocessing(MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode);
                    close(postProcLongEdgeImage)
                else
                    break;
                end
            end

            %(iii) Identification and sort of left and right edges
            while 1
                fprintf('Identification and sort of left and right edge... \n');
                [edgeL,edgeR] = leftrightedges(longestedge,jump_dist); % split edge into left and right and sort it
                %- Show left and right edge image
                postProcImage = figure;
                imshow(imRef)
                hold on
                plot(edgeL.x,edgeL.y,'.') %Plot left edge
                plot(edgeR.x,edgeR.y,'.') %Plot right edge
                set(postProcImage,'Name','Identification and sort of left and right edges','Units','normalized','Position',[0.62,0.3,0.3,0.3]);
                %- Ask for changes of identification and sort of left and right edges settings
                ansPostProcessing = input('Do you like to change identification and sort of left and right edges settings and redo this step (Y/N) [N]? ','s');
                if isempty(ansPostProcessing)
                    ansPostProcessing = 'N';
                end
                if ansPostProcessing == 'Y' %Open post-processing image settings
                    [MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode] = SettingsPostprocessing(MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode);
                    close(postProcImage)
                else
                    break;
                end
            end

            %4) Baseline determination
            while 1
                fprintf('Baseline determination... \n');
                %- Baseline determination
                if BaselineDetOpt == 0  %- Manual baseline detection (from Akbari and Antonini 2021 code)
                    fprintf('Manual baseline determination... \n');
                    [baseL,baseR,dropCenter] = ManualSelection(imRef); %Determination of reference points
                    if yAverageLR %--- Use of average y left and right baseline (dropcenter(2))
                        ybaseL = dropCenter(2);
                        ybaseR = dropCenter(2);
                    else %--- Use of  y left and right baseline separately
                        ybaseL = baseL(2);
                        ybaseR = baseR(2);
                    end
                    %-- Edge detection margin for CP search
                    indexLMin = find(abs((edgeL.y-(ybaseL-(marginBaselineDet/2)))) == min(abs(edgeL.y-(ybaseL-(marginBaselineDet/2)))),1);
                    indexLMax = find(abs((edgeL.y-(ybaseL+(marginBaselineDet/2)))) == min(abs(edgeL.y-(ybaseL+(marginBaselineDet/2)))),1);
                    indexRMin = find(abs((edgeR.y-(ybaseR-(marginBaselineDet/2)))) == min(abs(edgeR.y-(ybaseR-(marginBaselineDet/2)))),1);
                    indexRMax = find(abs((edgeR.y-(ybaseR+(marginBaselineDet/2)))) == min(abs(edgeR.y-(ybaseR+(marginBaselineDet/2)))),1);
                    if isempty(indexLMin) || isempty(indexLMax) || isempty(indexRMin) || isempty(indexRMax)
                        fprintf('Baseline cannot be identified! \n');
                        fprintf('Try to modify yDistWeigth or marginBaselineDet. \n');
                        fprintf('Suggestion1: Check that the entire edge is being identified. Review the edge detection method used. \n')
                        %-- Definition of CPs
                        x0L = edgeL.x(1);
                        y0L = edgeL.y(1);
                        x0R = edgeR.x(end);
                        y0R = edgeR.y(end);
                    else
                        %-- CP search
                        distL = sqrt((dropCenter(1)-edgeL.x(indexLMin:indexLMax)).^2 + ((ybaseL-edgeL.y(indexLMin:indexLMax)).^2).*yDistWeight); %Distance between edge points and drop baseline center
                        indexL = find(distL == min(distL)) + indexLMin - 1;
                        distR = sqrt((dropCenter(1)-edgeR.x(indexRMin:indexRMax)).^2 + ((ybaseR-edgeR.y(indexRMin:indexRMax)).^2).*yDistWeight); %Distance between edge points and drop baseline center
                        indexR = find(distR == min(distR)) + indexRMin - 1;
                        %-- Definition of CPs
                        x0L = edgeL.x(indexL);
                        y0L = edgeL.y(indexL);
                        x0R = edgeR.x(indexR);
                        y0R = edgeR.y(indexR);
                    end
                    if (isempty(x0L) || isempty(y0L)) && (isempty(x0R) || isempty(y0R))
                        fprintf('Left and right baseline points cannot be identified! \n');
                        fprintf('Try to modify yDistWeigth or marginBaselineDet. \n');
                        fprintf('Suggestion1: Check that the entire edge is being identified. Review the edge detection method used. \n')
                        %-- Definition of CPs
                        x0L = edgeL.x(1);
                        y0L = edgeL.y(1);
                        x0R = edgeR.x(end);
                        y0R = edgeR.y(end);
                    elseif isempty(x0L) || isempty(y0L)
                        fprintf('Left baseline point cannot be identified! \n');
                        fprintf('Try to modify yDistWeigth or marginBaselineDet. \n');
                        fprintf('Suggestion1: Check that the entire edge is being identified. Review the edge detection method used. \n')
                        %-- Definition of CPs
                        x0L = edgeL.x(1);
                        y0L = edgeL.y(1);
                    elseif isempty(x0R) || isempty(y0R)
                        fprintf('Right baseline point cannot be identified! \n');
                        fprintf('Try to modify yDistWeigth or marginBaselineDet. \n');
                        fprintf('Suggestion1: Check that the entire edge is being identified. Review the edge detection method used. \n')
                        %-- Definition of CPs
                        x0R = edgeR.x(1);
                        y0R = edgeR.y(1);
                    end
                elseif BaselineDetOpt == 1 %- Automatic baseline detection for reflective drops (from Andersen et al 2017 code)
                    fprintf('Automatic baseline detection for reflective drops... \n');
                    [x0L,y0L,indexL] = findreflection([edgeL.x,edgeL.y],nBaselineDet,cutBaselineDet); %find reflection
                    [x0R,y0R,indexR] = findreflection([edgeR.x,edgeR.y],nBaselineDet,cutBaselineDet);
                elseif BaselineDetOpt == 2  %- Automatic baseline detection for non-reflective drops (find baseline(z-plane)) (from Liu et al 2019 code)
                    fprintf('Automatic baseline detection for non-reflective drops... \n');
                    [p,pos0L,pos0R,pos0L_cut,pos0R_cut,index_left,index_right,index_cut_left,index_cut_right] = plane_mod([edgeL.x,edgeL.y],[edgeR.x,edgeR.y],percentBaselineDet);
                    %- Use of percent
                    x0L = pos0L_cut(:,1);
                    y0L = pos0L_cut(:,2);
                    x0R = pos0R_cut(:,1);
                    y0R = pos0R_cut(:,2);
                elseif BaselineDetOpt == 3 %- Automatic baseline detection for non-reflective drops (find baseline(z-plane)) (from Akbari and Antonini 2021 code)(mask method)
                    fprintf('Automatic baseline detection for reflective and non-reflective... \n');
                    [I_New,AXxL,AXyL,AXxR,AXyR,Apex_x,Apex_y,Left_edgeTx,Left_edgeDx,Right_edgeTx,Right_edgeDx,Num_Down_L, Num_Up_L,Num_Down_R,...
                        Num_Up_R,End_x_L,End_x_R] = ImageAnalysis2(imPreProc,edgeL,edgeR,fillBlankHorLines,marginMaskBaselineDet);
                    %Automatic identification of contact points using mask method
                    [b,alphac,alphamL,alphamR,x_CPL,y_CPL,x_CPR,y_CPR,indexCPL,indexCPR,Contact_angle_left,Contact_angle_right,Contact_angle,TiltAngle_S] ...
                        = CircleMask2(imPreProc,I_New,maskSize,AXxL,AXyL,AXxR,AXyR,Apex_x,Apex_y,Left_edgeTx,Left_edgeDx,Right_edgeTx,Right_edgeDx,...
                        Num_Down_L, Num_Up_L,Num_Down_R,Num_Up_R,End_x_L,End_x_R);
                    if (isempty(x_CPL) || isempty(y_CPL)) && (isempty(x_CPR) || isempty(y_CPR))
                        fprintf('Left and right baseline points cannot be identified! \n');
                        fprintf('Suggestion1: Check that the entire edge is continuous. Review the edge detection method used. \n')
                        %-- Definition of CPs
                        x0L = edgeL.x(1);
                        y0L = edgeL.y(1);
                        x0R = edgeR.x(end);
                        y0R = edgeR.y(end);
                    elseif isempty(x_CPL) || isempty(y_CPL)
                        fprintf('Left baseline point cannot be identified! \n');
                        fprintf('Suggestion1: Check that the entire edge is continuous. Review the edge detection method used. \n')
                        %-- Definition of CPs
                        x0L = edgeL.x(1);
                        y0L = edgeL.y(1);
                    elseif isempty(x_CPR) || isempty(y_CPR)
                        fprintf('Right baseline point cannot be identified! \n');
                        fprintf('Suggestion1: Check that the entire edge is continuous. Review the edge detection method used. \n')
                        %-- Definition of CPs
                        x0R = edgeR.x(1);
                        y0R = edgeR.y(1);
                    else
                        y0L = x_CPL;
                        x0L = y_CPL;
                        y0R = x_CPR;
                        x0R = y_CPR;
                    end
                end
                %- Show baseline and highlighe contact points
                baselineImage = figure;
                imshow(imRef)
                hold on
                plot(edgeL.x,edgeL.y,'.') %Plot left edge
                plot(edgeR.x,edgeR.y,'.') %Plot right edge
                t = linspace(-3,3);
                plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r--','LineWidth',2) %Plot baseline
                plot(x0L,y0L,'yx','MarkerSize',10,'LineWidth',2) %Plot left contact point
                plot(x0R,y0R,'yx','MarkerSize',10,'LineWidth',2) %Plot right contact point
                set(baselineImage,'Name','Baseline level and identification of contact points','Units','normalized','Position',[0,0.05,0.3,0.3]);
                %- Ask for changes of baseline determination
                ansBaselineDet = input('Do you like to change baseline determination settings and redo this step (Y/N) [N]? ','s');
                if isempty(ansBaselineDet)
                    ansBaselineDet = 'N';
                end
                if ansBaselineDet == 'Y' %Open baseline determination settings
                    [BaselineDetOpt,yDistWeight,marginBaselineDet,yAverageLR,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize] = SettingsBaselineDet(BaselineDetOpt,yDistWeight,marginBaselineDet,yAverageLR,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize);
                    close(baselineImage)
                else
                    break;
                end
            end
            fprintf('Edge detection and baseline determination routine completed! \n')
            ansContEdgeDetection = input('Do you like to proceed with the current settings for the analysis of the rest of the images (Y/N) [Y]? ','s');
            if isempty(ansContEdgeDetection)
                ansContEdgeDetection = 'Y';
            end
            if ansContEdgeDetection == 'Y'
                break;
            end
        else
            break;
        end
    end
    close all
    %----------------------------------------------------------------------
    %- Edge detection and baseline determination of all images
    %----------------------------------------------------------------------
    fprintf('Analysis of all images... \n')
    %- Predefine matrices prior to loop
    %listExcTiltImages = "";
    BaselineTilt = zeros(nImages,1); % baseline tilt
    EdgeCell = cell(1,nImages); % cell structure to store detected edges
    BaseVec = zeros(nImages,4); % matrix to store baseline/contact point coordinates [x0L,y0L,x0R,y0R]

    %- Settings to export edge detected images
    dirEdgeDetection = 'EdgeDetection'; %Creation of EdgeDetection directory
    mkdir(outputPath,dirEdgeDetection);
    plotEdgeDetectionPath = fullfile(outputPath,dirEdgeDetection);
    %contExcTilt = 0;
    %nExcTiltImages = 0;
    limBaselineTilt = input('Max. limit of baseline tilt (°) [1]: ');
    if isempty(limBaselineTilt)
        limBaselineTilt = 1;
    end
    plotIntEdgeDetImages = input('Plot interval of edge detected images [5]: ');
    if isempty(plotIntEdgeDetImages)
        plotIntEdgeDetImages = 5;
    end  

    %- Analysis of all images
    imPreProc = cell(1,length(fileNames));
    for i = 1:nImages
        fprintf('Edge detection image %d of %d... \n',i,nImages);
        %imProc = im{i};
        %-- Pre-processing of image
        fprintf('Pre-processing of image... \n');
        %imPreProc = imProc;
        imPreProc{i} = im{i};
        if EnhanceContrast
            if low_in == 0 && high_in == 1 && low_out == 0 && high_out == 1
                imPreProc{i} = imadjust(imPreProc{i}); %Adjust image intensity values or colormap. Saturates the bottom 1% and the top 1% of all pixel values. Increases the contrast of the output image
            else
                imPreProc{i} = imadjust(imPreProc{i},[low_in high_in],[low_out high_out]);
            end
        end
        if SmoothingFilter
            imPreProc{i} = imgaussfilt(imPreProc{i},sigmaSmoothFilter,"FilterSize",sizeSmoothFilter); %Apply Gaussian smoothing filter (Gaussian smoothing kernel). Sigma is the standard deviation of the Gaussian distribution (increasing sigma, increases the smoothing/blur effect).
        end

        %-- Edge detection
        if EdgeDetMethodOpt == 1 || EdgeDetMethodOpt == 2 || EdgeDetMethodOpt == 3 %Sobel, Prewitt or Roberts
            if EdgeDetMethodOpt == 1
                fprintf('Sobel edge detection method... \n');
            elseif EdgeDetMethodOpt == 2
                fprintf('Prewitt edge detection method... \n');
            elseif EdgeDetMethodOpt == 3
                fprintf('Roberts edge detection method... \n');
            end
            edgeBinImage = edge(imPreProc{i},edgeDetMetName,edgeDetThreshold,edgeDetDirection); %Edge binary image
        elseif EdgeDetMethodOpt == 4 %LoG
            fprintf('Laplacian of Gaussian edge detection method... \n');
            edgeBinImage = edge(imPreProc{i},edgeDetMetName,edgeDetThreshold,edgeDetSigma); %Edge binary image
        elseif EdgeDetMethodOpt == 5 %Canny
            fprintf('Canny edge detection method... \n');
            edgeBinImage = edge(imPreProc{i},edgeDetMetName,edgeDetThreshold,edgeDetSigma); %Edge binary image
        elseif EdgeDetMethodOpt == 6 %Zerocross
            fprintf('Zero-crossing edge detection method... \n');
            if edgeDetFilter %FIlter modification and creation
                edgeBinImage = edge(imPreProc{i},edgeDetMetName,edgeDetThreshold,h); %Edge binary image
            else
                edgeBinImage = edge(imPreProc{i},edgeDetMetName,edgeDetThreshold); %Edge binary image
            end
        elseif EdgeDetMethodOpt == 7 %Partial Area Effect
            fprintf('Partial area effect edge detection method... \n');
            [edges, RI] = subpixelEdges(imPreProc{i}, edgeDetDifThreshold,'Order',edgeDetOrder,'SmoothingIter',edgeDetSmoothingIter); %find edge
            %- Compute the edge binary image
            %-- Fast but less precision process
            %edgeBinImage = zeros(size(imPreProc));
            %edgeBinImage(edges.position) = 1;
            %-- Computes a high resolution edge binary image
            edgeBinImage = subpixelImage(edges,size(imPreProc{i}),1);
        end

        %-- Post-processing (morphological operations)
        edgeBinPostProcImage = edgeBinImage;
        if MorphOpen
            fprintf('Running open morphological operations... \n');
            if defaultopenP
                openP = round(size(edgeBinPostProcImage,1)/2);
            end
            edgeBinPostProcImage = bwareaopen(edgeBinPostProcImage,openP,openConn);
        elseif MorphBridge
            fprintf('Running bridge morphological operations... \n');
            edgeBinPostProcImage = bwmorph(edgeBinPostProcImage,'bridge',bridgeN);
        elseif MorphClose
            fprintf('Running close morphological operations... \n');
            edgeBinPostProcImage = bwmorph(edgeBinPostProcImage,'close',closeN);
        elseif MorphThin
            fprintf('Running thin morphological operations... \n');
            edgeBinPostProcImage = bwmorph(edgeBinPostProcImage,'thin',thinN);
        end

        %-- Post-processing (identification of the longest edge)
        fprintf('Identification of the longest edge... \n');
        if EdgeDetMethodOpt ~= 7 || MorphOpen || MorphBridge || MorphClose || MorphThin
            B = bwboundaries(edgeBinPostProcImage);
            edgeCoords = B{1}; %Initialize edgeCoords array
            for j = 2:length(B)
                edgeCoords = cat(1,edgeCoords,B{j});
            end
            edges = EdgePixel;
            edges.x = edgeCoords(:,2);
            edges.y = edgeCoords(:,1);
            edges.position = sub2ind(size(imPreProc{i}),edges.y,edges.x); %Convert subscripts to linear indices
            edges.curv = LineCurvature2D([edges.x,edges.y]);
            N = LineNormals2D([edges.x,edges.y]);
            edges.nx = N(:,1);
            edges.ny = N(:,2);
            edges.i0 = zeros(length(edges.position),1);
            edges.i1 = zeros(length(edges.position),1);
        end
        longestedge = findlongestedge(edges,size(edgeBinPostProcImage),margin); % Selection of the longest edge

        %-- Post-processing (identification and sort of left and right edges)
        fprintf('Identification and sort of left and right edges... \n');
        [edgeL,edgeR] = leftrightedges(longestedge,jump_dist); % split edge into left and right and sort it

        %-- Baseline determination
        if BaselineDetOpt == 0  %- Manual baseline detection (from Akbari and Antonini 2021 code)
            fprintf('Manual baseline determination... \n');
            if i == 1 %First image, determination of reference points
                 while 1
                    nImageManDet = input('Enter an image number [1]: ');
                    if isempty(nImageManDet)
                        nImageManDet = 1;
                    end
                    if nImageManDet > nImages
                        fprintf('Invalid image number. Please enter a lower number. \n');
                    else
                        break;
                    end
                end
                [baseL,baseR,dropCenter] = ManualSelection(im{nImageManDet});
           end
           if yAverageLR %--- Use of average y left and right baseline (dropcenter(2))
                ybaseL = dropCenter(2);
                ybaseR = dropCenter(2);
            else %--- Use of  y left and right baseline separately
                ybaseL = baseL(2);
                ybaseR = baseR(2);
            end
            %-- Edge detection margin for CP search
            indexLMin = find(abs((edgeL.y-(ybaseL-(marginBaselineDet/2)))) == min(abs(edgeL.y-(ybaseL-(marginBaselineDet/2)))),1);
            indexLMax = find(abs((edgeL.y-(ybaseL+(marginBaselineDet/2)))) == min(abs(edgeL.y-(ybaseL+(marginBaselineDet/2)))),1);
            indexRMin = find(abs((edgeR.y-(ybaseR-(marginBaselineDet/2)))) == min(abs(edgeR.y-(ybaseR-(marginBaselineDet/2)))),1);
            indexRMax = find(abs((edgeR.y-(ybaseR+(marginBaselineDet/2)))) == min(abs(edgeR.y-(ybaseR+(marginBaselineDet/2)))),1);
            if isempty(indexLMin) || isempty(indexLMax) || isempty(indexRMin) || isempty(indexRMax)
                fprintf('Baseline cannot be identified! \n');
                fprintf('Try to modify yDistWeigth or marginBaselineDet. \n');
                fprintf('Suggestion1: Check that the entire edge is being identified. Review the edge detection method used. \n')
                %-- Definition of CPs
                x0L = edgeL.x(1);
                y0L = edgeL.y(1);
                x0R = edgeR.x(end);
                y0R = edgeR.y(end);
            else
                %-- CP search
                distL = sqrt((dropCenter(1)-edgeL.x(indexLMin:indexLMax)).^2 + ((ybaseL-edgeL.y(indexLMin:indexLMax)).^2).*yDistWeight); %Distance between edge points and drop baseline center
                indexL = find(distL == min(distL)) + indexLMin - 1;
                distR = sqrt((dropCenter(1)-edgeR.x(indexRMin:indexRMax)).^2 + ((ybaseR-edgeR.y(indexRMin:indexRMax)).^2).*yDistWeight); %Distance between edge points and drop baseline center
                indexR = find(distR == min(distR)) + indexRMin - 1;
                %-- Definition of CPs
                x0L = edgeL.x(indexL);
                y0L = edgeL.y(indexL);
                x0R = edgeR.x(indexR);
                y0R = edgeR.y(indexR);
            end
            if (isempty(x0L) || isempty(y0L)) && (isempty(x0R) || isempty(y0R))
                fprintf('Left and right baseline points cannot be identified! \n');
                fprintf('Try to modify yDistWeigth or marginBaselineDet. \n');
                fprintf('Suggestion1: Check that the entire edge is being identified. Review the edge detection method used. \n')
                %-- Definition of CPs
                x0L = edgeL.x(1);
                y0L = edgeL.y(1);
                x0R = edgeR.x(end);
                y0R = edgeR.y(end);
            elseif isempty(x0L) || isempty(y0L)
                fprintf('Left baseline point cannot be identified! \n');
                fprintf('Try to modify yDistWeigth or marginBaselineDet. \n');
                fprintf('Suggestion1: Check that the entire edge is being identified. Review the edge detection method used. \n')
                %-- Definition of CPs
                x0L = edgeL.x(1);
                y0L = edgeL.y(1);
            elseif isempty(x0R) || isempty(y0R)
                fprintf('Right baseline point cannot be identified! \n');
                fprintf('Try to modify yDistWeigth or marginBaselineDet. \n');
                fprintf('Suggestion1: Check that the entire edge is being identified. Review the edge detection method used. \n')
                %-- Definition of CPs
                x0R = edgeR.x(1);
                y0R = edgeR.y(1);
            end
        elseif BaselineDetOpt == 1 %- Automatic baseline detection for reflective drops (from Andersen et al 2017 code)
            fprintf('Automatic baseline detection for reflective drops... \n');
            [x0L,y0L,indexL] = findreflection([edgeL.x,edgeL.y],nBaselineDet,cutBaselineDet); %find reflection
            [x0R,y0R,indexR] = findreflection([edgeR.x,edgeR.y],nBaselineDet,cutBaselineDet);
        elseif BaselineDetOpt == 2  %- Automatic baseline detection for non-reflective drops (find baseline(z-plane)) (from Liu et al 2019 code)
            fprintf('Automatic baseline detection for non-reflective drops... \n');
            [p,pos0L,pos0R,pos0L_cut,pos0R_cut,index_left,index_right,index_cut_left,index_cut_right] = plane_mod([edgeL.x,edgeL.y],[edgeR.x,edgeR.y],percentBaselineDet);
            %- Use of percent
            x0L = pos0L_cut(:,1);
            y0L = pos0L_cut(:,2);
            x0R = pos0R_cut(:,1);
            y0R = pos0R_cut(:,2);
        elseif BaselineDetOpt == 3 %- Automatic baseline detection for non-reflective drops (find baseline(z-plane)) (from Akbari and Antonini 2021 code)(mask method)
            fprintf('Automatic baseline detection for reflective and non-reflective... \n');
            [I_New,AXxL,AXyL,AXxR,AXyR,Apex_x,Apex_y,Left_edgeTx,Left_edgeDx,Right_edgeTx,Right_edgeDx,Num_Down_L, Num_Up_L,Num_Down_R,...
                Num_Up_R,End_x_L,End_x_R] = ImageAnalysis2(imPreProc{i},edgeL,edgeR,fillBlankHorLines,marginMaskBaselineDet);
            %Automatic identification of contact points using mask method
            [b,alphac,alphamL,alphamR,x_CPL,y_CPL,x_CPR,y_CPR,indexCPL,indexCPR,Contact_angle_left,Contact_angle_right,Contact_angle,TiltAngle_S] ...
                = CircleMask2(imPreProc{i},I_New,maskSize,AXxL,AXyL,AXxR,AXyR,Apex_x,Apex_y,Left_edgeTx,Left_edgeDx,Right_edgeTx,Right_edgeDx,...
                Num_Down_L, Num_Up_L,Num_Down_R,Num_Up_R,End_x_L,End_x_R);
            if (isempty(x_CPL) || isempty(y_CPL)) && (isempty(x_CPR) || isempty(y_CPR))
                fprintf('Left and right baseline points cannot be identified! \n');
                fprintf('Suggestion1: Check that the entire edge is continuous. Review the edge detection method used. \n')
                %-- Definition of CPs
                x0L = edgeL.x(1);
                y0L = edgeL.y(1);
                x0R = edgeR.x(end);
                y0R = edgeR.y(end);
            elseif isempty(x_CPL) || isempty(y_CPL)
                fprintf('Left baseline point cannot be identified! \n');
                fprintf('Suggestion1: Check that the entire edge is continuous. Review the edge detection method used. \n')
                %-- Definition of CPs
                x0L = edgeL.x(1);
                y0L = edgeL.y(1);
            elseif isempty(x_CPR) || isempty(y_CPR)
                fprintf('Right baseline point cannot be identified! \n');
                fprintf('Suggestion1: Check that the entire edge is continuous. Review the edge detection method used. \n')
                %-- Definition of CPs
                x0R = edgeR.x(1);
                y0R = edgeR.y(1);
            else
                y0L = x_CPL;
                x0L = y_CPL;
                y0R = x_CPR;
                x0R = y_CPR;
            end
        end

        %-- Store edge detection image
        if mod(i,plotIntEdgeDetImages) == 0 || i == 1
            fprintf('Storing edge detection image... \n');
            figureName = strcat('EdgeDetImage_',fileNames{i});
            figure('Name',figureName);
            %-- Show original image
            imshow(im{i})
            hold on
            %-- Plot baseline
            t = linspace(-3,3);
            plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r--','LineWidth',2,'DisplayName','Baseline')
            %-- Plot identified edges
            plot(edgeL.x,edgeL.y,':r','LineWidth',2,'DisplayName','Left edge')
            plot(edgeR.x,edgeR.y,':b','LineWidth',2,'DisplayName','Right edge')
            %-- Display contact points
            plot(x0L,y0L,'yx','MarkerSize',10,'LineWidth',2,'DisplayName','Left CP')
            plot(x0R,y0R,'yx','MarkerSize',10,'LineWidth',2,'DisplayName','Right CP')
            %-- Display legend
            legend('Location','north','Orientation','horizontal','FontSize',8)
            %-- Save figure
            fullfileFigureName = fullfile(plotEdgeDetectionPath,figureName);
            saveas(gcf,fullfileFigureName,'png')
            close(gcf)
        end

        %-- Check baseline tilt
        %{
        if abs(BaselineTilt(i)) >= limBaselineTilt
            fprintf('Baseline tilt > %.1f°. \n',limBaselineTilt);
            fprintf('Resulted image should be discarded! \n');
            contExcTilt = contExcTilt + 1; % Count images with an excess of abseline tilt
            nExcTiltImages(contExcTilt) = i;
        else
            fprintf('Baseline tilt = %.2f° \n',BaselineTilt(i));
        end
        %}
        %-- Store data
        BaselineTilt(i) = atan2(y0R-y0L,x0R-x0L) * 180/pi; %Store baseline tilt
        fprintf('Baseline tilt = %.2f° \n',BaselineTilt(i));
        EdgeCell{i,1} = edgeL; % store left detected edge
        EdgeCell{i,2} = edgeR; % store right detected edge
        BaseVec(i,:) = [x0L,y0L,x0R,y0R]; % store baseline coordinates
    end
    fprintf('Edge detection and baseline determination completed! \n');

    %-- Analyzing baseline tilt and separate images to fit and determine CA
    fprintf('Analyzing baseline tilt...\n')
    while 1
        %-- Report
        contExcTilt = 0;
        nExcTiltImages = 0;
        listExcTiltImages = "";
        for j = 1:length(BaselineTilt)
            if abs(BaselineTilt(j)) >= limBaselineTilt
                contExcTilt = contExcTilt + 1; % Count images with an excess of baseline tilt
                nExcTiltImages(contExcTilt) = j; % Store the image numbers that exceeded the baseline tilt limit  
            end
        end
        fprintf('- Report: \n');
        fprintf('Total analyzed images: %i \n',nImages);
        fprintf('Maximum allowed limit of baseline tilt: %.2f° \n',limBaselineTilt)
        fprintf('Number of analyzed images that should be discarded: %i \n',contExcTilt);
        percDiscIm = contExcTilt*100/nImages; % percent of discarded images
        fprintf('Percent of discarded images: %.2f%% \n',percDiscIm);
        fprintf('List of images that should be discarded: ');
        for j = 1:length(nExcTiltImages)
            if j == length(nExcTiltImages)
                fprintf('%i. \n',nExcTiltImages(j));
                listExcTiltImages = strcat(listExcTiltImages,string(nExcTiltImages(j)));
            else
                fprintf('%i, ',nExcTiltImages(j));
                listExcTiltImages = strcat(listExcTiltImages,string(nExcTiltImages(j)),",");
            end
        end
        ansBaselineTilt = input('Do you want to change baseline tilt limit and update the report (Y/N) [N]? ','s');
        if isempty(ansBaselineTilt)
            ansBaselineTilt = 'N';
        end
        if ansBaselineTilt == 'N'
            break;
        else
            limBaselineTilt = input('Max. limit of baseline tilt (°) [1]: ');
            if isempty(limBaselineTilt)
                limBaselineTilt = 1;
            end
            %listExcTiltImages = "";
        end
    end
    %-- Ask if the user wants to run again edge detection routine
    ansEndEdgeDetection = input('Do you like to redo all the edge detection routine (Y/N) [N]? ','s');
    if isempty(ansEndEdgeDetection)
        ansEndEdgeDetection = 'N';
    end
    if ansEndEdgeDetection == 'N'
        break;
    else
        %rmdir(plotEdgeDetectionPath)
    end
end
%- Saving report
fileOptName = strcat(outputName,'_EdgeDetReport',outputFmt);
fullfileOptName = fullfile(outputPath,fileOptName);
varnames = {'TotalNumberOfImages','LimitbaselineTilt[°]','NumberOfImagesDiscarded','PercentDiscardedImages','ListDiscardedImages'}; %Nome das variáveis no arquivo txt
if isempty(edgeDetThreshold)
    edgeDetThresholdValue = 0;
end
T = table(nImages,round(limBaselineTilt,4),contExcTilt,string(strcat(num2str(round(percDiscIm,2)),"%")),listExcTiltImages,'VariableNames',varnames);
writetable(rows2vars(T),fullfileOptName,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file

fprintf('Edge detection and baseline determination routine completed! \n');
%%
%---------------------------------------------------------------------------
%        5.       Average baselines found in all frames
%---------------------------------------------------------------------------
fprintf('---------------------------- 5. AVERAGE BASELINE -------------------------- \n');
ansAverageBaseVec = input('Do you like to average baseline vector found in all frames (Y/N) [N]? ','s');
if isempty(ansAverageBaseVec)
    ansAverageBaseVec = 'N';
end
if ansAverageBaseVec == 'Y'
    UseAverageBaseVec = 1;
else
    UseAverageBaseVec = 0;
end
if UseAverageBaseVec % Find the average baseline using all baseline coordinates
    xBase = [BaseVec(:,1);BaseVec(:,3)]; %baseline x coordinates
    yBase = [BaseVec(:,2);BaseVec(:,4)]; %baseline y coordinates
    B = AverageBaseline(xBase,yBase);
    %AverageBaseVec = bsxfun(@times,ones(nImages,4),[0,B(1),size(im,2),size(im,2)*B(2)+B(1)]);
    fprintf('Average baseline completed! \n');
end
%%
%---------------------------------------------------------------------------
%                   6.       Fit contact angles
%---------------------------------------------------------------------------
fprintf('------------------------ 6. FIT CONTACT ANGLES --------------------------- \n');
while 1
    % Predefine variables prior to loop
    TLvec = NaN(nImages,4);
    CircleCAS = NaN(nImages,2);
    PolyDropenCAS = NaN(nImages,2);
    PolyCAS = NaN(nImages,2);
    EllipseCAS = NaN(nImages,2);
    MaskCAS = NaN(nImages,2);
    contactDist = NaN(nImages,1);
    dispCP = NaN(nImages,2);
    dispCPmm = NaN(nImages,2);
    dispCPStop = NaN(nImages,2);
    dispCPStopmm = NaN(nImages,2);
    vel = NaN(nImages,2); %Velocity of the triple (contact) points
    volEstWasher = NaN(nImages,1); %Drop volume estimated from disk method
    volEstSph = NaN(nImages,1); %Drop volume estimated from spherical cap hypothesis
    xVel = NaN(nImages,2); %x coordinates for the calculation of triple line velocity
    yVel = NaN(nImages,2); %y coordinates for the calculation of triple line velocity
    volVel = NaN(nImages,1); %Experimental volume for the calculation of triple line velocity

    % Fitting settings
    ansFittingSettings = input('Do you like to change fitting options (Y/N) [N]? ','s');
    if isempty(ansFittingSettings)
        ansFittingSettings = 'N';
    end
    if ansFittingSettings == 'Y'
        [UseCircFit,n_circfit,PlotCircleImages,UsePolyFit,n_polyfit,poly_ord,PlotPolynomialImages,UsePolyFitDropen,n_polyfitDropen,poly_ordDropen,...
            PlotPolynomialDropenImages,UseEllipseFit,PlotEllipticImages,UseMaskMethod,fillBlankHorLines,marginMaskBaselineDet,maskSize,PlotMaskImages] = SettingsFittingCA(UseCircFit,...
            n_circfit,PlotCircleImages,UsePolyFit,n_polyfit,poly_ord,PlotPolynomialImages,UsePolyFitDropen,n_polyfitDropen,poly_ordDropen,PlotPolynomialDropenImages...
            ,UseEllipseFit,PlotEllipticImages,UseMaskMethod,fillBlankHorLines,marginMaskBaselineDet,maskSize,PlotMaskImages);
    end
    %- Plot interval
    plotIntFittingImages = input('Plot interval of fitting images [5]: ');
    if isempty(plotIntFittingImages)
        plotIntFittingImages = 5;
    end
    %- Create folders to export plotting images
    if PlotCircleImages
        dirCircFitting = 'CircleFitting';
        mkdir(outputPath,dirCircFitting);
        plotCircFittingPath = fullfile(outputPath,dirCircFitting);
    end
    if PlotPolynomialImages
        dirPolyFitting = 'PolyFitting';
        mkdir(outputPath,dirPolyFitting);
        plotPolyFittingPath = fullfile(outputPath,dirPolyFitting);
    end
    if PlotPolynomialDropenImages
        dirPolyDropenFitting = 'PolyDropenFitting';
        mkdir(outputPath,dirPolyDropenFitting);
        plotPolyDropenFittingPath = fullfile(outputPath,dirPolyDropenFitting);
    end
    if PlotEllipticImages
        dirEllipFitting = 'EllipseFitting';
        mkdir(outputPath,dirEllipFitting);
        plotEllipFittingPath = fullfile(outputPath,dirEllipFitting);
    end
    if PlotMaskImages
        dirMaskMethod = 'MaskMethod';
        mkdir(outputPath,dirMaskMethod);
        plotMaskMethodPath = fullfile(outputPath,dirMaskMethod);
    end

    % Fitting data to estimate contact diameter, drop volume and CA
    fprintf('Fitting data to calculate CA... \n');
    for i = 1:nImages
        fprintf('Analyzing image %d of %d... \n',i,nImages);
        if any(i == nExcTiltImages) %Discard image
            fprintf('Limit baseline tilt exceeded! \n')
            fprintf('Resulted image discarded! \n')
        else
            %- Get left and right detected edges
            edgeL = EdgeCell{i,1};
            edgeR = EdgeCell{i,2};

            %- Use of average baseline vector
            if UseAverageBaseVec
                %v = num2cell(AverageBaseVec(i,:));
                %[x0L,y0L,x0R,y0R] = v{:};
                x0L = BaseVec(i,1);
                y0L = B(1);
                x0R = BaseVec(i,3);
                y0R = B(1);
            else
                %[x0L,y0L,x0R,y0R] = BaseVec(i,:);
                x0L = BaseVec(i,1);
                y0L = BaseVec(i,2);
                x0R = BaseVec(i,3);
                y0R = BaseVec(i,4);
            end

            %{
        %- Calculate baseline tilt
        fprintf('Calculating baseline tilt... \n');
        if y0L ~= y0R
            TiltAngle_S = atan2(y0R-y0L,x0R-x0L) * 180/pi;
        end
            %}
            %- Defining x coordinates and tilting angle arrays for the calculation of triple line velocity
            xVel(i,1) = x0L; %Left x coordinate for the calculation of the triple line velocity
            xVel(i,2) = x0R; %Right x coordinate for the calculation of the triple line velocity
            yVel(i,1) = y0L; %Left y coordinate for the calculation of the triple line velocity
            yVel(i,2) = y0R; %Right y coordinate for the calculation of the triple line velocity
            volVel(i,1) = volExp(i); %Experimental volume for the calculation of triple line velocity

            % Calculate left and right contact points displacement
            fprintf('Calculating left and right contact points displacement... \n');
            dispCP(i,1) = abs(BaseVec(i,1) - BaseVec(1,1)); %Displacement of the left contact point in pixels
            dispCP(i,2) = abs(BaseVec(i,3) - BaseVec(1,3)); %Displacement of the right contact point in pixels
            dispCPmm(i,1) = dispCP(i,1)/scale; %Displacement of the left contact point in mm
            dispCPmm(i,2) = dispCP(i,2)/scale; %Displacement of the right contact point in mm
            if mod(i,nImagesStops) == 1 || nImagesStops == 1 %Update displacement reference
                dispCPStopRef(1) = BaseVec(i,1);
                dispCPStopRef(2) = BaseVec(i,3);
            end
            dispCPStop(i,1) = abs(BaseVec(i,1) - dispCPStopRef(1));
            dispCPStop(i,2) = abs(BaseVec(i,3) - dispCPStopRef(2));
            dispCPStopmm(i,1) = dispCPStop(i,1)/scale;
            dispCPStopmm(i,2) = dispCPStop(i,2)/scale;
            fprintf('Left CP displacement = %.2f mm \n', dispCPStopmm(i,1));
            fprintf('Right CP displacement = %.2f mm \n', dispCPStopmm(i,2));

            %- Calculate contact distance (diameter) in pixels
            fprintf('Estimating contact diameter... \n');
            contactDist(i) = sqrt((x0R-x0L)^2+(y0R-y0L)^2);
            fprintf('Contact diameter = %.2f pixels \n', contactDist(i))

            %- Circle Fitting (from DropenV01.m - Akbari and Antonini 2021)
            if UseCircFit
                fprintf('Using circle fitting to calculate CA... \n');
                %CircleData = circlefit([AXxL,AXyL],[AXxR,AXyR],[x0L,y0L],[x0R,y0R],n_circfit);
                CircleData = circlefit([edgeL.y,edgeL.x],[edgeR.y,edgeR.x],[x0L,y0L],[x0R,y0R],n_circfit);
                %- Store obtained data
                CircleCAS(i,1) = CircleData.CAL; % Store left CA obtained from circle fitting
                CircleCAS(i,2) = CircleData.CAR; % Store right CA obtained from circle fitting
                fprintf('Circle fitting CAL = %.2f° \n', CircleCAS(i,1));
                fprintf('Circle fitting CAR = %.2f° \n', CircleCAS(i,2));
            end
            %- Polynomial Fitting (from DropenV01.m - Akbari and Antonini 2021)
            if UsePolyFitDropen
                fprintf('Using polynomial fitting (Dropen) to calculate CA... \n');
                %PolyDataDropen = polynomialfit_Dropen([AXxL,AXyL],[AXxR,AXyR],[x0L,y0L],[x0R,y0R],n_polyfitDropen,poly_ordDropen);
                PolyDataDropen = polynomialfit_Dropen([edgeL.y,edgeL.x],[edgeR.y,edgeR.x],[x0L,y0L],[x0R,y0R],n_polyfitDropen,poly_ordDropen);
                %- Store obtained data
                PolyDropenCAS(i,1) = PolyDataDropen.CAL; % Store left CA obtained from polynomial fitting
                PolyDropenCAS(i,2) = PolyDataDropen.CAR; % Store right CA obtained from polynomial fitting
                fprintf('Polynomial fitting (Dropen) CAL = %.2f° \n', PolyDropenCAS(i,1));
                fprintf('Polynomial fitting (Dropen) CAR = %.2f° \n', PolyDropenCAS(i,2));
            end
            %- Polynomial fitting
            if UsePolyFit
                fprintf('Using polynomial fitting to calculate CA... \n');
                %PolyData = polynomialfit(edgeL,edgeR,[x0L,y0L],[x0R,y0R],n_polyfit,poly_ord); % Fit data to polynomial
                PolyData = polynomialfit(edgeL,edgeR,[x0L,y0L],[x0R,y0R],n_polyfit,poly_ord); % Fit data to polynomial
                %- Store obtained data
                if isempty(PolyData.CAL) || isempty(PolyData.CAR)
                    PolyData.CAL = NaN;
                    PolyData.CAR = NaN;
                end
                PolyCAS(i,1) = PolyData.CAL; % Store left CA obtained from polynomial fitting
                PolyCAS(i,2) = PolyData.CAR; % Store right CA obtained from polynomial fitting
                fprintf('Polynomial fitting CAL = %.2f° \n', PolyCAS(i,1));
                fprintf('Polynomial fitting CAR = %.2f° \n', PolyCAS(i,2));
            end
            %- Double sided elipptical fitting
            if UseEllipseFit
                fprintf('Using double-sided elliptical method to calculate CA... \n');
                EllipseData = EllipticFit(edgeL,edgeR,[x0L,y0L],[x0R,y0R]); % Fit data to ellipse
                %- Store obtained data
                EllipseCAS(i,1) = EllipseData.CAL; % Store left CA obtained from ellipse fitting
                EllipseCAS(i,2) = EllipseData.CAR; % Store right CA obtained from ellipse fitting
                TLvec(i,1:2)=EllipseData.TLL;
                TLvec(i,3:4)=EllipseData.TLR;
                fprintf('Ellipse fitting CAL = %.2f° \n', EllipseCAS(i,1));
                fprintf('Ellipse fitting CAR = %.2f° \n', EllipseCAS(i,2));
            end
            %- Mask method
            if UseMaskMethod
                fprintf('Using mask method to calculate CA... \n');
                [I_New,AXxL,AXyL,AXxR,AXyR,Apex_x,Apex_y,Left_edgeTx,Left_edgeDx,Right_edgeTx,Right_edgeDx,Num_Down_L, Num_Up_L,Num_Down_R,...
                    Num_Up_R,End_x_L,End_x_R] = ImageAnalysis2(imPreProc{i},edgeL,edgeR,fillBlankHorLines,marginMaskBaselineDet);
                [b,alphac,alphamL,alphamR,x_CPL,y_CPL,x_CPR,y_CPR,indexCPL,indexCPR,maskCAL,maskCAR,maskCA,TiltAngle_S] ...
                    = CircleMask2(imPreProc{i},I_New,maskSize,AXxL,AXyL,AXxR,AXyR,Apex_x,Apex_y,Left_edgeTx,Left_edgeDx,Right_edgeTx,Right_edgeDx,...
                    Num_Down_L, Num_Up_L,Num_Down_R,Num_Up_R,End_x_L,End_x_R);
                if BaselineDetOpt ~= 3
                    %Find index closest to the baseline points
                    distL = sqrt((x0L-alphamL(2,:)).^2 + (y0L-alphamL(1,:)).^2);
                    indexCPL = find(distL == min(distL),1);
                    distR = sqrt((x0R-alphamR(2,:)).^2 + (y0R-alphamR(1,:)).^2);
                    indexCPR = find(distR == min(distR),1);
                    %Calculation of contact angles
                    maskCAL = alphamL(5,indexCPL) + BaselineTilt(i);
                    maskCAR = alphamR(5,indexCPR) - BaselineTilt(i);
                end
                %- Store obtained data
                MaskCAS(i,1) = maskCAL; % Store left CA obtained from ellipse fitting
                MaskCAS(i,2) = maskCAR; % Store right CA obtained from ellipse fitting
                fprintf('Mask method CAL = %.2f° \n', MaskCAS(i,1));
                fprintf('Mask method CAR = %.2f° \n', MaskCAS(i,2));
            end

            %- Estimate drop volume
            fprintf('Estimating drop volume... \n');
            %-- Using disk method
            volEstWasher(i,1) = dropVolumeWasher(edgeL,edgeR,[x0L y0L],[x0R y0R],scale); %Drop volume in uL
            fprintf('Estimated volume (disk method) = %.2f \x3bcL \n', volEstWasher(i,1));
            %-- Considering a spherical drop
            if ~isempty(PolyData.TLL) && ~isempty(PolyData.TLR) && ~isempty(PolyData.CAL) && ~isempty(PolyData.CAR) %Check that polynomial fitting converged to results
                avgCA = (PolyCAS(i,1) + PolyCAS(i,2))/2; %Average CA considering polynomial fitting data
            else
                avgCA = (PolyDropenCAS(i,1) +  PolyDropenCAS(i,2))/2; %Average CA considering polynomial Dropen fitting data
            end
            volEstSph(i,1) = dropVolumeSph([x0L y0L],[x0R y0R],avgCA,scale); %Drop volume estimated by spherical cap hypothesis in uL
            fprintf('Estimated volume (spherical assumption) = %.2f \x3bcL \n', volEstSph(i,1));
            %}

            %- Plot fitted data
            %-- Plot an example image for every stop
            %if mod(i,nImagesStops) == 0
            if mod(i,plotIntFittingImages) == 0 || i == 1
                %-- Plot circle fitting data (from DropenV01.m - Akbari and Antonini 2021)
                if PlotCircleImages
                    fprintf('Plotting and exporting circle fitting image... \n');
                    figureName = strcat('CircleFit_',fileNames{i});
                    figure('Name',figureName);
                    %-- Show original image
                    imshow(im{i})
                    hold on
                    %-- Plot baseline
                    t = linspace(-3,3);
                    plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r--','LineWidth',2,'DisplayName','Baseline')
                    %-- Plot fitted edge
                    plot(CircleData.EvalCircL(:,1),CircleData.EvalCircL(:,2),':r','LineWidth',2,'DisplayName','Left circle fitted edge')
                    plot(CircleData.EvalCircR(:,1),CircleData.EvalCircR(:,2),':b','LineWidth',2,'DisplayName','Right circle fitted edge')
                    %-- Plot fitted data closer to contact points
                    radius = 40;
                    tilt = atand((y0R-y0L)/(x0R-x0L));
                    plot([CircleData.TLL(2),CircleData.TLL(2)+radius*cosd(CircleData.CAL-tilt)],[CircleData.TLL(1),CircleData.TLL(1)-radius*sind(CircleData.CAL-tilt)],'LineStyle','--','LineWidth', 2,'color','g','DisplayName','Fitted data closer to left CP')
                    plot([CircleData.TLR(2),CircleData.TLR(2)-radius*cosd(CircleData.CAR+tilt)],[CircleData.TLR(1),CircleData.TLR(1)-radius*sind(CircleData.CAR+tilt)],'LineStyle','--','LineWidth', 2,'color','c','DisplayName','Fitted data closer to right CP')
                    %-- Display contact points
                    plot(x0L,y0L,'mx','MarkerSize',10,'LineWidth',2,'DisplayName','Left CP')
                    plot(x0R,y0R,'yx','MarkerSize',10,'LineWidth',2,'DisplayName','Right CP')
                    %-- Display legend
                    legend('Location','northeastoutside','Orientation','vertical','FontSize',7)
                    %-- Save figure
                    %fullfileFigureName = fullfile(outputPath,figureName);
                    fullfileFigureName = fullfile(plotCircFittingPath,figureName);
                    saveas(gcf,fullfileFigureName,'png')
                    close(gcf)
                end
                %-- Plot polynomial fitting data (from DropenV01.m - Akbari and Antonini 2021)
                if PlotPolynomialDropenImages
                    fprintf('Plotting and exporting polynomial fitting (Dropen) image... \n');
                    figureName = strcat('PolyFitDropen_',fileNames{i});
                    figure('Name',figureName);
                    %-- Show original image
                    imshow(im{i})
                    hold on
                    %-- Plot baseline
                    t = linspace(-3,3);
                    plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r--','LineWidth',2,'DisplayName','Baseline')
                    %-- Plot fitted edge
                    plot(PolyDataDropen.EvalPolyL(:,1),PolyDataDropen.EvalPolyL(:,2),':r','LineWidth',2,'DisplayName','Left polynomial fitted edge')
                    plot(PolyDataDropen.EvalPolyR(:,1),PolyDataDropen.EvalPolyR(:,2),':b','LineWidth',2,'DisplayName','Right polynomial fitted edge')
                    %-- Plot fitted data closer to contact points
                    radius = 40;
                    tilt = atand((y0R-y0L)/(x0R-x0L));
                    plot([PolyDataDropen.TLL(2),PolyDataDropen.TLL(2)+radius*cosd(PolyDataDropen.CAL-tilt)],[PolyDataDropen.TLL(1),PolyDataDropen.TLL(1)-radius*sind(PolyDataDropen.CAL-tilt)],'LineStyle','--','LineWidth', 2,'color','g','DisplayName','Fitted data closer to left CP')
                    plot([PolyDataDropen.TLR(2),PolyDataDropen.TLR(2)-radius*cosd(PolyDataDropen.CAR+tilt)],[PolyDataDropen.TLR(1),PolyDataDropen.TLR(1)-radius*sind(PolyDataDropen.CAR+tilt)],'LineStyle','--','LineWidth', 2,'color','c','DisplayName','Fitted data closer to right CP')
                    %-- Display contact points
                    plot(x0L,y0L,'mx','MarkerSize',10,'LineWidth',2,'DisplayName','Left CP')
                    plot(x0R,y0R,'yx','MarkerSize',10,'LineWidth',2,'DisplayName','Right CP')
                    %-- Display legend
                    legend('Location','northeastoutside','Orientation','vertical','FontSize',7)
                    %-- Save figure
                    %fullfileFigureName = fullfile(outputPath,figureName);
                    fullfileFigureName = fullfile(plotPolyDropenFittingPath,figureName);
                    saveas(gcf,fullfileFigureName,'png')
                    close(gcf)
                end
                %-- Plot polynomial fitting data
                if PlotPolynomialImages
                    if ~isempty(PolyData.TLL) && ~isempty(PolyData.TLR) && ~isempty(PolyData.CAL) && ~isempty(PolyData.CAR) %Check that polynomial fitting converged to results
                        fprintf('Plotting and exporting polynomial fitting image... \n');
                        figureName = strcat('PolyFit_',fileNames{i});
                        figure('Name',figureName);
                        %-- Show original image
                        imshow(im{i})
                        hold on
                        %-- Plot baseline
                        t = linspace(-3,3);
                        plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r--','LineWidth',2,'DisplayName','Baseline')
                        %-- Plot fitted edge
                        plot(PolyData.EvalPolyL(:,1),PolyData.EvalPolyL(:,2),':r','LineWidth',2,'DisplayName','Left polynomial fitted edge')
                        plot(PolyData.EvalPolyR(:,1),PolyData.EvalPolyR(:,2),':b','LineWidth',2,'DisplayName','Right polynomial fitted edge')
                        %-- Plot fitted data closer to contact points
                        radius = 40;
                        tilt = atand((y0R-y0L)/(x0R-x0L));
                        plot([PolyData.TLL(1),PolyData.TLL(1)+radius*cosd(PolyData.CAL-tilt)],[PolyData.TLL(2),PolyData.TLL(2)-radius*sind(PolyData.CAL-tilt)],'LineStyle','--','LineWidth', 2,'color','g','DisplayName','Fitted data closer to left CP')
                        plot([PolyData.TLR(1),PolyData.TLR(1)-radius*cosd(PolyData.CAR+tilt)],[PolyData.TLR(2),PolyData.TLR(2)-radius*sind(PolyData.CAR+tilt)],'LineStyle','--','LineWidth', 2,'color','c','DisplayName','Fitted data closer to right CP')
                        %-- Display contact points
                        plot(x0L,y0L,'mx','MarkerSize',10,'LineWidth',2,'DisplayName','Left CP')
                        plot(x0R,y0R,'yx','MarkerSize',10,'LineWidth',2,'DisplayName','Right CP')
                        %-- Display legend
                        legend('Location','northeastoutside','Orientation','vertical','FontSize',7)
                        %-- Save figure
                        %fullfileFigureName = fullfile(outputPath,figureName);
                        fullfileFigureName = fullfile(plotPolyFittingPath,figureName);
                        saveas(gcf,fullfileFigureName,'png')
                        close(gcf)
                    end
                end
                %-- Plot ellipse fitting data
                if PlotEllipticImages
                    if ~isnan(max(EllipseData.TLL)) && ~isnan(max(EllipseData.TLR)) && ~isnan(EllipseData.CAL) && ~isnan(EllipseData.CAR) %Check that locate an ellipse
                        fprintf('Plotting and exporting double sided elliptical fitting image... \n');
                        figureName = strcat('EllipseFit_',fileNames{i});
                        figure('Name',figureName);
                        %-- Show original image
                        imshow(im{i})
                        hold on
                        %-- Plot baseline
                        t = linspace(-3,3);
                        plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r--','LineWidth',2,'DisplayName','Baseline')
                        %-- Plot fitted edge
                        tilt=atand((y0R-y0L)/(x0R-x0L));
                        tvec=linspace(0,2*pi)';
                        ellipserim=@(e,t) ones(size(t))*[e.X0_in,e.Y0_in]+cos(t)*[cos(-e.phi),sin(-e.phi)]*e.a+sin(t)*[-sin(-e.phi),cos(-e.phi)]*e.b;
                        EL=ellipserim(EllipseData.ellipse{1},tvec);
                        ER=ellipserim(EllipseData.ellipse{2},tvec);
                        plot(EL(:,1),EL(:,2),':r','LineWidth',2,'DisplayName','Left ellipse fitted edge');
                        plot(ER(:,1),ER(:,2),':b','LineWidth',2,'DisplayName','Right ellipse fitted edge');
                        %-- Plot fitted data closer to contact points
                        radius = 40;
                        plot([EllipseData.TLL(1),EllipseData.TLL(1)+radius*cosd(EllipseData.CAL-tilt)],[EllipseData.TLL(2),EllipseData.TLL(2)-radius*sind(EllipseData.CAL-tilt)],'LineStyle','--','LineWidth', 2,'color','g','DisplayName','Fitted data closer to left CP')
                        plot([EllipseData.TLR(1),EllipseData.TLR(1)-radius*cosd(EllipseData.CAR+tilt)],[EllipseData.TLR(2),EllipseData.TLR(2)-radius*sind(EllipseData.CAR+tilt)],'LineStyle','--','LineWidth', 2,'color','c','DisplayName','Fitted data closer to right CP')
                        %-- Display contact points
                        plot(x0L,y0L,'mx','MarkerSize',10,'LineWidth',2,'DisplayName','Left CP')
                        plot(x0R,y0R,'yx','MarkerSize',10,'LineWidth',2,'DisplayName','Right CP')
                        %-- Display legend
                        legend('Location','northeastoutside','Orientation','vertical','FontSize',7)
                        %-- Save figure
                        %fullfileFigureName = fullfile(outputPath,figureName);
                        fullfileFigureName = fullfile(plotEllipFittingPath,figureName);
                        saveas(gcf,fullfileFigureName,'png')
                        close(gcf)
                    end
                end
                %-- Plot mask method data
                if PlotMaskImages
                    %(i) Drop image
                    fprintf('Plotting and exporting mask method images... \n');
                    figureName = strcat('MaskMethod_',fileNames{i});
                    figure('Name',figureName);
                    %-- Show original image
                    imshow(im{i})
                    hold on
                    %-- Plot baseline
                    t = linspace(-3,3);
                    plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r--','LineWidth',2,'DisplayName','Baseline')
                    %-- Plot detected edges and contact points
                    plot(AXyL,AXxL,':r','LineWidth',2,'DisplayName','Left detected edge')
                    plot(AXyR,AXxR,':b','LineWidth',2,'DisplayName','Right detected edge')
                    %-- Plot fitted data closer to contact points
                    radius = 40;
                    tilt=atand((y0R-y0L)/(x0R-x0L));
                    plot([x0L,x0L+radius*cosd(maskCAL-tilt)],[y0L,y0L-radius*sind(maskCAL-tilt)],'LineStyle','--','LineWidth', 2,'color','g','DisplayName','Fitted data closer to left CP')
                    plot([x0R,x0R-radius*cosd(maskCAR+tilt)],[y0R,y0R-radius*sind(maskCAR+tilt)],'LineStyle','--','LineWidth', 2,'color','c','DisplayName','Fitted data closer to right CP')
                    %-- Display contact points
                    plot(x0L,y0L,'mx','MarkerSize',10,'LineWidth',2,'DisplayName','Left CP')
                    plot(x0R,y0R,'yx','MarkerSize',10,'LineWidth',2,'DisplayName','Right CP')
                    %-- Display legend
                    legend('Location','northeastoutside','Orientation','vertical','FontSize',7)
                    %-- Save figure
                    %fullfileFigureName = fullfile(outputPath,figureName);
                    fullfileFigureName = fullfile(plotMaskMethodPath,figureName);
                    saveas(gcf,fullfileFigureName,'png')
                    close(gcf)

                    %(ii) Arc length x local slope
                    figureName = strcat('MaskMethod_SxAlpha_',fileNames{i});
                    figure('Name',figureName);
                    %- Plot s x alpha
                    plot(alphamL(3,:),alphamL(5,:),'r--*','DisplayName','left edge');
                    hold on
                    plot(alphamR(3,:),alphamR(5,:),'b--*','DisplayName','right edge');
                    plot(alphamL(3,indexCPL),alphamL(5,indexCPL),'o','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerSize',18,'DisplayName','left CP'); %left contact point
                    plot(alphamR(3,indexCPR),alphamR(5,indexCPR),'o','MarkerEdgeColor',[0 0.4470 0.7410],'MarkerSize',18,'DisplayName','right CP'); %right contact point
                    xlabel('s');
                    ylabel('\alpha [°]')
                    ylim([min(min(alphamL(5,:)),min(alphamR(5,:))) (max(max(alphamL(5,:)),max(alphamR(5,:)))+5)])
                    legend('Location','southoutside','Orientation','horizontal')
                    %-- Save figure
                    %fullfileFigureName = fullfile(outputPath,figureName);
                    fullfileFigureName = fullfile(plotMaskMethodPath,figureName);
                    saveas(gcf,fullfileFigureName,'png')
                    close(gcf)
                end
            end
        end
    end
    %-- Ask if the user wants to run again fitting CA routine
    ansEndFittingCA = input('Do you like to redo all the fitting CA routine (Y/N) [N]? ','s');
    if isempty(ansEndFittingCA)
        ansEndFittingCA = 'N';
    end
    if ansEndFittingCA == 'N'
        break;
    else
        %{
        %- Delete folders created to export plotting images
        if PlotCircleImages
            rmdir(plotCircFittingPath)
        end
        if PlotPolynomialImages
            rmdir(plotPolyFittingPath)
        end
        if PlotPolynomialDropenImages
            rmdir(plotPolyDropenFittingPath)
        end
        if PlotEllipticImages
            rmdir(plotEllipFittingPath)
        end
        if PlotMaskImages
            rmdir(plotMaskMethodPath)
        end
        %}
    end
end
%- Calculate triple line velocity [um/°]
if nImages > 2
    fprintf('Calculating triple line velocity for all images... \n');
    for i = 2:nImages-1
        %- Three-point linear regression
        if (abs(xVel(i,1) - xVel(i-1,1)) == 0) || (abs(xVel(i,2) - xVel(i-1,2)) == 0)
            fL = 1;
            fR = 1;
        else
            fL = -(xVel(i,1) - xVel(i-1,1))/abs((xVel(i,1) - xVel(i-1,1))); %f for the left side(+1 for advancing movement and -1 for receding movement)
            fR = (xVel(i,2) - xVel(i-1,2))/abs((xVel(i,2) - xVel(i-1,2))); %f for the right side(+1 for advancing movement and -1 for receding movement)
        end
        %fL = -(xVel(i,1) - xVel(i-1,1))/abs((xVel(i,1) - xVel(i-1,1))); %f for the left side(+1 for advancing movement and -1 for receding movement)
        dxdtiltL = ((xVel(i-1,1)*volVel(i-1) + xVel(i,1)*volVel(i) + xVel(i+1,1)*volVel(i+1))*3 - (volVel(i-1) + volVel(i) + volVel(i+1))*(xVel(i-1,1) + xVel(i,1) + xVel(i+1,1))) / ((volVel(i-1)^2 + volVel(i)^2 + volVel(i+1)^2)*3 - (volVel(i-1) + volVel(i) + volVel(i+1))^2);
        dydtiltL = ((yVel(i-1,1)*volVel(i-1) + yVel(i,1)*volVel(i) + yVel(i+1,1)*volVel(i+1))*3 - (volVel(i-1) + volVel(i) + volVel(i+1))*(yVel(i-1,1) + yVel(i,1) + yVel(i+1,1))) / ((volVel(i-1)^2 + volVel(i)^2 + volVel(i+1)^2)*3 - (volVel(i-1) + volVel(i) + volVel(i+1))^2);
        vel(i,1) = 1000 * (fL * sqrt(dxdtiltL^2 + dydtiltL^2))/scale; %Triple line velocity for the left side in um/uL
        
        %fR = (xVel(i,2) - xVel(i-1,2))/abs((xVel(i,2) - xVel(i-1,2))); %f for the right side(+1 for advancing movement and -1 for receding movement)
        dxdtiltR = ((xVel(i-1,2)*volVel(i-1) + xVel(i,2)*volVel(i) + xVel(i+1,2)*volVel(i+1))*3 - (volVel(i-1) + volVel(i) + volVel(i+1))*(xVel(i-1,2) + xVel(i,2) + xVel(i+1,2))) / ((volVel(i-1)^2 + volVel(i)^2 + volVel(i+1)^2)*3 - (volVel(i-1) + volVel(i) + volVel(i+1))^2);      
        dydtiltR = ((yVel(i-1,2)*volVel(i-1) + yVel(i,2)*volVel(i) + yVel(i+1,2)*volVel(i+1))*3 - (volVel(i-1) + volVel(i) + volVel(i+1))*(yVel(i-1,2) + yVel(i,2) + yVel(i+1,2))) / ((volVel(i-1)^2 + volVel(i)^2 + volVel(i+1)^2)*3 - (volVel(i-1) + volVel(i) + volVel(i+1))^2);
        vel(i,2) = 1000 * (fR * sqrt(dxdtiltR^2 + dydtiltR^2))/scale; %Triple line velocity for the right side in um/°
    end
end

fprintf('Fitting completed! \n');
%%
%---------------------------------------------------------------------------
%        7.       Saving and plotting results
%---------------------------------------------------------------------------
%
fprintf('------------------------ 7. SAVING AND PLOTTING RESULTS --------------------------- \n');
close all
fprintf('Pease wait, generating and storing results... \n');
% Predefine variables prior to loop
%{
volAvg = zeros(nStops,1);
volStd = zeros(nStops,1);
volEstWasherAvg = zeros(nStops,1);
volEstWasherStd = zeros(nStops,1);
x0Avg = zeros(nStops,2);
x0Std = zeros(nStops,2);
dispCPAvg = zeros(nStops,2);
dispCPStd = zeros(nStops,2);
dispCPmmAvg = zeros(nStops,2);
dispCPmmStd = zeros(nStops,2);
CircleCASAvg = zeros(nStops,2);
CircleCASStd = zeros(nStops,2);
PolyDropenCASAvg = zeros(nStops,2);
PolyDropenCASStd = zeros(nStops,2);
PolyCASAvg = zeros(nStops,2);
PolyCASStd = zeros(nStops,2);
EllipseCASAvg = zeros(nStops,2);
EllipseCASStd = zeros(nStops,2);
MaskCASAvg = zeros(nStops,2);
MaskCASStd = zeros(nStops,2);
%}
volAvg = NaN(nStops,1);
volStd = NaN(nStops,1);
volEstWasherAvg = NaN(nStops,1);
volEstWasherStd = NaN(nStops,1);
BaselineTiltAvg = NaN(nStops,1);
BaselineTiltStd = NaN(nStops,1);
x0Avg = NaN(nStops,2);
x0Std = NaN(nStops,2);
dispCPAvg = NaN(nStops,2);
dispCPStd = NaN(nStops,2);
dispCPmmAvg = NaN(nStops,2);
dispCPmmStd = NaN(nStops,2);
xVelAvg = NaN(nStops,2);
xVelStd = NaN(nStops,2);
yVelAvg = NaN(nStops,2);
yVelStd = NaN(nStops,2);
velAvg = NaN(nStops,2);
velStd = NaN(nStops,2);
CircleCASAvg = NaN(nStops,2);
CircleCASStd = NaN(nStops,2);
PolyDropenCASAvg = NaN(nStops,2);
PolyDropenCASStd = NaN(nStops,2);
PolyCASAvg = NaN(nStops,2);
PolyCASStd = NaN(nStops,2);
EllipseCASAvg = NaN(nStops,2);
EllipseCASStd = NaN(nStops,2);
MaskCASAvg = NaN(nStops,2);
MaskCASStd = NaN(nStops,2);
tCaptVec = linspace(0,tCapt,nImagesStops); %Capture interval time in s

% Analyze results within each nStop
for i = 1:nStops
    %- Calculate mean and standard deviation within each nStop
    volAvg(i,1) = mean(volExp((nImagesStops*(i-1)+1):i*nImagesStops,1));
    volStd(i,1) = std(volExp((nImagesStops*(i-1)+1):i*nImagesStops,1));
    volEstWasherAvg(i,1) = mean(volEstWasher((nImagesStops*(i-1)+1):i*nImagesStops,1));
    volEstWasherStd(i,1) = std(volEstWasher((nImagesStops*(i-1)+1):i*nImagesStops,1));
    volEstSphAvg(i,1) = mean(volEstSph((nImagesStops*(i-1)+1):i*nImagesStops,1));
    volEstSphStd(i,1) = std(volEstSph((nImagesStops*(i-1)+1):i*nImagesStops,1));
    BaselineTiltAvg(i,1) = mean(BaselineTilt((nImagesStops*(i-1)+1):i*nImagesStops,1));
    BaselineTiltStd(i,1) = std(BaselineTilt((nImagesStops*(i-1)+1):i*nImagesStops,1));   
    x0Avg(i,1) = mean(BaseVec((nImagesStops*(i-1)+1):i*nImagesStops,1));
    x0Std(i,1) = std(BaseVec((nImagesStops*(i-1)+1):i*nImagesStops,1));
    x0Avg(i,2) = mean(BaseVec((nImagesStops*(i-1)+1):i*nImagesStops,3));
    x0Std(i,2) = std(BaseVec((nImagesStops*(i-1)+1):i*nImagesStops,3));
    dispCPAvg(i,1) = mean(dispCP((nImagesStops*(i-1)+1):i*nImagesStops,1));
    dispCPStd(i,1) = std(dispCP((nImagesStops*(i-1)+1):i*nImagesStops,1));
    dispCPAvg(i,2) = mean(dispCP((nImagesStops*(i-1)+1):i*nImagesStops,2));
    dispCPStd(i,2) = std(dispCP((nImagesStops*(i-1)+1):i*nImagesStops,2));
    dispCPmmAvg(i,1) = mean(dispCPmm((nImagesStops*(i-1)+1):i*nImagesStops,1));
    dispCPmmStd(i,1) = std(dispCPmm((nImagesStops*(i-1)+1):i*nImagesStops,1));
    dispCPmmAvg(i,2) = mean(dispCPmm((nImagesStops*(i-1)+1):i*nImagesStops,2));
    dispCPmmStd(i,2) = std(dispCPmm((nImagesStops*(i-1)+1):i*nImagesStops,2));
    xVelAvg(i,1) = mean(xVel((nImagesStops*(i-1)+1):i*nImagesStops,1));
    xVelStd(i,1) = std(xVel((nImagesStops*(i-1)+1):i*nImagesStops,1));
    xVelAvg(i,2) = mean(xVel((nImagesStops*(i-1)+1):i*nImagesStops,2));
    xVelStd(i,2) = std(xVel((nImagesStops*(i-1)+1):i*nImagesStops,2));
    yVelAvg(i,1) = mean(yVel((nImagesStops*(i-1)+1):i*nImagesStops,1));
    yVelStd(i,1) = std(yVel((nImagesStops*(i-1)+1):i*nImagesStops,1));
    yVelAvg(i,2) = mean(yVel((nImagesStops*(i-1)+1):i*nImagesStops,2));
    yVelStd(i,2) = std(yVel((nImagesStops*(i-1)+1):i*nImagesStops,2));
    velAvg(i,1) = mean(vel((nImagesStops*(i-1)+1):i*nImagesStops,1));
    velStd(i,1) = std(vel((nImagesStops*(i-1)+1):i*nImagesStops,1));
    velAvg(i,2) = mean(vel((nImagesStops*(i-1)+1):i*nImagesStops,2));
    velStd(i,2) = std(vel((nImagesStops*(i-1)+1):i*nImagesStops,2));
    CircleCASAvg(i,1) = mean(CircleCAS((nImagesStops*(i-1)+1):i*nImagesStops,1));
    CircleCASStd(i,1) = std(CircleCAS((nImagesStops*(i-1)+1):i*nImagesStops,1));
    CircleCASAvg(i,2) = mean(CircleCAS((nImagesStops*(i-1)+1):i*nImagesStops,2));
    CircleCASStd(i,2) = std(CircleCAS((nImagesStops*(i-1)+1):i*nImagesStops,2));
    PolyDropenCASAvg(i,1) = mean(PolyDropenCAS((nImagesStops*(i-1)+1):i*nImagesStops,1));
    PolyDropenCASStd(i,1) = std(PolyDropenCAS((nImagesStops*(i-1)+1):i*nImagesStops,1));
    PolyDropenCASAvg(i,2) = mean(PolyDropenCAS((nImagesStops*(i-1)+1):i*nImagesStops,2));
    PolyDropenCASStd(i,2) = std(PolyDropenCAS((nImagesStops*(i-1)+1):i*nImagesStops,2));
    PolyCASAvg(i,1) = mean(PolyCAS((nImagesStops*(i-1)+1):i*nImagesStops,1));
    PolyCASStd(i,1) = std(PolyCAS((nImagesStops*(i-1)+1):i*nImagesStops,1));
    PolyCASAvg(i,2) = mean(PolyCAS((nImagesStops*(i-1)+1):i*nImagesStops,2));
    PolyCASStd(i,2) = std(PolyCAS((nImagesStops*(i-1)+1):i*nImagesStops,2));
    EllipseCASAvg(i,1) = mean(EllipseCAS((nImagesStops*(i-1)+1):i*nImagesStops,1));
    EllipseCASStd(i,1) = std(EllipseCAS((nImagesStops*(i-1)+1):i*nImagesStops,1));
    EllipseCASAvg(i,2) = mean(EllipseCAS((nImagesStops*(i-1)+1):i*nImagesStops,2));
    EllipseCASStd(i,2) = std(EllipseCAS((nImagesStops*(i-1)+1):i*nImagesStops,2));
    MaskCASAvg(i,1) = mean(MaskCAS((nImagesStops*(i-1)+1):i*nImagesStops,1));
    MaskCASStd(i,1) = std(MaskCAS((nImagesStops*(i-1)+1):i*nImagesStops,1));
    MaskCASAvg(i,2) = mean(MaskCAS((nImagesStops*(i-1)+1):i*nImagesStops,2));
    MaskCASStd(i,2) = std(MaskCAS((nImagesStops*(i-1)+1):i*nImagesStops,2));

    %- Ploting results for each nStops (x Position [pixels] && x displacement [mm] && vel [um/°] && CA [°] x time [s])
    plotName = strcat('Plot_xPosition&&xdispCPmm&&vel&&CAxTime_V(uL)_',num2str(round(volExp((nImagesStops*(i-1)+1),1),2),'%.2f'));
    figure('Name',plotName,'Units','normalized','Position',[0.25,0.05,0.5,0.85])
    tile = tiledlayout(3,1); %Create a tiled chart layout
    tile.TileSpacing = 'compact';
    tile.Padding = 'compact';
    %-- Top plot (x Position [pixels] x time [s])
    ax1 = nexttile;
    plot(ax1,tCaptVec,BaseVec((nImagesStops*(i-1)+1):i*nImagesStops,1)','--*','color','blue','DisplayName','Left CP');
    hold on
    plot(ax1,tCaptVec,BaseVec((nImagesStops*(i-1)+1):i*nImagesStops,3)','--x','color','red','DisplayName','Right CP');
    xlabel(ax1,'Time [s]')
    ylabel(ax1,'x Position [pixels]')
    %ylim(ax1,[0,1000]) %Limits of the x Axis
    %yticks(linspace(0,1000,5))
    lgdTop = legend;
    lgdTop.NumColumns = 1;
    lgdTop.Location = 'northeastoutside';
    %-- Medium plot ((x Displacement [mm] && velocity [um/°]) x time [s])
    ax2 = nexttile;
    plot(ax2,tCaptVec,dispCPStopmm((nImagesStops*(i-1)+1):i*nImagesStops,1)','--*','color','blue','DisplayName','xDisp left CP');
    hold on
    plot(ax2,tCaptVec,dispCPStopmm((nImagesStops*(i-1)+1):i*nImagesStops,2)','--x','color','red','DisplayName','xDisp right CP');
    xlabel(ax2,'Time [s]')
    ylabel(ax2,'x Displacement [mm]')
    yyaxis right
    plot(ax2,volExp,vel(:,1)',':*','color','cyan','DisplayName','vel left CP');
    hold on
    plot(ax2,volExp,vel(:,2)',':x','color','magenta','DisplayName','vel right CP');
    ylabel(ax2,strcat('velocity [',char(181),'m/',char(181),'L]'),'Color','k')
    ax2.YColor = 'k';
    %ylim(ax2,[0,6]) %Limits of the x Axis
    %yticks(linspace(0,1000,5))
    lgdTop = legend;
    lgdTop.NumColumns = 1;
    lgdTop.Location = 'northeastoutside';
    %-- Bottom plot (CA [°] x time [s])
    ax3 = nexttile;
    if UseCircFit
        plot(ax3,tCaptVec,CircleCAS((nImagesStops*(i-1)+1):i*nImagesStops,1)','--*','color','blue','DisplayName','CircleFit CAL');
        hold on
        plot(ax3,tCaptVec,CircleCAS((nImagesStops*(i-1)+1):i*nImagesStops,2)','--*','color','red','DisplayName','CircleFit CAR');
    end
    if UsePolyFitDropen
        plot(ax3,tCaptVec,PolyDropenCAS((nImagesStops*(i-1)+1):i*nImagesStops,1)','--+','color','blue','DisplayName','PolyFitDropen CAL');
        hold on
        plot(ax3,tCaptVec,PolyDropenCAS((nImagesStops*(i-1)+1):i*nImagesStops,2)','--+','color','red','DisplayName','PolyFitDropen CAR');
    end
    if UsePolyFit
        plot(ax3,tCaptVec,PolyCAS((nImagesStops*(i-1)+1):i*nImagesStops,1)','--square','color','blue','DisplayName','PolyFit CAL');
        hold on
        plot(ax3,tCaptVec,PolyCAS((nImagesStops*(i-1)+1):i*nImagesStops,2)','--square','color','red','DisplayName','PolyFit CAR');
    end
    if UseEllipseFit
        plot(ax3,tCaptVec,EllipseCAS((nImagesStops*(i-1)+1):i*nImagesStops,1)','--o','color','blue','DisplayName','EllipseFit CAL');
        hold on
        plot(ax3,tCaptVec,EllipseCAS((nImagesStops*(i-1)+1):i*nImagesStops,2)','--o','color','red','DisplayName','EllipseFit CAR');
    end
    if UseMaskMethod
        plot(ax3,tCaptVec,MaskCAS((nImagesStops*(i-1)+1):i*nImagesStops,1)','--^','color','blue','DisplayName','MaskMet CAL');
        hold on
        plot(ax3,tCaptVec,MaskCAS((nImagesStops*(i-1)+1):i*nImagesStops,2)','--^','color','red','DisplayName','MaskMet CAR');
    end
    xlabel(ax3,'Time [s]')
    ylabel(ax3,'CA [°]')
    %ylim(ax3,[70,85]) %Limits of the CA axis
    lgdBot = legend;
    lgdBot.NumColumns = 1;
    lgdBot.Location = 'northeastoutside';
    %-- Saving figure
    filePlotName = strcat(plotName,'.tif');
    fullfilePlotName = fullfile(outputPath,filePlotName);
    saveas(gcf,fullfilePlotName)

    %- Creating and saving table for each nStops
    fileOutputStopName = strcat(outputName,'_V(uL)_',num2str(round(volExp((nImagesStops*(i-1)+1),1),2),'%.2f'),outputFmt);
    fullfileOutputStopName = fullfile(outputPath,fileOutputStopName);
    varnames = {'Name','V[uL]','VestDisk[uL]','VestSph[uL]','BaselineTilt[°]','x0L[pxs]','x0R[pxs]','xDispCPL[pxs]','xDispCPR[pxs]','xDispCPL[mm]','xDispCPR[mm]','velL[um/uL]','velR[um/uL]','circleCAleft[°]','circleCAright[°]','polyDropenCAleft[°]','polyDropenCAright[°]','polyCAleft[°]','polyCAright[°]','ellipseCAleft[°]','ellipseCAright[°]','maskCAleft[°]','maskCAright[°]'}; %Nome das variáveis no arquivo txt
    T = table(fileNames((nImagesStops*(i-1)+1):i*nImagesStops)',round(volExp((nImagesStops*(i-1)+1):i*nImagesStops,1),2),round(volEstWasher((nImagesStops*(i-1)+1):i*nImagesStops,1),2),round(volEstSph((nImagesStops*(i-1)+1):i*nImagesStops,1),2),round(BaselineTilt((nImagesStops*(i-1)+1):i*nImagesStops,1),2),round(BaseVec((nImagesStops*(i-1)+1):i*nImagesStops,1),4),round(BaseVec((nImagesStops*(i-1)+1):i*nImagesStops,3),4),round(dispCPStop((nImagesStops*(i-1)+1):i*nImagesStops,1),4),round(dispCPStop((nImagesStops*(i-1)+1):i*nImagesStops,2),4),round(dispCPStopmm((nImagesStops*(i-1)+1):i*nImagesStops,1),4),round(dispCPStopmm((nImagesStops*(i-1)+1):i*nImagesStops,2),4),round(vel((nImagesStops*(i-1)+1):i*nImagesStops,1),4),round(vel((nImagesStops*(i-1)+1):i*nImagesStops,2),4),round(CircleCAS((nImagesStops*(i-1)+1):i*nImagesStops,1),4),round(CircleCAS((nImagesStops*(i-1)+1):i*nImagesStops,2),4),round(PolyDropenCAS((nImagesStops*(i-1)+1):i*nImagesStops,1),4),round(PolyDropenCAS((nImagesStops*(i-1)+1):i*nImagesStops,2),4),round(PolyCAS((nImagesStops*(i-1)+1):i*nImagesStops,1),4),round(PolyCAS((nImagesStops*(i-1)+1):i*nImagesStops,2),4),round(EllipseCAS((nImagesStops*(i-1)+1):i*nImagesStops,1),4),round(EllipseCAS((nImagesStops*(i-1)+1):i*nImagesStops,2),4),round(MaskCAS((nImagesStops*(i-1)+1):i*nImagesStops,1),4),round(MaskCAS((nImagesStops*(i-1)+1):i*nImagesStops,2),4),'VariableNames',varnames);
    writetable(T,fullfileOutputStopName,'Delimiter',' ','WriteVariableNames',true) %Write table to .txt file
end
%- Ploting results considering all data (xPos && dispCP && CA x Volume)
plotName = strcat('Plot_xPos&&xdispCPmm&&vel&&CAxVolume');
figure('Name',plotName,'Units','normalized','Position',[0.25,0.05,0.5,0.85])
tile = tiledlayout(3,1); %Create a tiled chart layout
tile.TileSpacing = 'compact';
tile.Padding = 'compact';
%-- Top plot (x Position [pixels] x Volume [uL])
ax1 = nexttile;
errorbar(ax1,volAvg(:,1)',x0Avg(:,1)',x0Std(:,1)','--*','color','blue','MarkerSize',6,'DisplayName','Left CP');
hold on
errorbar(ax1,volAvg(:,1)',x0Avg(:,2)',x0Std(:,2)','--x','color','red','MarkerSize',6,'DisplayName','Right CP');
xlabel(ax1,strcat('Volume [',char(181),'L]'))
ylabel(ax1,'x Position [pixels]')
%ylim(ax1,[0,1000]) %Limits of the x Axis
%yticks(linspace(0,1000,5))
lgdTop = legend;
lgdTop.NumColumns = 1;
lgdTop.Location = 'northeastoutside';
%-- Medium plot ((x Displacement [mm] && vel [um/uL]) x Volume [uL])
ax2 = nexttile;
errorbar(ax2,volAvg(:,1)',dispCPmmAvg(:,1)',dispCPmmStd(:,1)','--*','color','blue','MarkerSize',6,'DisplayName','xDisp left CP');
hold on
errorbar(ax2,volAvg(:,1)',dispCPmmAvg(:,2)',dispCPmmStd(:,2)','--x','color','red','MarkerSize',6,'DisplayName','xDisp right CP');
ylabel(ax2,'Displacement [pixel]')
xlabel(ax2,strcat('Volume [',char(181),'L]'))
yyaxis right
errorbar(ax2,volAvg(:,1)',velAvg(:,1)',velStd(:,1)',':*','color','cyan','MarkerSize',6,'DisplayName','vel left CP');
hold on
errorbar(ax2,volAvg(:,1)',velAvg(:,2)',velStd(:,2)',':x','color','magenta','MarkerSize',6,'DisplayName','vel right CP');
ylabel(ax2,strcat('velocity [',char(181),'m/',char(181),'L]'),'Color','k')
ax2.YColor = 'k';
%ylim(ax1,[0,1000]) %Limits of the x Axis
%yticks(linspace(0,1000,5))
lgdMed = legend;
lgdMed.NumColumns = 1;
lgdMed.Location = 'northeastoutside';
%-- Bottom plot (CA [°] x Volume [uL])
ax3 = nexttile;
if UseCircFit
errorbar(ax3,volAvg(:,1)',CircleCASAvg(:,1)',CircleCASStd(:,1)','--*','color','blue','MarkerSize',6,'DisplayName','CircleFit CAL');
hold on
errorbar(ax3,volAvg(:,1)',CircleCASAvg(:,2)',CircleCASStd(:,2)','--*','color','red','MarkerSize',6,'DisplayName','CircleFit CAR');
end
if UsePolyFitDropen
errorbar(ax3,volAvg(:,1)',PolyDropenCASAvg(:,1)',PolyDropenCASStd(:,1)','--+','color','blue','MarkerSize',6,'DisplayName','PolyFitDropen CAL');
hold on
errorbar(ax3,volAvg(:,1)',PolyDropenCASAvg(:,2)',PolyDropenCASStd(:,2)','--+','color','red','MarkerSize',6,'DisplayName','PolyFitDropen CAR');
end
if UsePolyFit
errorbar(ax3,volAvg(:,1)',PolyCASAvg(:,1)',PolyCASStd(:,1)','--square','color','blue','MarkerSize',6,'DisplayName','PolyFit CAL');
hold on
errorbar(ax3,volAvg(:,1)',PolyCASAvg(:,2)',PolyCASStd(:,2)','--square','color','red','MarkerSize',6,'DisplayName','PolyFit CAR');
end
if UseEllipseFit
errorbar(ax3,volAvg(:,1)',EllipseCASAvg(:,1)',EllipseCASStd(:,1)','--o','color','blue','MarkerSize',5,'DisplayName','EllipFit CAL');
hold on
errorbar(ax3,volAvg(:,1)',EllipseCASAvg(:,2)',EllipseCASStd(:,2)','--o','color','red','MarkerSize',5,'DisplayName','EllipFit CAR');
end
if UseMaskMethod
errorbar(ax3,volAvg(:,1)',MaskCASAvg(:,1)',MaskCASStd(:,1)','--^','color','blue','MarkerSize',5,'DisplayName','MaskMet CAL');
hold on
errorbar(ax3,volAvg(:,1)',MaskCASAvg(:,2)',MaskCASStd(:,2)','--^','color','red','MarkerSize',5,'DisplayName','MaskMet CAR');
end
xlabel(ax3,strcat('Volume [',char(181),'L]'))
ylabel(ax3,'CA [°]')
%ylim(ax2,[70,85]) %Limits of the CA axis
lgdBot = legend;
lgdBot.NumColumns = 1;
lgdBot.Location = 'northeastoutside';
%-- Saving figure
filePlotName = strcat(plotName,'.tif');
fullfilePlotName = fullfile(outputPath,filePlotName);
saveas(gcf,fullfilePlotName)

%- Ploting analysis of drop volume (Experimental Volume x Estimated Volume)
plotName = strcat('Plot_ExpVolumexEstVolume');
figure('Name',plotName)
errorbar(volAvg(:,1)',volEstWasherAvg(:,1)',volEstWasherStd(:,1)','--*','color','blue','MarkerSize',6,'DisplayName','Disk method');
hold on
errorbar(volAvg(:,1)',volEstSphAvg(:,1)',volEstSphStd(:,1)','--*','color','green','MarkerSize',6,'DisplayName','Spherical cap');
if length(volAvg) ~= 1
    equalLine = linspace(floor(min([min(volAvg),min(volEstWasherAvg),min(volEstSphAvg)])),ceil(max([max(volAvg),max(volEstWasherAvg),max(volEstSphAvg)])),10);
    plot(equalLine,equalLine,'--r','DisplayName','Expected behavior')
end
xlabel(strcat('Experimental Volume [',char(181),'L]'))
ylabel(strcat('Estimated Volume [',char(181),'L]'))
if length(volAvg) ~= 1
    xlim([floor(min([min(volAvg),min(volEstWasherAvg),min(volEstSphAvg)])) ceil(max([max(volAvg),max(volEstWasherAvg),max(volEstSphAvg)]))])
    ylim([floor(min([min(volAvg),min(volEstWasherAvg),min(volEstSphAvg)])) ceil(max([max(volAvg),max(volEstWasherAvg),max(volEstSphAvg)]))])
end
lgd = legend;
lgd.NumColumns = 1;
lgd.Location = 'northwest';
%-- Saving figure
filePlotName = strcat(plotName,'.tif');
fullfilePlotName = fullfile(outputPath,filePlotName);
saveas(gcf,fullfilePlotName)

%- Creating and saving table with average and standard deviation values for each stop
fileOutputAvgName = strcat(outputName,'_AvgData',outputFmt);
fullfileOutputAvgName = fullfile(outputPath,fileOutputAvgName);
varnames = {'V_avg[uL]','V_std[uL]','VestDisk_avg[uL]','VestDisk_std[uL]','VestSph_avg[uL]','VestSph_std[uL]','BaselineTilt_avg[°]','BaselineTilt_std[°]','x0L_avg[pxs]','x0L_std[pxs]','x0R_avg[pxs]','x0R_std[pxs]','xDispCPL_avg[pxs]','xDispCPL_std[pxs]','xDispCPR_avg[pxs]','xDispCPR_std[pxs]','xDispCPL_avg[mm]','xDispCPL_std[mm]','xDispCPR_avg[mm]','xDispCPR_std[mm]','velL_avg[um/uL]','velL_std[um/uL]','velR_avg[um/uL]','velR_std[um/uL]','circleCAleft_avg[°]','circleCAleft_std[°]','circleCAright_avg[°]','circleCAright_std[°]','polyDropenCAleft_avg[°]','polyDropenCAleft_std[°]','polyDropenCAright_avg[°]','polyDropenCAright_std[°]','polyCAleft_avg[°]','polyCAleft_std[°]','polyCAright_avg[°]','polyCAright_std[°]','ellipseCAleft_avg[°]','ellipseCAleft_std[°]','ellipseCAright_avg[°]','ellipseCAright_std[°]','maskCAleft_avg[°]','maskCAleft_std[°]','maskCAright_avg[°]','maskCAright_std[°]'}; %Nome das variáveis no arquivo txt
T = table(round(volAvg(:,1),2),round(volStd(:,1),2),round(volEstWasherAvg(:,1),2),round(volEstWasherStd(:,1),2),round(volEstSphAvg(:,1),2),round(volEstSphStd(:,1),2),round(BaselineTiltAvg(:,1),2),round(BaselineTiltStd(:,1),2),x0Avg(:,1),round(x0Std(:,1),4),round(x0Avg(:,2),4),round(x0Std(:,2),4),round(dispCPAvg(:,1),4),round(dispCPStd(:,1),4),round(dispCPAvg(:,2),4),round(dispCPStd(:,2),4),round(dispCPmmAvg(:,1),4),round(dispCPmmStd(:,1),4),round(dispCPmmAvg(:,2),4),round(dispCPmmStd(:,2),4),round(velAvg(:,1),4),round(velStd(:,1),4),round(velAvg(:,2),4),round(velStd(:,2),4),round(CircleCASAvg(:,1),4),round(CircleCASStd(:,1),4),round(CircleCASAvg(:,2),4),round(CircleCASStd(:,2),4),round(PolyDropenCASAvg(:,1),4),round(PolyDropenCASStd(:,1),4),round(PolyDropenCASAvg(:,2),4),round(PolyDropenCASStd(:,2),4),round(PolyCASAvg(:,1),4),round(PolyCASStd(:,1),4),round(PolyCASAvg(:,2),4),round(PolyCASStd(:,2),4),round(EllipseCASAvg(:,1),4),round(EllipseCASStd(:,1),4),round(EllipseCASAvg(:,2),4),round(EllipseCASStd(:,2),4),round(MaskCASAvg(:,1),4),round(MaskCASStd(:,1),4),round(MaskCASAvg(:,2),4),round(MaskCASStd(:,2),4),'VariableNames',varnames);
writetable(T,fullfileOutputAvgName,'Delimiter',' ','WriteVariableNames',true) %Write table to .txt file

%- Creating and saving table with all the results
outputFileName = strcat(outputName,'_AllData',outputFmt);
fullfileOutputName = fullfile(outputPath,outputFileName);
varnames = {'Name','V[uL]','VestDisk[uL]','VestSph[uL]','BaselineTilt[°]','x0L[pxs]','y0L[pxs]','x0R[pxs]','y0R[pxs]','xDispCPL[pxs]','xDispCPR[pxs]','xDispCPL[mm]','xDispCPR[mm]','velL[um/uL]','velR[um/uL]','circleCAleft[°]','circleCAright[°]','polyDropenCAleft[°]','polyDropenCAright[°]','polyCAleft[°]','polyCAright[°]','ellipseCAleft[°]','ellipseCAright[°]','maskCAleft[°]','maskCAright[°]'}; %Nome das variáveis no arquivo txt
T = table(fileNames',round(volExp(:,1),2),round(volEstWasher(:,1),2),round(volEstSph(:,1),2),round(BaselineTilt(:,1),2),round(BaseVec(:,1),4),round(BaseVec(:,2),4),round(BaseVec(:,3),4),round(BaseVec(:,4),4),round(dispCP(:,1),4),round(dispCP(:,2),4),round(dispCPmm(:,1),4),round(dispCPmm(:,2),4),round(vel(:,1),4),round(vel(:,2),4),round(CircleCAS(:,1),4),round(CircleCAS(:,2),4),round(PolyDropenCAS(:,1),4),round(PolyDropenCAS(:,2),4),round(PolyCAS(:,1),4),round(PolyCAS(:,2),4),round(EllipseCAS(:,1),4),round(EllipseCAS(:,2),4),round(MaskCAS(:,1),4),round(MaskCAS(:,2),4),'VariableNames',varnames);
writetable(T,fullfileOutputName,'Delimiter',' ','WriteVariableNames',true) %Write table to .txt file

%- Creating and saving table with settings
fileOptName = strcat(outputName,'_Settings',outputFmt);
fullfileOptName = fullfile(outputPath,fileOptName);
varnames = {'Date','scale[pxs/mm]','EnhanceContrast','low_in','high_in','low_out','high_out','SmoothingFilter','sigmaSmoothFilter','sizeSmoothFilter','EdgeDetMethodOpt','edgeDetMetName','edgeDetDirection','edgeDetSigma','edgeDetThreshold','edgeDetFilter','edgeDetFilterTypeOpt','edgeDetFilterTypeName','edgeDetFilterSize','edgeDetFilterSigma','edgeDetFilterShape','edgeDetDifThreshold','edgeDetOrder','edgeDetSmoothingIter','margin','jump_dist','filterMode','MorphOpen','openP','defaultopenP','openConn','MorphBridge','bridgeN','MorphClose','closeN','MorphThin','thinN','BaselineDetOpt','yDistWeight','marginBaselineDet','nBaselineDet','cutBaselineDet','percentBaselineDet','fillBlankHorLines','marginMaskBaselineDet','maskSize','limBaselineTilt','plotIntEdgeDetImage','UseAverageBaseVec','UseCircFit','n_circfit','PlotCircleImages','UsePolyFit','n_polyfit','poly_ord','PlotPolynomialImages','UsePolyFitDropen','n_polyfitDropen','poly_ordDropen','PlotPolynomialDropenImages','UseEllipseFit','PlotEllipticImages','UseMaskMethod','PlotMaskImages','plotIntFittingImages'}; %Nome das variáveis no arquivo txt
if isempty(edgeDetThreshold)
    edgeDetThresholdValue = 0;
end
T = table(datetime,round(scale,4),EnhanceContrast,low_in,high_in,low_out,high_out,SmoothingFilter,sigmaSmoothFilter,sizeSmoothFilter,EdgeDetMethodOpt,string(edgeDetMetName),string(edgeDetDirection),edgeDetSigma,edgeDetThresholdValue,edgeDetFilter,string(edgeDetFilterTypeOpt),string(edgeDetFilterTypeName),edgeDetFilterSize,edgeDetFilterSigma,edgeDetFilterShape,edgeDetDifThreshold,edgeDetOrder,edgeDetSmoothingIter,margin,jump_dist,string(filterMode),MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,BaselineDetOpt,yDistWeight,marginBaselineDet,nBaselineDet,cutBaselineDet,percentBaselineDet,fillBlankHorLines,marginMaskBaselineDet,maskSize,limBaselineTilt,plotIntEdgeDetImage,UseAverageBaseVec,UseCircFit,n_circfit,PlotCircleImages,UsePolyFit,n_polyfit,poly_ord,PlotPolynomialImages,UsePolyFitDropen,n_polyfitDropen,poly_ordDropen,PlotPolynomialDropenImages,UseEllipseFit,PlotEllipticImages,UseMaskMethod,PlotMaskImages,plotIntFittingImages,'VariableNames',varnames);
writetable(rows2vars(T),fullfileOptName,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file

%- Update of saved variables
%routVarNames = {'EnhanceContrast','low_in','high_in','low_out','high_out','SmoothingFilter','sigmaSmoothFilter','sizeSmoothFilter','EdgeDetMethodOpt','edgeDetMetName','edgeDetDirection','edgeDetSigma','edgeDetThreshold','edgeDetFilter','edgeDetFilterTypeOpt','edgeDetFilterTypeName','edgeDetFilterSize','edgeDetFilterSigma','edgeDetFilterShape','edgeDetDifThreshold','edgeDetOrder','edgeDetSmoothingIter','margin','jump_dist','filterMode','MorphOpen','openP','defaultopenP','openConn','MorphBridge','bridgeN','MorphClose','closeN','MorphThin','thinN','BaselineDetOpt','yDistWeight','marginBaselineDet','yAverageLR','nBaselineDet','cutBaselineDet','percentBaselineDet','fillBlankHorLines','marginMaskBaselineDet','maskSize','limBaselineTilt','plotIntEdgeDetImage','UseAverageBaseVec','UseCircFit','n_circfit','PlotCircleImages','UsePolyFit','n_polyfit','poly_ord','PlotPolynomialImages','UsePolyFitDropen','n_polyfitDropen','poly_ordDropen','PlotPolynomialDropenImages','UseEllipseFit','PlotEllipticImages','UseMaskMethod','PlotMaskImages','plotIntFittingImages'};
save('VarIntSessDrop.mat',routVarNames{:});

fprintf('Evaluation completed! Files sucessfully created and saved! \n');
end

