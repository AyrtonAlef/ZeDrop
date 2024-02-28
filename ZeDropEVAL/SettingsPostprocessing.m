function [MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode] = SettingsPostprocessing(MorphOpen,openP,defaultopenP,openConn,MorphBridge,bridgeN,MorphClose,closeN,MorphThin,thinN,margin,jump_dist,filterMode)
%SETTINGSPOSPROCESSING Enable the user to edit pos-processing image settings
%   INPUT/OUTPUT:
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

%Post-processing settings
%- Current settings
fprintf('-------------------------------------------------------------- \n');
fprintf('- Current post-processing settings: \n');
fprintf(strcat("MorphOpen = ",num2str(MorphOpen)," \n"));
fprintf(strcat("openP = ",num2str(openP)," \n"));
fprintf(strcat("openConn = ",num2str(openConn)," \n"));
fprintf(strcat("MorphBridge = ",num2str(MorphBridge)," \n"));
fprintf(strcat("bridgeN = ",num2str(bridgeN)," \n"));
fprintf(strcat("MorphClose = ",num2str(MorphClose)," \n"));
fprintf(strcat("closeN = ",num2str(closeN)," \n"));
fprintf(strcat("MorphThin = ",num2str(MorphThin)," \n"));
fprintf(strcat("thinN = ",num2str(thinN)," \n"));
fprintf(strcat("margin = ",num2str(margin)," \n"));
fprintf(strcat("jump_dist = ",num2str(jump_dist)," \n"));
fprintf(strcat("filterMode = ",filterMode," \n"));
fprintf('-------------------------------------------------------------- \n');

%- Edit settings
PostProcessingOpt = 0;
while PostProcessingOpt ~= 1 && PostProcessingOpt ~= 2 && PostProcessingOpt ~= 3
    while 1
        fprintf('- Edit post-processing settings: \n');
        fprintf('1. Morphological operations \n');
        fprintf('2. Identification of the longest edge \n');
        fprintf('3. Identification and sort of left and right edges \n');
        PostProcessingOpt = input('Select an option: ');
        if ~isempty(PostProcessingOpt)
            break;
        end
    end
end
if PostProcessingOpt == 1 %Morphological operation settings
    fprintf('Morphological operations settings... \n')
    %- Open operation
    ansMorphOpen = input('Do you like to enable OPEN operation (Y/N) [N]? ','s');
    if isempty(ansMorphOpen)
        ansMorphOpen = 'N';
    end
    if ansMorphOpen == 'N'
        MorphOpen = 0;
    else
        MorphOpen = 1;
        openP = input('Maximum number of pixels in objects [100]: ');
        if isempty(openP)
            openP = 100;
        end
        if openP ~= 100
            defaultopenP = 0;
        end
        while 1
            openConn = input('Pixel connectivity (4/8) [8]: ');
            if isempty(openConn)
                openConn = 8;
            end
            if openConn ~= 4 && openConn ~= 8
                fprintf('Invalid value. Please, enter a new value. \n')
            else
                break;
            end
        end
    end
    %- Bridge operation
    ansMorphBridge = input('Do you like to enable BRIDGE operation (Y/N) [N]? ','s');
    if isempty(ansMorphBridge)
        ansMorphBridge = 'N';
    end
    if ansMorphBridge == 'N'
        MorphBridge = 0;
    else
        MorphBridge = 1;
        bridgeN = input('Number of times to perform bridge operation [Inf]: ');
        if isempty(bridgeN)
            bridgeN = Inf;
        end
    end
    %- Close operation
    ansMorphClose = input('Do you like to enable CLOSE operation (Y/N) [N]? ','s');
    if isempty(ansMorphClose)
        ansMorphClose = 'N';
    end
    if ansMorphClose == 'N'
        MorphClose = 0;
    else
        MorphClose = 1;
        closeN = input('Number of times to perform close operation [Inf]: ');
        if isempty(closeN)
            closeN = Inf;
        end
    end
    %- Thin operation
    ansMorphThin = input('Do you like to enable THIN operation (Y/N) [N]? ','s');
    if isempty(ansMorphThin)
        ansMorphThin = 'N';
    end
    if ansMorphThin == 'N'
        MorphThin = 0;
    else
        MorphThin = 1;
        thinN = input('Number of times to perform thin operation [Inf]: ');
        if isempty(thinN)
            thinN = Inf;
        end
    end
elseif PostProcessingOpt == 2 %Identification of the longest edge settings
    fprintf('Identification of longest edge settings... \n')
    margin = input('Margin in pixels [10]: ');
    if isempty(margin)
        margin = 10;
    end
elseif PostProcessingOpt == 3 %Sort of left and right edges settings
    fprintf('Identification and sort of left and right edges settings... \n')
    defaultJumpDist = margin; %Default value for jump distance;
    strJumpDist = strcat("Jump distance in pixels [",num2str(defaultJumpDist),"]: ");
    jump_dist = input(strJumpDist);
    if isempty(jump_dist)
        jump_dist = defaultJumpDist;
    end
    filterMode = input("Filter mode (H - Horizontal/ V - Vertical) [V]: ",'s');
    if isempty(filterMode) || filterMode ~= 'H'
        filterMode = 'V';
    end
end
end

