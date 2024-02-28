function [scale] = setScaleFactor(ID)
%SCALEFACTOR Determination of the scale factor [pixels/mm].
%   INPUT:
% ID: image
%   OUTPUT:
% scale: scale factor in pixels/mm

figure
resp = "Y";
while resp ~= "N"
    Ibase=ID;
    imshow(Ibase);
    % Select left point
    fprintf("Please, select the left point ... \n")
    text(1,-10,'Left click for zoom with pan and right click to select the LEFT POINT',...
        'FontSize',8,'color','w','BackgroundColor',[.7 .7 .7])
    [xleft, yleft] = ginput(1);
    hold on
    plot(xleft, yleft,'rx','MarkerSize',15)
    fprintf("Left point seleceted! \n")

    % Select right point
    fprintf("Please, select the right point ... \n")
    text(1,-10,'Left click for zoom with pan and right click to select the RIGHT POINT',...
        'FontSize',8,'color','w','BackgroundColor',[.7 .7 .7])
    [xright, yright] = ginput(1);
    plot(xright, yright,'yx','MarkerSize',15)
    fprintf("Right point seleceted! \n")

    % Confirmation of the selection
    x = [xleft xright];
    ymed = (yleft + yright)/2;
    y = [ymed ymed];
    plot(x, y,'--y','LineWidth',2); %Plot horizontal line between selected points
    resp = input('Do you like to remake selections (Y/N) [N]? ','s');
    if isempty(resp)
        resp = "N";
    end
end
close(gcf)
% Calculation of the scale factor
pixels = xright - xleft; %calculation of dimensions in pixels 
mm = input('Enter size of scale in mm: '); %function that allows user to enter dimension in mm
scale = abs(pixels/mm);
fprintf('Scale factor determined: %.2f pxs/mm. \n',scale);
end

