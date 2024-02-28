function [baseL,baseR,dropCenter] = ManualSelection(ID)
%MANUALSELECTION Manual selection of baseline and vertical drop axis of 
% symmetry coordinates
%   INPUT:
% ID: image
%   OUTPUT:
% x_CPL and y_CPL: coordinates of the left contact point
% x_CPR and y_CPR: coordinates of the right contact point

resp = "Y";
while resp ~= "N"
    fprintf('- Determination of reference points... \n');
    manualSelImage = figure;
    set(manualSelImage,'Name','Manual baseline determination');
    movegui(manualSelImage,'north');
    %Ibase=ID;
    imshow(ID);
    hold on
    % identify first contact point
    fprintf("Please, select the LEFT CONTACT POINT... \n")
    fprintf("Left click for ZOOM, when zoom is okay press any key to proceed... \n")
    text(1,-10,'Selection of LEFT BASELINE POINT (ZOOM MODE)',...
        'FontSize',8,'color','w','BackgroundColor',[.7 .7 .7])
    zoom on; %Enable zoom option
    pause(); %Wait any key to proceed
    zoom off; %Escap the zoom mode
    %text(1,-10,'Left click for zoom with pan and right click to select the LEFT BASELINE POINT',...
        %'FontSize',8,'color','w','BackgroundColor',[.7 .7 .7])
    fprintf("Left click to select the LEFT BASELINE POINT... \n")
    [x_base1, y_base1] = ginput(1);
    %hold on
    plot(x_base1, y_base1,'yx','MarkerSize',15)
    fprintf("Point seleceted! \n")
    zoom out; % Go to the original size of your image

    % identify second contact point
    fprintf("Please, select the RIGHT CONTACT POINT... \n")
    fprintf("Left click for ZOOM, when zoom is okay press any key to proceed... \n")
    text(1,-10,'Selection of RIGHT BASELINE POINT (ZOOM MODE)',...
        'FontSize',8,'color','w','BackgroundColor',[.7 .7 .7])
    zoom on; %Enable zoom option
    pause(); %Wait any key to proceed
    zoom off; %Escap the zoom mode
    %text(1,-10,'Left click for zoom with pan and right click to select the RIGHT BASELINE POINT',...
        %'FontSize',8,'color','w','BackgroundColor',[.7 .7 .7])
    fprintf("Left click to select the RIGHT BASELINE POINT... \n")
    [x_base2, y_base2] = ginput(1);
    plot(x_base2, y_base2,'yx','MarkerSize',15)
    fprintf("Point seleceted! \n")
    zoom out; % Go to the original size of your image

    % identify drop axis of symmetry
    fprintf("Please, select a point to identify the DROP VERTICAL AXIS OF SYMMETRY... \n")
    fprintf("Left click to select the DROP VERTICAL AXIS OF SYMMETRY... \n")
    text(1,-10,'Left click to select the DROP VERTICAL AXIS OF SYMMETRY',...
        'FontSize',8,'color','w','BackgroundColor',[.7 .7 .7])
    [x_AxisSym, y_AxisSym] = ginput(1);
    plot(x_AxisSym, y_AxisSym,'rx','MarkerSize',15)
    fprintf("Point seleceted! \n")

    % Plotting lines
    res_v = size(ID,1);
    res_h = size(ID,2);
    plot([x_AxisSym;x_AxisSym],[1;res_v],'--r','LineWidth',2); %Plot needle vertical axis of symmetry
    y_base = (y_base1 + y_base2)/2;
    plot([1;res_h],[y_base,y_base],'--y','LineWidth',2); %Plot baseline

    %Calculation of the baseline tilt
    tiltAngle = atan2(y_base2-y_base1,x_base2-x_base1) * 180/pi;
    fprintf("Baseline tilt angle: %.2f° \n",tiltAngle);

    resp = input('Do you like to make another manual identification (Y/N) [N]? ','s');
    if isempty(resp)
        resp = "N";
    end
    close(manualSelImage)
end
baseL = [x_base1, y_base1];
baseR = [x_base2, y_base2];
dropCenter = [x_AxisSym,y_base];

%{
% baseline points
x_baseline = (1:size(ID,2));
y_baseline = round((y_base2 - y_base1) / (x_base2 - x_base1) .* (x_baseline - x_base1) + y_base1);
y_baseline(y_baseline>size(ID,1)) = size(ID,1);
plot(x_baseline, y_baseline,'r')

x_CPL=round(y_base1,0);
y_CPL=round(x_base1,0);
x_CPR=round(y_base2,0);
y_CPR=round(x_base2,0);
%}

end

