function [x_AxisSym,y_Int,x_IntLeft,x_IntRight] = ManualSelectionNeedle(ID)
%MANUALSELECTIONNEEDLE Manual selection of needle position (vertical center
% of the needle) and horizontal interface between needle and drop
%   INPUT:
% ID: image
%   OUTPUT:
% x_AxisSym and y_AxisSym: coordinates of the vertical axis of symmetry of the needle
% x_IntLeft and y_IntLeft: coordinates of the left interface point between needle and drop
% x_IntRight and y_IntRight: coordinates of the right interface point between needle and drop

%fprintf("----------------- MANUAL IDENTIFICATION OF NEEDLE POSITION --------------- \n")
figure
resp = "Y";
while resp ~= "N"
    %Ibase=ID;
    imshow(ID);
    % identify needle vertical axis of symmetry
    fprintf("Please, select a point in the vertical axis of symmetry of the needle ... \n")
    text(1,-10,'Left click for zoom with pan and right click to select the VERTICAL AXIS OF SYMMETRY',...
        'FontSize',8,'color','w','BackgroundColor',[.7 .7 .7])
    [x_AxisSym, y_AxisSym] = ginput(1);
    hold on
    plot(x_AxisSym, y_AxisSym,'rx','MarkerSize',15)
    fprintf("Point seleceted! \n")

    % identify left interface point between needle and drop
    fprintf("Please, select the left interface point between needle and drop ... \n")
    text(1,-10,'Left click for zoom with pan and right click to select the LEFT INTERFACE POINT',...
        'FontSize',8,'color','w','BackgroundColor',[.7 .7 .7])
    [x_IntLeft, y_IntLeft] = ginput(1);
    plot(x_IntLeft, y_IntLeft,'yx','MarkerSize',15)
    fprintf("Point seleceted! \n")

    % identify right interface point between needle and drop
    fprintf("Please, select the right interface point between needle and drop ... \n")
    text(1,-10,'Left click for zoom with pan and right click to select the RIGHT INTERFACE POINT',...
        'FontSize',8,'color','w','BackgroundColor',[.7 .7 .7])
    [x_IntRight, y_IntRight] = ginput(1);
    plot(x_IntRight, y_IntRight,'yx','MarkerSize',15)
    fprintf("Point seleceted! \n")

    % Plotting lines
    res_v = size(ID,1);
    res_h = size(ID,2);
    plot([x_AxisSym;x_AxisSym],[1;res_v],'--r','LineWidth',2); %Plot needle vertical axis of symmetry
    plot([x_IntLeft;x_IntLeft],[1;res_v],'--y','LineWidth',2); %Plot vertical line perpassing the interface between needle and drop left profile
    plot([x_IntRight;x_IntRight],[1;res_v],'--y','LineWidth',2); %Plot vertical line perpassing the interface between needle and drop right profile
    y_Int = (y_IntLeft + y_IntRight)/2;
    plot([1;res_h],[y_Int,y_Int],'--y','LineWidth',2); %Plot horizontal line perpassing the interface between needle and drop
%{
% baseline points
x_baseline = (1:size(Ibase,2));
y_baseline = round((y_base2 - y_base1) / (x_base2 - x_base1) .* (x_baseline - x_base1) + y_base1);
y_baseline(y_baseline>size(Ibase,1)) = size(Ibase,1);
plot(x_baseline, y_baseline,'r')
x_center=round(y_base1,0);
y_center=round(x_base1,0);
x_left=round(y_base2,0);
y_left=round(x_base2,0);
x_right=round(y_base3,0);
y_right=round(x_base3,0);
%}
    resp = input('Do you like to make another manual identification (Y/N) [N]? ','s');
    if isempty(resp)
        resp = "N";
    end
    %close(gcf)
end
close(gcf)
end

