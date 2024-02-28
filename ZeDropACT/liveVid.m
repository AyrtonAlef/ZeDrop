function liveVid(vidobj,gridRow,gridCol)
%LIVEVID Preview of video
%   INPUT
% vidobj - video object
% gridrow - number of rows in the grid
% gridcol - number of columns in the grid

f = figure('Name', 'Live video','NumberTitle','off'); %Set live video figure
%uicontrol('String', 'Close', 'Callback', 'close(gcf)'); %Create close button
vidResValues = vidobj.VideoResolution; %Get resolution of the video object
vidSizeRatio = vidResValues(2)/vidResValues(1); %Calculate video object size ratio
nBands = vidobj.NumberOfBands;
fsize = 1000; %Width of the video that will appear on the figure
f.Position(3:4) = [fsize fsize*vidSizeRatio]; %Redimension live video figure
movegui('center') %Centralize live video figure in the screen
hImage = image( zeros(vidResValues(2), vidResValues(1), nBands) ); 
preview(vidobj, hImage); %Preview live video object
%%{
%-Create grid over live video
rows = vidResValues(2);
columns = vidResValues(1);
hold on;
%stepSizeY = vidResValues(2)/2; % Distance between lines in Y direction
%stepSizeX = vidResValues(1)/2; %Distance between lines in X direction
stepSizeY = vidResValues(2)/gridRow; % Distance between lines in Y direction
stepSizeX = vidResValues(1)/gridCol; %Distance between lines in X direction
for row = 1 : stepSizeY : rows
    yline(row, '-.r', 'LineWidth', 1);
end
for col = 1 : stepSizeX : columns
    xline(col, '-.r', 'LineWidth', 1);
end
%}
end

