%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Summary of static drop results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear all
%Read files
fprintf('Select the directory where the static drop data (.txt) is located... \n');
filePath = uigetdir('C:\'); %Select the directory where the data is located
folderinfo = dir(filePath); %Structure containing information about all the files in the filePath
fileNames = {folderinfo(3:end).name}; %Get only the name of the interested files
fullfileName = fullfile(filePath,fileNames);
outputName = input('Enter a name for the exported data files: ','s');
outputFmt = '.txt';
outputFileName = strcat(outputName,outputFmt);
fullfileOutputName = fullfile(filePath,outputFileName);
%nImages = input('Enter the number of images in each stop [20]: ');
%if isempty(nImages)
    %nImages = 20;
%end
nFiles = size(fullfileName,2);
%Create variables
name = cell(length(fullfileName),1);
VestWasher = zeros(length(fullfileName),1);
VestSph = zeros(length(fullfileName),1);
xOL_pxs = zeros(length(fullfileName),1);
x0R_pxs = zeros(length(fullfileName),1);
circleCAleft = zeros(length(fullfileName),1);
circleCAright = zeros(length(fullfileName),1);
polyDropenCAleft = zeros(length(fullfileName),1);
polyDropenCAright = zeros(length(fullfileName),1);
polyCAleft = zeros(length(fullfileName),1);
polyCAright = zeros(length(fullfileName),1);
ellipseCAleft = zeros(length(fullfileName),1);
ellipseCAright = zeros(length(fullfileName),1);
maskCAleft = zeros(length(fullfileName),1);
maskCAright = zeros(length(fullfileName),1);

%Reading files
fprintf('Reading files... \n');
for k = 1 : length(fullfileName)
    fileID = fopen(fullfileName{k}); %Open file
    C = textscan(fileID,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',' '); %Scan data from . gam file
    name{k,1} = C{1}{2}; %Name of the profile
    VestWasher(k,1) = str2num(C{3}{2});
    VestSph(k,1) = str2num(C{4}{2});
    xOL_pxs(k,1) = str2num(C{6}{2});
    x0R_pxs(k,1) = str2num(C{7}{2});
    circleCAleft(k,1) = str2num(C{14}{2});
    circleCAright(k,1) = str2num(C{15}{2});
    polyDropenCAleft(k,1) = str2num(C{16}{2});
    polyDropenCAright(k,1) = str2num(C{17}{2});
    polyCAleft(k,1) = str2num(C{18}{2});
    polyCAright(k,1) = str2num(C{19}{2});
    ellipseCAleft(k,1) = str2num(C{20}{2});
    ellipseCAright(k,1) = str2num(C{21}{2});
    maskCAleft(k,1) = str2num(C{22}{2});
    maskCAright(k,1) = str2num(C{23}{2});
    fclose(fileID); %Close file
end

%Creating and saving table for each series of volume measurement
fprintf('Creating and saving data... \n');
%for k = 1:length(fullfileName)
    tableName = strcat(outputName,outputFmt);
    fullfileTableName = fullfile(filePath,tableName);
    varnames = {'Name','VestWasher[uL]','VestSph[uL]','x0L[pxs]','x0R[pxs]','circleCAleft[°]','circleCAright[°]','polyDropenCAleft[°]','polyDropenCAright[°]','polyCAleft[°]','polyCAright[°]','ellipseCAleft[°]','ellipseCAright[°]','maskCAleft[°]','maskCAright[°]'}; %Variable names in .txt file
    T = table(name(:,1),VestWasher(:,1),VestSph(:,1),xOL_pxs(:,1),x0R_pxs(:,1),circleCAleft(:,1),circleCAright(:,1),polyDropenCAleft(:,1),polyDropenCAright(:,1),polyCAleft(:,1),polyCAright(:,1),ellipseCAleft(:,1),ellipseCAright(:,1),maskCAleft(:,1),maskCAright(:,1),'VariableNames',varnames);
    writetable(T,fullfileTableName,'Delimiter',' ','WriteVariableNames',true) %Write table to .txt file
%end
fprintf('Evaluation completed! Files sucessfully created and saved! \n');