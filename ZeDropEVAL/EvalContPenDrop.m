function EvalContPenDrop()
%EVALCONTPENDROP Evaluation of continuous pendant drop quasi-static 
% experiment (results from Giada).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used for the analysis of continuous pendant drop quasi-static 
% experiment. Reads .gam results from Giada and organize them in tables.
% Analyzes .gam files resulting from the Giada software organizing them 
% into a single file for later export to Excel.
% Created by Ayrton Pereira at 02/18/23.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
fprintf('---- EVALUATION OF CONTINUOUS PENDANT DROP EXPERIMENT (RESULTS FROM GIADA) ----\n');
%Read files
fprintf('Select the directory where the .gam data is located... \n');
filePath = uigetdir('C:\'); %Select the directory where the data is located
folderinfo = dir(filePath); %Structure containing information about all the files in the filePath
fileNames = {folderinfo(3:end).name}; %Get only the name of the interested files
fullfileName = fullfile(filePath,fileNames);
outputName = input('Enter a name for the exported data files: ','s');
outputFmt = '.txt';
outputFileName = strcat(outputName,outputFmt);
fullfileOutputName = fullfile(filePath,outputFileName);


%Create variables
X0 = zeros(1,length(fullfileName));
Y0 = zeros(1,length(fullfileName));
R0 = zeros(1,length(fullfileName));
beta = zeros(1,length(fullfileName));
gamma = zeros(1,length(fullfileName));
Vgon = zeros(1,length(fullfileName));
Vgiada = zeros(1,length(fullfileName));
A = zeros(1,length(fullfileName));

%Reading .gam files
for k = 1 : length(fullfileName)
    fileID = fopen(fullfileName{k});
    C = textscan(fileID,'%s %s %d %f %f %f %f %f %f %f');
    name{k} = C{1}{1}; %Name of the profile
    X0(k) = C{4}; %X-coordinate of the apex in pxs
    Y0(k) = C{5}; %Y-coordinate of the apex in pxs
    R0(k) = C{6}; %Radius of curvature at the apex in pxs
    beta(k) = C{7}; %Bond number
    gamma(k) = C{8}; %Surface tension of the liquid in mJ/m
    Vgon(k) = str2double(extractAfter(C{1}{1},"V(uL)_")); %Drop volume obtained directly from the goniometer routine
    %Vgon(k) = str2double(extractBetween(C{1}{1},"V(uL)_","_")); %Drop volume obtained directly from the goniometer routine
    Vgiada(k) = C{10}; %drop volume in uL obtained from giada
    A(k) = C{9}; %Drop surface area in mm²
end

%Creating and saving table
varnames = {'Name','X0[pxs]','Y0[pxs]','R0[pxs]','beta','gamma[mN/m]','Vgon[uL]','Vgiada[uL]','A[mm²]'}; %Variable names in .txt file
T = table(name',X0',Y0',R0',beta',gamma',Vgon',Vgiada',A','VariableNames',varnames);
writetable(T,fullfileOutputName,'Delimiter',' ','WriteVariableNames',true) %Write table to .txt file

fprintf('Evaluation completed! Files sucessfully created and saved! \n');
end

