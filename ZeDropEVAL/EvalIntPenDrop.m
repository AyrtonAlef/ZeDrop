function EvalIntPenDrop()
%EVALINTPENDROP Evaluation of intermittent pendant drop quasi-static 
% experiment (results from Giada).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used for the analysis of intermittent pendant drop quasi-static 
% experiment. Reads .gam results from Giada and organize them in tables.
% Analyzes .gam files resulting from the Giada software organizing them 
% into a single file for later export to Excel.
% Created by Ayrton Pereira at 02/18/23.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
fprintf('---- EVALUATION OF INTERMITTENT PENDANT DROP EXPERIMENT (RESULTS FROM GIADA) ----\n');
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
nImages = input('Enter the number of images in each stop [20]: ');
if isempty(nImages)
    nImages = 20;
end

%Create variables
name = cell((length(fullfileName)./nImages),nImages);
X0 = zeros((length(fullfileName)./nImages),nImages);
Y0 = zeros((length(fullfileName)./nImages),nImages);
R0 = zeros((length(fullfileName)./nImages),nImages);
beta = zeros((length(fullfileName)./nImages),nImages);
gamma = zeros((length(fullfileName)./nImages),nImages);
Vgon = zeros((length(fullfileName)./nImages),nImages);
Vgiada = zeros((length(fullfileName)./nImages),nImages);
A = zeros((length(fullfileName)./nImages),nImages);

X0Avg = zeros((length(fullfileName)./nImages),1);
X0Stdev = zeros((length(fullfileName)./nImages),1);
Y0Avg = zeros((length(fullfileName)./nImages),1);
Y0Stdev = zeros((length(fullfileName)./nImages),1);
R0Avg = zeros((length(fullfileName)./nImages),1);
R0Stdev = zeros((length(fullfileName)./nImages),1);
betaAvg = zeros((length(fullfileName)./nImages),1);
betaStdev = zeros((length(fullfileName)./nImages),1);
gammaAvg = zeros((length(fullfileName)./nImages),1);
gammaStdev = zeros((length(fullfileName)./nImages),1);
VgonAvg = zeros((length(fullfileName)./nImages),1);
VgonStdev = zeros((length(fullfileName)./nImages),1);
VgiadaAvg = zeros((length(fullfileName)./nImages),1);
VgiadaStdev = zeros((length(fullfileName)./nImages),1);
AreaAvg = zeros((length(fullfileName)./nImages),1);
AreaStdev = zeros((length(fullfileName)./nImages),1);

%Reading .gam files
fprintf('Reading .gam files... \n');
i = 1;
j = 1;
for k = 1 : length(fullfileName)
    fileID = fopen(fullfileName{k}); %Open .gam file
    C = textscan(fileID,'%s %s %d %f %f %f %f %f %f %f'); %Scan data from . gam file
    name{i,j} = C{1}{1}; %Name of the profile
    X0(i,j) = C{4}; %X-coordinate of the apex in pxs
    Y0(i,j) = C{5}; %Y-coordinate of the apex in pxs
    R0(i,j) = C{6}; %Radius of curvature at the apex in pxs
    beta(i,j) = C{7}; %Bond number
    gamma(i,j) = C{8}; %Surface tension of the liquid in mJ/m
    Vgon(i,j) = str2double(extractBetween(C{1}{1},"V(uL)_","_")); %Drop volume obtained directly from the goniometer routine
    Vgiada(i,j) = C{10}; %drop volume in uL obtained from giada
    A(i,j) = C{9}; %Drop surface area in mm²
    j = j + 1;
    if mod(k,nImages) == 0 %Multiple of nImages
        % Calculating mean and standard deviation for each volume
        X0Avg(i,1) = mean(X0(i,:));
        X0Stdev(i,1) = std(X0(i,:));
        Y0Avg(i,1) = mean(Y0(i,:));
        Y0Stdev(i,1) = std(Y0(i,:));
        R0Avg(i,1) = mean(R0(i,:));
        R0Stdev(i,1) = std(R0(i,:));
        betaAvg(i,1) = mean(beta(i,:));
        betaStdev(i,1) = std(beta(i,:));
        gammaAvg(i,1) = mean(gamma(i,:));
        gammaStdev(i,1) = std(gamma(i,:));
        VgonAvg(i,1) = mean(Vgon(i,:));
        VgonStdev(i,1) = std(Vgon(i,:));
        VgiadaAvg(i,1) = mean(Vgiada(i,:));
        VgiadaStdev(i,1) = std(Vgiada(i,:));
        AreaAvg(i,1) = mean(A(i,:));
        AreaStdev(i,1) = std(A(i,:));
        % Increase i and reset j to read another volume
        i = i + 1;
        j = 1;
    end
    fclose(fileID); %Close file
end

%Creating and saving table for each series of volume measurement
fprintf('Creating and saving data for each stop... \n');
for k = 1:(length(fullfileName)/nImages)
    tableName = strcat(outputName,'_V(uL)_',num2str(Vgon(k,1),'%.2f'),outputFmt);
    fullfileTableName = fullfile(filePath,tableName);
    varnames = {'Name','X0[pxs]','Y0[pxs]','R0[pxs]','beta','gamma[mN/m]','Vgon[uL]','Vgiada[uL]','A[mm²]'}; %Variable names in .txt file
    T = table(name(k,:)',X0(k,:)',Y0(k,:)',R0(k,:)',beta(k,:)',gamma(k,:)',Vgon(k,:)',Vgiada(k,:)',A(k,:)','VariableNames',varnames);
    writetable(T,fullfileTableName,'Delimiter',' ','WriteVariableNames',true) %Write table to .txt file
end
%Creating and saving table for series measurement
fprintf('Creating and saving data for all data... \n');
varnames = {'X0_avg[pxs]','X0_stdev[pxs]','Y0_avg[pxs]','Y0_stdev[pxs]','R0_avg[pxs]','R0_stdev[pxs]','beta_avg','beta_stdev','gamma_avg[mN/m]','gamma_stdev[mN/m]','Vgon_avg[uL]','Vgon_stdev[uL]','Vgiada_avg[uL]','Vgiada_stdev[uL]','A_avg[mm²]','A_stdev[mm²]'}; %Variable names in .txt file
Tsum = table(X0Avg(:,1),X0Stdev(:,1),Y0Avg(:,1),Y0Stdev(:,1),R0Avg(:,1),R0Stdev(:,1),betaAvg(:,1),betaStdev(:,1),gammaAvg(:,1),gammaStdev(:,1),VgonAvg(:,1),VgonStdev(:,1),VgiadaAvg(:,1),VgiadaStdev(:,1),AreaAvg(:,1),AreaStdev(:,1),'VariableNames',varnames);
writetable(Tsum,fullfileOutputName,'Delimiter',' ','WriteVariableNames',true) %Write table to .txt file

fprintf('Evaluation completed! Files sucessfully created and saved! \n');
end

