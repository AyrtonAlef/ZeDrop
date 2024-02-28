function [incMesaTotal] = routineIntIncDrop(vidAdap,vidID,vidFmt,imFmt,smInc,passoMotor,nEntradasSemFim,nDentesCoroa,wMin,wMax,framerate,incMesaTotal)
%ROUTINEINTINCDROP Execution of intermittent inclined drop experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execution of intermittent inclined drop experiment. Performs sucessive
%   and intermittent tilts on the stage until reaching a pre-established
%   maximum tilt. Between tilts, drop images are captured allowing further
%   analysis. The routine requires the user to define the initial stage 
%   tilt (incIni), the maximum stage tilt (incMax), the stage tilt rate 
%   (taxaInc), the waiting time for drop equilibrium (tEspera), the 
%   time of image capture during a stop (tCapt), the number of captured 
%   images during a stop (nImages) and the number of stops along the 
%   intermittent tilt increase of the drop (nStops). The routine allows 
%   choosing the directory where the images will be saved, as well as 
%   defining the name of the images. At the end, a text file is exported 
%   with all relevant parameters of the experiment. *Requires the 
%   connection of the tilting stage system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% vidAdap - Character vector that specifies the name of the adaptor used to communicate with the device
% vidID - Numeric scalar value that identifies a particle device available through the specified adaptor 
% vidFmt - Character vector that specifies a particular video format supported by the device
% imFmt - Format of the captured image 
% smInc - Syringe motor control
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% nEntradasSemFim - Worm number of starts
% nDentesCoroa - Number of worm gear teeth
% wMin - Minimum stable rotation speed (RPM) in which the motor can reach. Corresponds to a sm.RPM = 2 (Microstep)
% wMax - Maximum stable rotation speed (RPM) in which the motor can reach
% framerate - Actual camera frame rate
% incMesaTotal - Current stage tilt [°]
%   OUTPUT
% incMesaTotal - Current stage tilt [°]

fprintf('-------------------- INCLINED DROP - INTERMITTENT EXPERIMENT -------------------\n');
% Experiment parameters (entered by the user)
fprintf('Experiment parameters: \n');
taxaIncMax = (wMax*360*nEntradasSemFim)/(nDentesCoroa*60); %Maximum stage tilt rate in °/s
taxaIncMin = (wMin*360*nEntradasSemFim)/(nDentesCoroa*60); %Minimum stage tilt rate in °/s
incStart = incMesaTotal; %Current stage tilt [°] 
incIni = input("- Initial tilt (°): "); %Initial stage tilt before starting the measurement series
prompttaxaIncIni = strcat("- Initial stage tilt rate (°/s) (max. ",num2str(taxaIncMax,'%.2f')," °/s): ");
while 1
    taxaIncIni = input(prompttaxaIncIni); %Initial stage tilt rate in °/s
    if taxaIncIni > taxaIncMax %Check that the stage tilt does not exceed the maximum motor rotation speed
        fprintf("Unable to proceed with the routine. Please decrease the initital tilt rate. \n")
    else
        break;
    end
end
incMax = input("- Maximum tilt (°): "); %Maximum tilt reached by the stage
prompttaxaInc = strcat("- Stage tilt rate (°/s) (max. ",num2str(taxaIncMax,'%.2f')," °/s): ");
while 1
    %taxaInc = input("- Stage tilt rate (°/s) (max. %.2f °/s): "); %Stage tilt rate in °/s
    taxaIncEns = input(prompttaxaInc); %Stage tilt rate in °/s
    if taxaIncEns > taxaIncMax %Check that the stage tilt rate does not exceed the maximum motor rotation speed
        fprintf("Unable to continue the routine. Please lower the stage tilt rate. \n")
    else
        break;
    end
end
sentidoInc = input("- Tilt direction (clockwise (CW)/ counterclockwise (CCW)) [CW]: ", "s");
if isempty(sentidoInc) || (sentidoInc ~= "CW") || (sentidoInc ~= "CCW")
    sentidoInc = "CW";
end
if sentidoInc == "CW"
    sentido = -1;
    sentidoRet = "CCW";
elseif sentidoInc == "CCW"
    sentido = 1;
    sentidoRet = "CW";
end
tEspera = input("- Waiting time for drop equilibrium (s): "); %Waiting time for drop equilibrium
tCapt = input("- Capture time during a stop (s): "); %Capture time during a stop
nImages = input("- Number of captured images during a stop: "); %Number of captured images during a stop
nStops = input("- Number of stops: "); %Number of stops

% Experiment parameters (calculated)
incVarIni = abs(incIni - incMesaTotal); %Initial tilt variation
incVarEns = abs(incMax - incIni)/nStops; %Tilt variation in each stop
incVarRet = abs(incMax - incStart); %Tilt variation to return the stage to its initial tilt position
tCaptInterval = tCapt/nImages; %Time interval between captured images in s
%volCaptInterval = taxaDepGota*tCaptInterval; %Variação de volume entre imagens em uL
totalFrame = framerate*tCapt; %Total number of frames acquired during the experiment
frameInterval = totalFrame/nImages; %Interval between frames to acquire the desired number of images
erroIncTotal = 0; %Error in total stage tilt [°]
clockEnsaio = 0; %Total time to execute the routine

%- Drive motor settings
nPassosIncIni = inc2passo(incVarIni,passoMotor,nEntradasSemFim,nDentesCoroa); %No of steps corresponding to the initial stage tilt 
nPassosIncEns = inc2passo(incVarEns,passoMotor,nEntradasSemFim,nDentesCoroa); %No of steps corresponding to the stage tilt in each stop

%- Image capture settings
fprintf('Path selection to export video and images... \n');
filePath = uigetdir('C:\'); %Path selection to export video and images
name = input('Enter a name for the video and image series: ','s');
vidDiskFmt = '.avi'; %Video format ('.avi','.mp4','.mj2')
%imFmt = '.tif'; %Images format ('.jpeg','.png','.tif')
vidComp = 'Motion JPEG AVI'; %Video compression ('Archival','Motion JPEG AVI','Motion JPEG 2000','MPEG-4','Uncompressed AVI','Indexed AVI' e 'Grayscale AVI')
vidQual = 90; %Video quality (related with the compression ratio)

%- Begin of experiment
fprintf('--------------------- EXECUTION OF TILTING ROUTINE ----------------------\n');
fprintf('Press any key to proceed with the tilting routine... \n');
pause;
tic
close all %Close all open windows
%delete(vidobj) %Delete any video objects that may have been created
clear vidobj vidsrc %Clear video related objects
clockEnsaio1 = toc;
clockEnsaio = clockEnsaio + clockEnsaio1; %Update experiment total time

%- Start of measurement series
%for k = 1:nStops
for k = 0:nStops
    tic
    %- Create video object
    [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

   %- Update of image names
    fileName = strcat(name,'_',num2str(k),'_alpha(°)_',num2str((incIni+k*incVarEns),'%.2f'));
    fullvidName = strcat(fileName,vidDiskFmt);
    fullfilevidPath = fullfile(filePath,fullvidName);

    %- Preparation of video object
    diskLogger = VideoWriter(fullfilevidPath,vidComp);
    diskLogger.FrameRate = framerate;%Matching the video frame rate to the actual camera frame rate
    diskLogger.Quality = vidQual; %Definition of the video quality
    vidobj.DiskLogger = diskLogger;
    %vidobj.FramesPerTrigger = framerate*tCapt; %Total frames to acquire
    vidobj.FramesPerTrigger = nImages; %Total frames to acquire
    vidobj.FrameGrabInterval = round(frameInterval); %Frame interval
    clockPrepVid = toc;
    if k == 0
        clockPrepIni = 0; %Time for the initial drop preparation
        clockPrepIni = clockPrepIni + clockPrepVid;
    else
        clockStop(k) = 0; %Time for the execution of each stop
        clockStop(k) = clockStop(k) + clockPrepVid; %Update time for stop execution
    end
    clockEnsaio = clockEnsaio + clockPrepVid; %Update experiment total time
    
    %- Stage tilt
    if k == 0
        %- Initital stage tilt
        tic
        fprintf('Please wait, initial stage tilt... \n');
        taxaInc = taxaIncIni;
        nPassosInc = nPassosIncIni;
    else
        tic
        %- Stop/measurement execution 
        if k >= 1
            fprintf("----------------------------------------------- \n");
            fprintf("Execution of stop/measurement %d of %d \n",k,nStops);
        end
        %- Additional stage tilt
        fprintf('Please wait, additional stage tilt... \n');
        taxaInc = taxaIncEns;
        nPassosInc = nPassosIncEns;
    end
    if taxaInc <= taxaIncMin %Stage tilt rate less than minimum motor rotation speed. Motor will work in pulses
        smInc.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a tilt rate = 0,0818 uL/s
        wMotorEns = incs2RPM(taxaInc,nEntradasSemFim,nDentesCoroa); %Motor rotation speed (RPM) in order to achieved the desired tilt rate (°/s)
        tadpasso = tadicionalpasso(passoMotor,wMotorEns,wMin);  %Time added to each step (in seconds) to reach the desired motor rotational speed
        for i = 1:floor(nPassosInc) %Execution of an additional tilt on the stage
            move(smInc,sentido*1); %Executes 1 rotation step on the motor
            pause(tadpasso);
        end
    elseif (taxaIncMin < taxaInc) && (taxaInc <= taxaIncMax)
        wRPM = fincs2RPM(taxaInc,nEntradasSemFim,nDentesCoroa); %Calculation of the rotational speed [RPM] necessary to be sent to the motor to reach the desired tilt rate [°/s]
        smInc.RPM = floor(wRPM);
        move(smInc,sentido*floor(nPassosInc)); %Motor movement
    end
    release(smInc); %Motor release
    if k == 0
        fprintf('Initial stage tilt completed! \n');
        erroIncIni = passo2inc(nPassosInc - floor(nPassosInc),passoMotor,nEntradasSemFim,nDentesCoroa); %Error in tilt during initital stage tilt
        erroIncTotal = erroIncTotal + erroIncIni;
        incMesaTotal = incMesaTotal - sentido*incVarIni; %Update total stage tilt
        clockIncIni = toc;
        clockPrepIni = clockPrepIni + clockIncIni; %Update time for the initial drop preparation
        clockEnsaio = clockEnsaio + clockPrepIni; %Update experiment total time
    else
        fprintf('Additional stage tilt completed! \n');
        erroIncEns = passo2inc(nPassosInc - floor(nPassosInc),passoMotor,nEntradasSemFim,nDentesCoroa); %Error in tilt during a stop
        erroIncTotal = erroIncTotal + erroIncEns;
        incMesaTotal = incMesaTotal - sentido*incVarEns; %Update total stage tilt
        clockIncEns = toc;
        clockStop(k) = clockStop(k) + clockIncEns; %Update time for the execution of a stop
        clockEnsaio = clockEnsaio + clockIncEns; %Update experiment total time
    end
    
    %- Waiting time for drop equilibrium
    tic
    fprintf('Please wait, waiting time for drop reach equilibrium... \n');
    pause(tEspera);
    fprintf('Waiting time completed! \n');
    clockEspGota = toc;
    if k == 0
        clockPrepInc = clockPrepInc + clockEspGota; %Update time for the initial drop preparation
    else
        clockStop(k) = clockStop(k) + clockEspGota; %Update time for the execution of a stop
        clockEnsaio = clockEnsaio + clockEspGota; %Update experiment total time
    end

    %- Start video recording
    tic
    fprintf('Start video recording. \n');
    start(vidobj) %Start video recording
    wait(vidobj) %Wait until all frames are captured
    stop(vidobj);
    % Wait for all frames to be written to disk
    while vidobj.FramesAcquired ~= vidobj.DiskLoggerFrameCount
        pause(.1);
    end
    fprintf('Video recording completed! \n');
    %toc
    if k == 0
        clockVidIni = toc;
        clockPrepIni = clockPrepIni + clockVidIni; %Update time for the initial drop preparation
        clockEnsaio = clockEnsaio + clockVidIni; %Update experiment total time
    else
        clockVid(k) = toc;
        clockStop(k) = clockStop(k) + clockVid(k); %Update time for the execution of a stop
        clockEnsaio = clockEnsaio + clockVid(k); %Update experiment total time
    end
        
    %- Showing an example image
    tic
    close all %Close all open windows
    frame = getsnapshot(vidobj); %Capture an image
    imshow(frame); %Show image
    movegui(gcf,'center') %Centralize image window
    set(gcf,'Name',strcat(fileName,'_Example'))
    clockEnsaio2 = toc;
    if k == 0
        clockPrepIni = clockPrepIni + clockEnsaio2; %Update time for the initial drop preparation
        clockEnsaio = clockEnsaio + clockPrepIni; %Update experiment total time
    else
        clockStop(k) = clockStop(k) + clockEnsaio2; %Update time for the execution of a stop
        clockEnsaio = clockEnsaio + clockStop(k); %Update experiment total time
    end

    %- Converting video to images
    tic
    fprintf('Please wait, converting video to images... \n');
    %vid2frame(filePath,fullvidName,fileName,imFmt,tCaptInterval);
    fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
    vidFile = VideoReader(fullfilevidPath); %Read Video File
    %i = 1;
    for m = 1: nImages %Get an image from every video frame
        frame = read(vidFile,m);
        framegray = rgb2gray(frame);
        fullimName = strcat(fileName,'_',num2str(m),imFmt);
        fullfileimPath = fullfile(filePath,fullimName);
        imwrite(framegray,fullfileimPath);
        %i = i + 1; %Interval of image name
    end
    delete(fullfilevidPath) %Delete video
    fprintf('Conversion of video to images completed! \n');
    if k == 0
        clockImagesIni = toc;
        clockPrepIni = clockPrepIni + clockImagesIni; %Update time for the initial drop preparation
        clockEnsaio = clockEnsaio + clockImagesIni; %Update experiment total time
    else
        clockImages(k) = toc;
        clockStop(k) = clockStop(k) + clockImages(k); %Update time for the execution of a stop
        clockEnsaio = clockEnsaio + clockImages(k); %Update experiment total time
    end

    %- Deleting video object
    tic
    delete(vidobj)
    clear vidobj vidsrc
    clockEnsaio3 = toc;
    if k == 0
        clockPrepIni = clockPrepIni + clockEnsaio3;
    else
        clockStop(k) = clockStop(k) + clockEnsaio3; %Update time for the execution of a stop
    end
    clockEnsaio = clockEnsaio + clockEnsaio3; %Update experiment total time 
end

%- Return to the starting stage tilt 
tic
fprintf('Please wait, return stage to the starting tilt ... \n');
moverIncMesa(smInc,passoMotor,nEntradasSemFim,nDentesCoroa,incVarRet,sentidoRet,"alta");
fprintf('Return of stage completed! \n');
incMesaTotal = incMesaTotal - incVarRet; %Update of total stage tilt
clockIncRet = toc;
clockEnsaio = clockEnsaio + clockIncRet; %Update experiment total time 

%- Export/save experiment description
filenametxt = "Description.txt";
fullfilenametxt = fullfile(filePath,filenametxt);
varnames = {'Datetime','incIni[°]','taxaIncIni[°/s]','incMax[°]','taxaIncEns[°/s]','sentidoInc','nStops','tEspera[s]','tCapt[s]','tCaptInterval[s]','erroIncIni[°]','erroIncStop[°]','erroIncTotal[°]','clockPrepIni[s]','clockVid[s]','clockImages[s]','clockStop[s]','clockEnsaio[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
T = table(datetime,incIni,taxaIncIni,incMax,taxaInc,sentidoInc,nStops,tEspera,tCapt,tCaptInterval,round(erroIncIni,4),round(erroIncEns,4),round(erroIncTotal,4),round(clockPrepIni,2),round(mean(clockVid),2),round(mean(clockImages),2),round(mean(clockStop),2),round(clockEnsaio,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file

pause(5);
end

