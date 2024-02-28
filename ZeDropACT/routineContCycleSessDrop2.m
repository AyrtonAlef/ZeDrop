function [volSerTotal] = routineContCycleSessDrop2(vidAdap,vidID,vidFmt,imFmt,smSeringa,passoMot,passoSer,wMin,framerate,limTotalFrames,volSerTotal)
%ROUTINECONTCYCLESESSDROP Execution of continuous sessile drop experiment 
% cycles (needle in). Includes drop advancing and receding.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execution of continuous sessile drop experiment cycles (needle in). 
%   Includes drop advancing and receding. Images are captured concurrently
%   with the continuous increase/decrease of drop volume. After each cycle,
%   the initial drop is reformed to start a new cycle and the maximum
%   volume to be reached by the drop is increased. The routine requires the
%   user to define the initial drop volume (volGotaStart), the initial drop
%   volume before advancing (volGotaIni), the maximum volume reached by the
%   drop after the execution of the last cycle, the flow rate during 
%   advancing (taxaDepAdv), the flow rate during receding (taxaDepRec), the
%   waiting time for drop equilibrium before advancing (tEsperaAdv), the 
%   number of images captured during advancing (nImagesAdv), the waiting 
%   time for drop equilibrium during receding (tEsperaRec), the number of 
%   images captured during receding (nImagesRec), and the number of 
%   measurement cycles (nCycles).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% vidAdap - Character vector that specifies the name of the adaptor used to communicate with the device
% vidID - Numeric scalar value that identifies a particle device available through the specified adaptor 
% vidFmt - Character vector that specifies a particular video format supported by the device
% imFmt - Format of the captured image 
% smSerigna - Syringe motor control
% passoMot - Number of steps required for one revolution in the motor [steps/rev]
% passoSer - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
% wMin - Minimum stable rotation speed (RPM) in which the motor can reach. Corresponds to a sm.RPM = 2 (Microstep)
% framerate - Actual camera frame rate
% volSerTotal - Total liquid volume in the syringe in uL
%   OUTPUT
% volSerTotal - Total liquid volume in the syringe in uL

fprintf('------- SESSILE DROP (NEEDLE IN) - CYCLES OF CONTINUOUS EXPERIMENTS (ADVANCING AND RECEDING) --------\n');
% Experiment parameters (entered by the user)
fprintf('Experiment parameters: \n');
volGotaStart = input("- Initial drop volume (uL): "); %Initial drop volume
%volGotaVarAdvPrep = input("- Variação de volume para preparação da gota para o avanço (uL): "); %Volume variation for drop preparation for advancing
volGotaIni = input("- Initital drop volume before advancing (uL): "); %Initital drop volume before advancing
%volGotaIni = volGotaStart + volGotaVarAdvPrep; %Initial drop volume before the execution of a cycle 
taxaDepAdvIni = 0.1; %Flow rate in uL/s during preparation for advancing
while 1
    volGotaMax = input("- Maximum drop volume (uL): "); %Maximum drop volume reached during the run of the last cycle
    if volGotaMax <= volGotaIni
        fprintf('Invalid maximum drop volume. Please enter a larger volume. \n');
    else
        break;
    end
end
while 1
    taxaDepAdv = input("- Flow rate during advancing (uL/s) (max. 1 uL/s): "); %Flow rate in uL/s during advancing
    if taxaDepAdv > 1 %Check that the flow rate does not exceed the maximum rotation speed limit
        fprintf("Unable to continue the routine. Please lower the flow rate. \n")
    else
        break;
    end
end
%volGotaVarRecPrep = input("- Volume vairation for drop preparation for receding (uL): "); %Volume vairation for drop preparation for receding
%taxaDepRecPrep = 0.1; %Flow rate in uL/s during drop preparation for receding
while 1
    taxaDepRec = input("- Flow rate during receding (uL/s) (max. 1 uL/s): "); %Flow rate in uL/s during receding
    if taxaDepRec > 1 %Check that the flow rate does not exceed the maximum rotation speed limit
        fprintf("Unable to continue the routine. Please lower the flow rate. \n")
    else
        break;
    end
end
tEsperaAdv = input("- Waiting time for drop equilibrium before advancing (s): "); %Waiting time for drop equilibrium before advancing
nImagesAdv = input("- Number of images captured during advancing: "); %Number of images captured during advancing
tEsperaRec = input("- Waiting time for drop equilibrium during receding (s): "); %Waiting time for drop equilibrium during receding
nImagesRec = input("- Number of images captured during receding: "); %Number of images captured during receding
nCycles = input("- Number of measurement cycles (advancing and receding): "); %Number of measurement cycles (advancing and receding)

% Experimental parameters (calculated)
volGotaVarAdvIni = volGotaIni - volGotaStart; %Initial drop volume variation (increase) before advancing 
volGotaIncCycle = (volGotaMax-volGotaIni)/nCycles; %Maximum volume variation between cycles
%volGotaVarCycle = linspace(volGotaIncCycle,(volGotaMax-volGotaIni),nCycles); %Volume variation that will be executed in each cycle
volGotaVarAdvCycle = linspace(volGotaIncCycle,(volGotaMax-volGotaIni),nCycles); %Volume variation that will be executed during advancing in each cycle
volGotaVarRecCycle = volGotaVarAdvCycle + volGotaIni; %Volume variation that will be executed during receding in each cycle
volGotaMaxCycle = volGotaIni + volGotaVarAdvCycle; %Maximum drop volume in each cycle
%volGotaRecPrepCycle = volGotaMaxCycle + volGotaVarRecPrep; %Maximum preparation drop volume for receding in each cycle

tEnsaioAdv = volGotaVarAdvCycle/taxaDepAdv; %Total experiment time in s during advancing in each cycle
totalFrameAdv = ceil(framerate*tEnsaioAdv); %Total number of frames acquired during advancing in each cycle
frameIntervalAdv = round(totalFrameAdv/nImagesAdv); %Frame interval during advancing in each cycle
tCaptIntervalAdv = frameIntervalAdv/framerate; %Time interval between captured frames during advancing in each cycle in s
volCaptIntervalAdv = taxaDepAdv*tCaptIntervalAdv; %Volume variation between captured frames during advancing in each cycle in uL
tEnsaioRec = volGotaVarRecCycle/taxaDepRec; %Total experiment time in s during receding in each cycle
totalFrameRec = ceil(framerate*tEnsaioRec); %Total number of frames acquired during receding in each cycle
frameIntervalRec = round(totalFrameRec/nImagesRec);  %Frame interval during receding in each cycle
tCaptIntervalRec = frameIntervalRec/framerate; %Time interval between captured frames during receding in each cycle in s
volCaptIntervalRec = taxaDepRec*tCaptIntervalRec; %Volume variation between captured frames during receding in each cycle in uL

erroVolAdvStart = 0; %Error in the volume during the initial drop formation
erroVolAdv = zeros(1,nCycles); %Error in the accumulated volume in each cycle during advancing
erroVolRec = zeros(1,nCycles); %Error in the accumulated volume in each cycle during receding
erroVolCycle = zeros(1,nCycles); %Error in the accumulated volume after the end of each cycle
erroVolTotal = 0; %Error in the total drop volume after the experiment
clockAdvStart = 0; %Tempo gasto durante a formação da gota em cada ciclo
clockAdvCycle = zeros(1,nCycles); %Time spent during advancing in each cycle
clockRecCycle = zeros(1,nCycles); %Time spent during receding in each cycle
clockCycle = zeros(1,nCycles); %Time spent during each cycle
clockEnsaio = 0; %Time spent during the execution of the entire experiment

if volGotaVarAdvCycle(end) >= volSerTotal  %Check that there is enough volume in the syringe to perform the procedure
    fprintf("Unable to complete the procedure. Please increase the liquid volume in the syringe. \n")
else
    %- Drive motor settings
    %nPassosvolAdvPrep = vol2passo(volGotaVarAdvPrep,passoMot,passoSer); %No of steps corresponding to the volume released during preparation for advancing
    %nPassosvolRecPrep = vol2passo(volGotaVarRecPrep,passoMot,passoSer); %No of steps corresponding to the volume released during preparation for receding
    nPassosvolAdvStart = vol2passo(volGotaVarAdvIni,passoMot,passoSer); %No of steps corresponding to the initial drop formation
    nPassosvolAdvIni = vol2passo(volGotaVarAdvIni,passoMot,passoSer); %No of steps corresponding to the volume released during preparation for advancing
    %nPassosvolCycle = zeros(1,nCycles);
    nPassosvolAdvCycle = zeros(1,nCycles);
    nPassosvolRecCycle = zeros(1,nCycles);

    %- Image capture settings
    fprintf('Path selection to export video and images... \n');
    filePath = uigetdir('C:\'); %Path selection to export video and images
    name = input('Enter a name for the video and image series: ','s');
    vidDiskFmt = '.avi'; %Video format ('.avi','.mp4','.mj2')
    %imFmt = '.tif'; %Images format ('.jpeg','.png','.tif')
    vidComp = 'Motion JPEG AVI'; %Video compression ('Archival','Motion JPEG AVI','Motion JPEG 2000','MPEG-4','Uncompressed AVI','Indexed AVI' e 'Grayscale AVI')
    vidQual = 90; %Video quality (related with the compression ratio)
    
    %- Begin of experiment
    fprintf('------------- EXECUTION OF MEASUREMENT CYCLES (ADVANCING AND RECEDING) --------------\n');
    fprintf('Press any key to proceed with the cycle experiment... \n');
    pause
    iEns = 1;
    for k = 1:nCycles
        %- Begin of the measurement cycle
        tic
        fprintf('----------------------------------------------------------------------------\n');
        fprintf('Execution of measurement %d of %d \n',k,nCycles);
        close all %Close all open windows
        %delete(vidobj) %Delete any video objects that may have been created
        clear vidobj vidsrc %Clear video related objects
        clockIni= toc;
        clockEnsaio = clockEnsaio + clockIni; %Update total advancing time
        
        %1) ADVANCING ROUTINE
        fprintf('Execution of advancing routine... \n');
        if k ~= 1 %Reformation of initital drop
            tic
            smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
            fprintf('Please wait, initial drop formation...\n');
            move(smSeringa,floor(nPassosvolAdvStart)); %Drive motor connected to the syringe
            release(smSeringa); %Motor release
            fprintf('Initial drop formation completed!\n');
            erroVolAdvStart = -passo2vol(nPassosvolAdvStart- floor(nPassosvolAdvStart),passoMot,passoSer); %Error in the volume of initital drop formation (uL)
            erroVolAdv(k) = erroVolAdv(k) + erroVolAdvStart;
            volSerTotal = volSerTotal - volGotaStart; %Update total liquid volume in the syringe
            pause(2); %Wait for the execution of the next step
            clockAdvStart = toc; %Time for initial drop formation in s
            clockAdvCycle(k) = clockAdvCycle(k) + clockAdvStart; %Update total advancing time
        end
        %- Initial drop increase (preparation for advancing)
        tic
        smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
        %fprintf('Please wait, inital drop increase (preparation for advancing)...\n');
        fprintf('Please wait, inital drop increase ...\n');
        wMotorEns = uLs2RPM(taxaDepAdvIni,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
        tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
        for i = 1:floor(nPassosvolAdvIni) %Initial drop increase
            move(smSeringa,1); %Executes 1 rotation step on the motor
            pause(tadpasso);
        end
        release(smSeringa); %Motor release
        fprintf('Initial drop increase completed! \n');
        %erroVolAdvPrep1 = -passo2vol(nPassosvolAdvPrep - floor(nPassosvolAdvPrep),passoMot,passoSer); %Error in the initial drop volume formed before advancing (uL)
        %erroVolAdv(k) = erroVolAdv(k) + erroVolAdvPrep1;
        erroVolAdvIni = -passo2vol(nPassosvolAdvIni - floor(nPassosvolAdvIni),passoMot,passoSer); %Error in the initial drop volume formed before advancing (uL)
        erroVolAdv(k) = erroVolAdv(k) + erroVolAdvIni;
        %volSerTotal = volSerTotal - volGotaVarAdvPrep; %Update total liquid volume in the syringe
        volSerTotal = volSerTotal - volGotaVarAdvIni; %Update total liquid volume in the syringe
        %clockAdvPrep1 = toc; %Tempo para o aumento inicial da gota em s antes do avanço
        %clockAdvCycle(k) = clockAdvCycle(k) + clockAdvPrep1; %Atualização do tempo total de avanço
        clockAdvIni = toc;  %Time for initial drop increase before advancing in s
        clockAdvCycle(k) = clockAdvCycle(k) + clockAdvIni; %Update total advancing time
    
        %- Waiting time for drop equilibrium
        tic
        fprintf('Please wait, waiting time for drop reach equilibrium...\n');
        pause(tEsperaAdv); %Waiting time for drop reach equilibrium
        fprintf('Waiting time completed! \n');
        clockAdvEspera = toc; %Waiting time for drop reach equilibrium
        clockAdvCycle(k) = clockAdvCycle(k) + clockAdvEspera; %Update total advancing time
    
        %- Increase of drop volume
        tic
        %-- Create video object
        [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

        %- Update of image names
        fileName = strcat(name,"_Avanco_",num2str(k));
        fullvidName = strcat(fileName,vidDiskFmt);
        fullfilevidPath = fullfile(filePath,fullvidName);
        diskLogger = VideoWriter(fullfilevidPath,vidComp);
        diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
        diskLogger.Quality = vidQual; %Definition of the video quality
        vidobj.DiskLogger = diskLogger;
        if totalFrameAdv(k) <= limTotalFrames %Continuous frame capture during video recording
            vidobj.FramesPerTrigger = totalFrameAdv(k); %Total frames to acquire
            vidobj.FrameGrabInterval = 1; %Frame interval
        else %Capturing spaced frames during video recording
            vidobj.FramesPerTrigger = nImagesAdv; %Total frames to acquire
            vidobj.FrameGrabInterval = frameIntervalAdv(k); %Frame interval
        end
        clockVidCreateAdv = toc;
        clockAdvCycle(k) = clockAdvCycle(k) + clockVidCreateAdv; %Update total advancing time

        %-- Start video recording
        tic
        fprintf('Start video recording. \n');
        start(vidobj) %Start video recording
        %-- Continuous increase in drop volume
        %nPassosvolCycle(k) = vol2passo(volGotaVarCycle(k),passoMot,passoSer); %No of steps corresponding to the volume released/withdrawal during advancing or receding in each cycle
        nPassosvolAdvCycle(k) = vol2passo(volGotaVarAdvCycle(k),passoMot,passoSer); %No of steps corresponding to the volume released/withdrawal during advancing or receding in each cycle
        fprintf('Please wait, continuous increase in drop volume... \n');
        if taxaDepAdv <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
            smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
            wMotorEns = uLs2RPM(taxaDepAdv,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
            tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed   
            for i = 1:floor(nPassosvolAdvCycle(k)) %Continuous drop volume release
                move(smSeringa,1); %Executes 1 rotation step on the motor
                pause(tadpasso);
            end
        elseif (0.32 < taxaDepAdv) && (taxaDepAdv <= 1)
            wRPM = fuLs2RPM(taxaDepAdv); %Calculation of the necessary rotational speed to be sent to the motor in RPM
            taxaDepAdv = fRPM2uLs(ceil(wRPM)); %Update flow rate
            % Update variables
            tEnsaioAdv(k) = volGotaVarAdvCycle(k)/taxaDepAdv;
            totalFrameAdv(k) = ceil(framerate*tEnsaioAdv(k));
            frameIntervalAdv(k) = round(totalFrameAdv(k)/nImagesAdv);
            tCaptIntervalAdv(k) = frameIntervalAdv(k)/framerate;
            volCaptIntervalAdv(k) = taxaDepAdv*tCaptIntervalAdv(k);
            % Syringe plunger movement
            smSeringa.RPM = round(wRPM);
            move(smSeringa,floor(nPassosvolAdvCycle(k))); %Motor movement
        end
        release(smSeringa); %Motor release
        fprintf('Increase of drop volume completed! \n');
        erroVolDepAdv = passo2vol(nPassosvolAdvCycle(k) - floor(nPassosvolAdvCycle(k)),passoMot,passoSer); %Error in the volume released during the experiment
        erroVolAdv(k) = erroVolAdv(k) + erroVolDepAdv; %Update error on total drop volume
        volSerTotal = volSerTotal - volGotaVarAdvCycle(k); %Update total liquid volume in the syringe
        %-- Finishing video recording
        wait(vidobj) %Wait until all frames are captured
        stop(vidobj);
        %--- Wait for all frames to be written to disk
        while vidobj.FramesAcquired ~= vidobj.DiskLoggerFrameCount
            pause(.1);
        end
        fprintf('Video recording completed! \n');
        clockAdv = toc;
        clockAdvCycle(k) = clockAdvCycle(k) + clockAdv; %Update total advancing time
        
        %-- Showing an example image
        tic
        close all %Close all open windows
        frame = getsnapshot(vidobj); %Capture an image
        imshow(frame); %Show image
        movegui(gcf,'center') %Centralize image window
        clockExAdv = toc;
        clockAdvCycle(k) = clockAdvCycle(k) + clockExAdv; %Update total advancing time
        
        %- Converting video to images
        tic
        fprintf('Please wait, converting video to images... \n');
        %iEns = 1;
        if totalFrameAdv(k) <= limTotalFrames %Continuous frame capture during video recording
            fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
            vidFile = VideoReader(fullfilevidPath); %Read Video File
            iAdv = 1;
            for m = 1: round(frameIntervalAdv(k)): totalFrameAdv(k) %Getting a image from every video second
                frame = read(vidFile,m);
                framegray = rgb2gray(frame);
                if iEns < 10 %Modifying the file name prefix to arrange files in ascending volume order
                    prefixo = "00";
                elseif (iEns >= 10) && (iEns < 100)
                    prefixo = "0";
                elseif (iEns >= 100) && (iEns < 1000)
                    prefixo = "";
                end
                fullimName = strcat(prefixo,num2str(iEns),'_',fileName,'_V(uL)_',num2str((volGotaIni + iAdv*volCaptIntervalAdv(k)),'%.2f'),imFmt);
                fullfileimPath = fullfile(filePath,fullimName);
                imwrite(framegray,fullfileimPath);
                iAdv = iAdv + 1;
                iEns = iEns + 1; %Interval of image name
            end
        else %Capturing spaced frames during video recording
            fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
            vidFile = VideoReader(fullfilevidPath); %Read Video File
            for m = 1: nImagesAdv %Getting a image from every video frame
                frame = read(vidFile,k);
                framegray = rgb2gray(frame);
                if iEns < 10 %Modifying the file name prefix to arrange files in ascending volume order
                    prefixo = "00";
                elseif (iEns >= 10) && (iEns < 100)
                    prefixo = "0";
                elseif (iEns >= 100) && (iEns < 1000)
                    prefixo = "";
                end
                fullimName = strcat(prefixo,num2str(iEns),'_',fileName,'_V(uL)_',num2str((volGotaIni + m*volCaptIntervalAdv(k)),'%.2f'),imFmt);
                fullfileimPath = fullfile(filePath,fullimName);
                imwrite(framegray,fullfileimPath);
                iEns = iEns + 1; %Interval of image name
            end
        end
        delete(fullfilevidPath) %Delete video
        fprintf('Conversion of video to images completed! \n');
        clockImagesAdv = toc;
        clockAdvCycle(k) = clockAdvCycle(k) + clockImagesAdv; %Update total advancing time
        
        %- Deleting video object
        tic
        delete(vidobj)
        clear vidobj vidsrc
        clockExcVideoAdv = toc;
        clockAdvCycle(k) = clockAdvCycle(k) + clockExcVideoAdv; %Update total advancing time
        
        %{
        %- Increase of drop volume (preparation for receding)
        tic
        smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
        fprintf('Please wait, increasing drop volume (preparation for receding)... \n');
        wMotorEns = uLs2RPM(taxaDepRecPrep,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
        tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
        for i = 1:floor(nPassosvolRecPrep) %Increase of drop volume
            move(smSeringa,1); %Executes 1 rotation step on the motor
            pause(tadpasso);
        end
        release(smSeringa); %Motor release
        fprintf('Increase of drop volume completed! \n');
        erroVolRecPrep1 = -passo2vol(nPassosvolRecPrep - floor(nPassosvolRecPrep),passoMot,passoSer); %Error in the preparation drop volume before receding (uL)
        erroVolAdv(k) = erroVolAdv(k) + erroVolRecPrep1;
        volSerTotal = volSerTotal - volGotaVarRecPrep; %Update total liquid volume in the syringe
        clockRecPrep1 = toc; %Time of preparaing drop volume before receding in s
        clockAdvCycle(k) = clockAdvCycle(k) + clockRecPrep1; %Update total advancing time
        %}

        %2) RECEDING ROUTINE
        tic
        pause(1)
        fprintf('Execution of receding routine... \n');
        close all %Close all open windows
        %delete(vidobj) %Delete any video objects that may have been created
        clear vidobj vidsrc %Clear video related objects
        clockRecuo1= toc;
        clockRecCycle(k) = clockRecCycle(k) + clockRecuo1; %Update total receding time
        
        %{
        %- Initial withdrawal of drop volume (preparation for receding)
        tic
        fprintf('Please wait, initial withdrawal of drop volume (preparation for receding)...\n');
        smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
        wMotorEns = uLs2RPM(taxaDepRecPrep,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
        tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
        for i = 1:floor(nPassosvolRecPrep) %Initital withdrawal of drop volume
            move(smSeringa,-1); %Executes 1 rotation step on the motor
            pause(tadpasso);
        end
        release(smSeringa); %Motor release
        fprintf('Initial withdrawal of drop volume completed! \n');
        erroVolRecPrep2 = passo2vol(nPassosvolRecPrep - floor(nPassosvolRecPrep),passoMot,passoSer); %Error in the initial drop volume before receding (uL)
        erroVolRec(k) = erroVolRec(k) + erroVolRecPrep2;
        volSerTotal = volSerTotal + volGotaVarRecPrep; %Update total liquid volume in the syringe
        clockRecPrep2 = toc; %Time for initial drop receding in s before receding
        clockRecCycle(k) = clockRecCycle(k) + clockRecPrep2; %Update total receding time
        %}
        %%{
        %- Waiting time for drop reach equilibrium
        tic
        fprintf('Please wait, waiting time for drop reach equilibrium...\n');
        pause(tEsperaRec); %Waiting time for drop reach equilibrium
        fprintf('Waiting time completed! \n');
        clockRecEspera = toc; %Waiting time for drop reach equilibrium
        clockRecCycle(k) = clockRecCycle(k) + clockRecEspera; %Update total receding time
        %}

        %- Decrease of drop volume
        tic
        %-- Create video object
        [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

        %- Update of image names
        fileName = strcat(name,"_Recuo_",num2str(k));
        fullvidName = strcat(fileName,vidDiskFmt);
        fullfilevidPath = fullfile(filePath,fullvidName);
        diskLogger = VideoWriter(fullfilevidPath,vidComp);
        diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
        diskLogger.Quality = vidQual; %Definition of the video quality
        vidobj.DiskLogger = diskLogger;
        if totalFrameRec(k) <= limTotalFrames %Continuous frame capture during video recording
            vidobj.FramesPerTrigger = totalFrameRec(k); %Total frames to acquire
            vidobj.FrameGrabInterval = 1; %Frame interval
        else %Capturing spaced frames during video recording
            vidobj.FramesPerTrigger = nImagesRec; %Total frames to acquire
            vidobj.FrameGrabInterval = frameIntervalRec(k); %Frame interval
        end
        clockVidCreateRec = toc;
        clockRecCycle(k) = clockRecCycle(k) + clockVidCreateRec; %Update total receding time

        %-- Start video recording
        tic
        fprintf('Start video recording. \n');
        start(vidobj) %Start video recording
        %-- Continuous decrease in drop volume
        fprintf('Please wait, continuous decrease in drop volume... \n');
        if taxaDepRec <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
            smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
            wMotorEns = uLs2RPM(taxaDepRec,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
            tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed   
            for i = 1:floor(nPassosvolRecCycle(k)) %Continuous drop volume withdrawal
                move(smSeringa,-1); %Executes 1 rotation step on the motor
                pause(tadpasso);
            end
        elseif (0.32 < taxaDepRec) && (taxaDepRec <= 1)
            wRPM = fuLs2RPM(taxaDepRec); %Calculation of the necessary rotational speed to be sent to the motor in RPM
            taxaDepRec = fRPM2uLs(ceil(wRPM)); %Update flow rate
            % Update variables
            tEnsaioRec(k) = volGotaVarRecCycle(k)/taxaDepRec;
            totalFrameRec(k) = ceil(framerate*tEnsaioRec(k));
            frameIntervalRec(k) = round(totalFrameRec(k)/nImagesRec);
            tCaptIntervalRec(k) = frameIntervalRec(k)/framerate;
            volCaptIntervalRec(k) = taxaDepRec*tCaptIntervalRec(k);
            % Syringe plunger movement
            smSeringa.RPM = round(wRPM);
            move(smSeringa,-floor(nPassosvolRecCycle(k))); %Motor movement
        end
        release(smSeringa); %Motor release
        fprintf('Continuous decrease in drop volume completed! \n');
        erroVolAdqRec = passo2vol(nPassosvolRecCycle(k) - floor(nPassosvolRecCycle(k)),passoMot,passoSer); %Error in the volume withdrawal during the experiment
        erroVolRec(k) = erroVolRec(k) + erroVolAdqRec; %Update of total error in the volume during receding
        volSerTotal = volSerTotal + volGotaVarRecCycle(k); %Update total liquid volume in the syringe
        %-- Finishing video recording
        wait(vidobj) %Wait until all frames are captured
        stop(vidobj);
        %--- Wait for all frames to be written to disk
        while vidobj.FramesAcquired ~= vidobj.DiskLoggerFrameCount
            pause(.1);
        end
        fprintf('Video recording completed! \n');
        clockRec = toc;
        clockRecCycle(k) = clockRecCycle(k) + clockRec; %Update total receding time
      
        %-- Showing an example image
        tic
        close all %Close all open windows
        frame = getsnapshot(vidobj); %Capture an image
        imshow(frame); %Show image
        movegui(gcf,'center') %Centralize image window
        clockExRec = toc;
        clockRecCycle(k) = clockRecCycle(k) + clockExRec; %Update total receding time
   
        %- Converting video to images
        tic
        fprintf('Please wait, converting video to images... \n');
        if totalFrameRec(k) <= limTotalFrames %Continuous frame capture during video recording
            fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
            vidFile = VideoReader(fullfilevidPath); %Read Video File
            iRec = 1;
            for m = 1: round(frameIntervalRec(k)): totalFrameRec(k) %Getting a image from every video second
                frame = read(vidFile,m);
                framegray = rgb2gray(frame);
                if iEns < 10 %Modifying the file name prefix to arrange files in ascending volume order
                    prefixo = "00";
                elseif (iEns >= 10) && (iEns < 100)
                    prefixo = "0";
                elseif (iEns >= 100) && (iEns < 1000)
                    prefixo = "";
                end
                fullimName = strcat(prefixo,num2str(iEns),'_',fileName,'_V(uL)_',num2str((volGotaMaxCycle(k) - iRec*volCaptIntervalRec(k)),'%.2f'),imFmt);
                fullfileimPath = fullfile(filePath,fullimName);
                imwrite(framegray,fullfileimPath);
                iRec = iRec + 1; 
                iEns = iEns + 1; %Interval of image name
            end
        else %Capturing spaced frames during video recording
            fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
            vidFile = VideoReader(fullfilevidPath); %Read Video File
            for m = 1: nImagesRec %Getting a image from every video frame
                frame = read(vidFile,m);
                framegray = rgb2gray(frame);
                if iEns < 10 %Modifying the file name prefix to arrange files in ascending volume order
                    prefixo = "00";
                elseif (iEns >= 10) && (iEns < 100)
                    prefixo = "0";
                elseif (iEns >= 100) && (iEns < 1000)
                    prefixo = "";
                end
                fullimName = strcat(prefixo,num2str(iEns),'_',fileName,'_V(uL)_',num2str((volGotaMaxCycle(k) - m*volCaptIntervalRec(k)),'%.2f'),imFmt);
                fullfileimPath = fullfile(filePath,fullimName);
                imwrite(framegray,fullfileimPath);
                iEns = iEns + 1; %Interval of image name
            end
        end
        delete(fullfilevidPath) %Delete video
        fprintf('Conversion of video to images completed! \n');
        clockImagesRec = toc;
        clockRecCycle(k) = clockRecCycle(k) + clockImagesRec; %Update total receding time

        %- Deleting video object
        tic
        delete(vidobj)
        clear vidobj vidsrc
        clockExcVideoRec = toc;
        clockRecCycle(k) = clockRecCycle(k) + clockExcVideoRec; %Update total receding time
        
        %{
        %- Final withdrawal of drop volume (preparation for the next cycle)
        tic
        fprintf('Please wait, final withdrawal of drop volume (preparation for the next cycle)...\n');
        smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
        wMotorEns = uLs2RPM(taxaDepAdvPrep,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
        tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
        for i = 1:floor(nPassosvolAdvPrep) %Final withdrawal of drop volume
            move(smSeringa,-1); %Executes 1 rotation step on the motor
            pause(tadpasso);
        end
        release(smSeringa); %Motor release
        fprintf('Final withdrawal of drop volume completed! \n');
        erroVolAdvPrep2 = passo2vol(nPassosvolAdvPrep - floor(nPassosvolAdvPrep),passoMot,passoSer); %Error in the final drop volume (uL)
        erroVolRec(k) = erroVolRec(k) + erroVolAdvPrep2;
        volSerTotal = volSerTotal + volGotaVarRecPrep; %Update total liquid volume in the syringe
        clockRecPrep2 = toc; %Time for the final withdrawal of drop volume in s
        clockRecCycle(k) = clockRecCycle(k) + clockRecPrep2; %Update total receding time
        %}

        %Calculation of time and volume error during each cycle
        erroVolCycle(k) = erroVolAdv(k) + erroVolRec(k);
        clockCycle(k) = clockAdvCycle(k) + clockRecCycle(k);

        %- Export/save description of each cycle
        filenametxt = strcat('Description_Cycle_',num2str(k),'.txt');
        fullfilenametxt = fullfile(filePath,filenametxt);
        varnames = {'Datetime','volGotaStart[uL]','volGotaIni[uL]','taxaDepAdvIni[uL]','volGotaMax[uL]','taxaDepAdv[uL/s]','tEsperaAdv[s]','tEnsaioAdv[s]','nImagesAdv','totalFrameAdv','frameIntervalAdv','tCaptIntervalAdv[s]','volCaptIntervalAdv[uL]','taxaDepRec[uL/s]','tEsperaRec[s]','tEnsaioRec[s]','nImagesRec','totalFrameRec','frameIntervalRec','tCaptIntervalRec[s]','volCaptIntervalRec[uL]','erroVolAdvStart[uL]','erroVolAdvIni[uL]','erroVolDepAdv[uL]''erroVolAdv[uL]','erroVolAdqRec[uL]','erroVolRec[uL]','erroVolCycle[uL]','clockAdvStart[s]','clockAdvIni[s]','clockAdv[s]','clockImagesAdv[s]','clockAdvCycle[s]','clockRec[s]','clockRecCycle[s]','clockCycle[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
        T = table(datetime,volGotaStart,volGotaIni,taxaDepAdvIni,volGotaMaxCycle(k),taxaDepAdv,tEsperaAdv,tEnsaioAdv(k),nImagesAdv,totalFrameAdv(k),frameIntervalAdv(k),tCaptIntervalAdv(k),volCaptIntervalAdv(k),taxaDepRec,tEsperaRec,tEnsaioRec(k),nImagesRec,totalFrameRec(k),frameIntervalRec(k),tCaptIntervalRec(k),volCaptIntervalRec(k),round(erroVolAdvStart(k),4),round(erroVolAdvIni(k),4),round(erroVolDepAdv(k),4),round(erroVolAdv(k),4),round(erroVolAdqRec(k),4),round(erroVolRec(k),4),round(erroVolCycle(k),4),round(clockAdvStart,2),round(clockAdvIni,2),round(clockAdv,2),round(clockImagesAdv,2),round(clockAdvCycle(k),2),round(clockRec,2),round(clockRecCycle(k),2),round(clockCycle(k),2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
        writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file
        pause(1);
    end   
    % Calculation of total time and volume error
    erroVolTotal = sum(erroVolCycle); %Total error in final drop volume including advancing and receding
    clockEnsaio = sum(clockCycle); %Total time spent on the entire routine including advancing and receding

    %- Export/save experiment description
    filenametxt = "Description.txt";
    fullfilenametxt = fullfile(filePath,filenametxt);
    varnames = {'Datetime','volGotaStart[uL]','volGotaIni[uL]','taxaDepAdvIni[uL]','volGotaMax[uL]','taxaDepAdv[uL/s]','tEsperaAdv[s]','nImagesAdv','taxaDepRec[uL/s]','tEsperaRec[s]','nImagesRec','nCycles','erroVolTotal[uL]','clockEnsaio[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
    T = table(datetime,volGotaStart,volGotaIni,taxaDepAdvIni,volGotaMax,taxaDepAdv,tEsperaAdv,nImagesAdv,taxaDepRec,tEsperaRec,nImagesRec,nCycles,round(erroVolTotal,4),round(clockEnsaio,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
    writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file
    pause(5);
    end
end

