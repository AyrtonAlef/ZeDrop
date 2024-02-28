function [volSerTotal] = routineIntSessDrop(vidAdap,vidID,vidFmt,imFmt,smSeringa,smSeringaZ,passoMot,passoSer,passoFuso,wMin,framerate,volSerTotal)
%ROUTINEINTSESSDROP Execution of a intermittent sessile drop experiment 
% (needle out). Includes advancing and receding routines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execution of a intermittent sessile drop experiment (needle out). 
%   Includes advancing and receding routines. Stops during advancing and 
%   receding are performed to remove the needle from the drop and capture 
%   images. Initially, the drop volume is increased until the indicated 
%   volume is reached (volGotaAdvIni) to carry out the advancing routine. 
%   After completing the advancing routine, a new increase in the drop 
%   volume is carried out and subsequent decrease of the drop volume, at 
%   high and low flow rate, until the receding drop volume is reached to 
%   perform the receding routine (volGotaRecIni) . The routine requires the
%   user to define the initial drop volume (volGotaStart), initial drop 
%   volume before advancing (volGotaAdvIni), maximum drop volume after 
%   advancing (volGotaAdvMax), the flow rate during advancing (taxaDepAdv),
%   the syringe displacement (dispZSer), the waiting time for drop 
%   equilibrium during advancing (tEsperaAdv), the capture time during a 
%   stop in advancing (tCaptAdv), the number of captured images during 
%   advancing (nImagesAdv), the number of stops during advancing 
%   (nStopsAdv), the preparation drop volume for receding (volGotaRecPrep),
%   the initial drop volume before receding (volGotaRecIni), the flow rate 
%   during receding (taxaDepRec), the waiting time for drop equilibrium 
%   during receding (tEsperaRec), the capture time during a stop in 
%   receding (tCaptRec), the number of captured images during receding 
%   (nImagesRec), number of stops during receding (nStopsRec). The routine 
%   allows choosing the directory where the images will be exported, as 
%   well as defining the name of the images. At the end, a 
%   text file is exported with all relevant parametes of the experiment. 
%   Prior to start the experiment, the needle is required to be inserted 
%   into the sessile drop. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% vidAdap - Character vector that specifies the name of the adaptor used to communicate with the device
% vidID - Numeric scalar value that identifies a particle device available through the specified adaptor 
% vidFmt - Character vector that specifies a particular video format supported by the device
% imFmt - Format of the captured image 
% smSerigna - Syringe motor control
% smSeringaZ - Syringe displacement motor control
% passoMot - Number of steps required for one revolution in the motor [steps/rev]
% passoSer - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
% passoFuso - Spindle pitch [mm/rev]
% wMin - Minimum stable rotation speed (RPM) in which the motor can reach. Corresponds to a sm.RPM = 2 (Microstep)
% framerate - Actual camera frame rate
% volSerTotal - Total liquid volume in the syringe in uL
%   OUTPUT
% volSerTotal - Total liquid volume in the syringe in uL

fprintf('------------------ SESSILE DROP (NEEDLE OUT) - INTERMITTENT EXPERIMENT (ADVANCING AND RECEDING) --------------------\n');
% Experiment parameters (entered by the user)
fprintf('1) Experiment parameters during advancing: \n');
volGotaStart = input("-  Initial drop volume (uL): "); %Initial drop volume
while 1
    volGotaAdvIni = input("- Initial drop volume before advancing (uL): "); %Initial drop volume for the execution of the advancing contact angle measurement protocol 
    if volGotaAdvIni < volGotaStart
        fprintf('Invalid initial drop volume before advancing. Please enter a larger volume. \n');
    else
        break;
    end
end
taxaDepAdvIni = 0.1; %Flow rate in uL/s during the initial drop advance
while 1
    volGotaAdvMax = input("- Maximum drop volume after advancing (uL): "); %Maximum drop volume reached during the execution of the advancing contact angle measurement protocol
    if volGotaAdvMax <= volGotaAdvIni
        fprintf('Invalid maximum drop volume after advancing. Please enter a larger volume. \n');
    else
        break;
    end
end
while 1
    taxaDepAdv = input("- Flow rate during advancing (uL/s) (max. 1 uL/s): "); %Drop flow rate in uL/s during advancing
    if taxaDepAdv > 1 %Check that the flow rate does not exceed the maximum rotation speed
        fprintf("Unable to proceed with the routine. Please lower the flow rate. \n")
    else
        break;
    end
end
dispZSer = input("- Syringe displacement (mm): "); %Syringe displacement during the execution of advancing and receding contact angle measurement protocol. Displacement required for the needle to move out of the camera's field of view
tEsperaAdv = input("- Waiting time for drop equilibrium during advancing (s): "); %Waiting time for drop equilibrium during advancing
tCaptAdv = input("- Capture time during a stop in advancing (s): "); %Capture time during a stop in advancing
nImagesAdv = input("- Number of captured images during a stop in advancing: "); %Number of captured images during a stop in advancing
nStopsAdv = input("-  Number of stops during advancing: "); %Number of stops during advancing

fprintf('2) Experiment parameters during receding: \n');
while 1
    volGotaRecPrep = input("- Preparation drop volume for receding (uL): "); %Preparation drop volume for receding
    if volGotaRecPrep < volGotaAdvMax
        fprintf('Invalid preparation drop volume for receding. Please enter a larger volume. \n');
    else
        break;
    end
end
while 1
    volGotaRecIni = input("- Initial drop volume before receding (uL): "); %Initial drop volume before receding
    if volGotaRecIni > volGotaRecPrep
        fprintf('Invalid initial drop volume before receding. Please enter a smaller volume. \n');
    elseif volGotaRecIni < volGotaAdvIni
        fprintf('Invalid initial drop volume before receding. Please enter a larger volume. \n');
    else
        break;
    end
end
taxaDepRecIni = 0.1; %Initial flow rate before receding in uL/s
while 1
    volGotaRecEnd = input("- Final drop volume after receding (uL): "); %Final drop volume after receding protocol
    if volGotaRecEnd > volGotaRecIni
        fprintf('Invalid final drop volume after receding. Please enter a smaller volume. \n');
    else
        break;
    end
end
while 1
    taxaDepRec = input("- Flow rate during receding (uL/s) (max. 1 uL/s): "); %Flow rate during receding in uL/s
    if taxaDepRec > 1 %Check that the flow rate does not exceed the maximum motor rotational speed limit
        fprintf("Unable to proceed with the routine. Please lower the flow rate. \n")
    else
        break;
    end
end
tEsperaRec = input("- Waiting time for drop equilibrium during receding (s): "); %Waiting time for drop equilibrium during receding
tCaptRec = input("- Capture time during a stop in receding (s): "); %Capture time during a stop in receding
nImagesRec = input("- Number of captured images during a stop in receding: "); %Number of captured images during a stop in receding
nStopsRec = input("- Number of stops during receding: "); %Number of stops during receding

% Experimental parameters (calculated)
volGotaVarAdvIni = volGotaAdvIni - volGotaStart; %Volume variation during intial drop increase before advancing
volGotaVarAdv = (volGotaAdvMax - volGotaAdvIni)/nStopsAdv; %Volume variation in each stop during advancing
volGotaVarRecPrep = volGotaRecPrep - volGotaAdvMax; %Volume variation during preparation (drop increase) before receding
volGotaVarRecIni = abs(volGotaRecPrep - volGotaRecIni); %Volume variation during the initial step of receding (drop decrease)
volGotaVarRecIniAlta = floor(volGotaVarRecIni*0.8); %Volume variation during the initial step of receding (drop decrease) at high flow rate
volGotaVarRecIniBaixa= volGotaVarRecIni - volGotaVarRecIniAlta; %Volume variation during the initial step of receding (drop decrease) at low flow rate
volGotaVarRec = abs(volGotaRecIni - volGotaRecEnd)/nStopsRec; %Volume variation in each stop during receding
tCaptIntervalAdv = tCaptAdv/nImagesAdv; %Acquisition time between images during advancing in s
tCaptIntervalRec = tCaptRec/nImagesRec; %Acquisition time between images during receding in s
totalFrameAdv = ceil(framerate*tCaptAdv); %Number of total frames acquired during advancing
totalFrameRec = ceil(framerate*tCaptRec); %Number of total frames acquired during receding
frameIntervalAdv = round(totalFrameAdv/nImagesAdv); %Frame interval during advancing
frameIntervalRec = round(totalFrameRec/nImagesRec); %Frame interval during receding
erroVolTotalAdv = 0;
erroVolTotalRec = 0;
clockPrepVidAdv = zeros(1,nStopsAdv); %Time spent in video preparation in each stop during advancing
clockAdvGota = zeros(1,nStopsAdv); %Time spent in addition volume release in each stop during advancing
clockDepSerAdv = zeros(1,nStopsAdv); %Time spent in syringe departure from the drop in each stop during advancing
clockEspAdv = zeros(1,nStopsAdv); %Time spent in waiting for drop equilibrium before video recording in each stop during advancing
clockVidAdv = zeros(1,nStopsAdv); %Time spent in video recording in each stop during advancing
clockImagesAdv = zeros(1,nStopsAdv); %Time spent in conversion video to images in each stop during advancing
clockRetSerAdv = zeros(1,nStopsAdv); %Time spent in syringe return to the drop in each stop during advancing
clockStopAdv = zeros(1,nStopsAdv); %Total time spent in each stop during advancing
clockPrepVidRec = zeros(1,nStopsRec); %Time spent in video preparation in each stop during receding
clockRecGota = zeros(1,nStopsRec); %Time spent in addition volume withdrawal in each stop during receding
clockDepSerRec = zeros(1,nStopsRec); %Time spent in syringe departure from the drop in each stop during receding
clockEspRec = zeros(1,nStopsRec); %Time spent in waiting for drop equilibrium before video recording in each stop during receding
clockVidRec = zeros(1,nStopsRec); %Time spent in video recording in each stop during receding
clockImagesRec = zeros(1,nStopsRec); %Time spent in conversion video to images in each stop during receding
clockRetSerRec = zeros(1,nStopsRec); %Time spent in syringe return to the drop in each stop during receding
clockStopRec = zeros(1,nStopsRec); %Total time spent in each stop during receding
clockAvanco = 0; %Total time spent during advancing
clockRecuo = 0; %Total time spent during receding
clockEnsaio = 0; %Total time spent in the experiment
%volCaptInterval = taxaDepGota*tCaptInterval; %Volume variation between images in uL

if volGotaRecPrep >= volSerTotal  %Check that there is enough volume in the syringe to perfrom the procedure
    fprintf("Unable to proceed with the routine. Please increase the liquid volume in the syringe. \n")
else
    %- Configuração de movimentação do motor
    nPassosvolDepAdvIni = vol2passo(volGotaVarAdvIni,passoMot,passoSer); %No of steps corresponding to the released volume during preparation for advancing
    nPassosvolDepAdv = vol2passo(volGotaVarAdv,passoMot,passoSer); %No of steps corresponding to the released volume between each measurement during advancing
    nPassosvolDepRecPrep = vol2passo(volGotaVarRecPrep,passoMot,passoSer); %No of steps corresponding to the volume released during the preparation for receding
    nPassosvolAdqRecIniAlta = vol2passo(volGotaVarRecIniAlta,passoMot,passoSer); %No of steps corresponding to the volume withdrawal before receding in high flow rate
    nPassosvolAdqRecIniBaixa = vol2passo(volGotaVarRecIniBaixa,passoMot,passoSer); %No of steps corresponding to the volume withdrawal before receding in low flow rate
    nPassosvolAdqRec = vol2passo(volGotaVarRec,passoMot,passoSer); %No of steps corresponding to the volume withdrawal in each stop during receding
    
    %- Image capture settings
    fprintf('Path selection to export video and images... \n');
    filePath = uigetdir('C:\'); %Path selection to export video and images
    name = input('Enter a name for the video and image series: ','s');
    vidDiskFmt = '.avi'; %Video format ('.avi','.mp4','.mj2')
    %imFmt = '.tif'; %Images format ('.jpeg','.png','.tif')
    vidComp = 'Motion JPEG AVI'; %Video compression ('Archival','Motion JPEG AVI','Motion JPEG 2000','MPEG-4','Uncompressed AVI','Indexed AVI' e 'Grayscale AVI')
    vidQual = 90; %Video quality (related with the compression ratio)
    
    %1) ADVANCING ROUTINE
    %- Begin of experiment
    fprintf('--------------------- EXECUTION OF DROP ADVANCING ROUTINE ----------------------\n');
    fprintf('Press any key to proceed with the advancing routine... \n');
    pause;
    tic
    close all %Close all open windows
    %delete(vidobj) %Delete any video objects that may have been created
    clear vidobj vidsrc %Clear video related objects
    clockAvanco1= toc;
    clockAvanco = clockAvanco + clockAvanco1; %Update total advancing time
    
    %- Initial drop increase
    tic
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, inital drop increase...\n');
    wMotorEns = uLs2RPM(taxaDepAdvIni,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
    tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
    for i = 1:floor(nPassosvolDepAdvIni) %Initial drop increase
        move(smSeringa,1); %Executes 1 rotation step on the motor
        pause(tadpasso);
    end
    release(smSeringa); %Mototr release
    fprintf('Initial drop increase completed! \n');
    erroVolDepAdvIni = -passo2vol(nPassosvolDepAdvIni - floor(nPassosvolDepAdvIni),passoMot,passoSer); %Error in the initial drop volume formed before advancing (uL)
    erroVolTotalAdv = erroVolTotalAdv + erroVolDepAdvIni;
    volSerTotal = volSerTotal - volGotaVarAdvIni; %Update total liquid volume in the syringe
    clockAdvGotaIni = toc; %Time for initial drop increase in s before advancing
    clockAvanco = clockAvanco + clockAdvGotaIni; %Update total advancing time

    %- Waiting time for drop equilibrium
    tic
    fprintf('Please wait, waiting time for drop reach equilibrium...\n');
    pause(tEsperaAdv); %Waiting time for drop equilibrium
    fprintf('Waiting time completed! \n');
    clockEspAdvGotaIni = toc; %Time spent to wait drop to reach equilibrium before advancing
    clockAvanco = clockAvanco + clockEspAdvGotaIni; %Update total advancing time

    %- Advancing measurement series
    iEns = 1; %Counter for all captured images during advancing
    for k = 1:nStopsAdv
        tic
        %- Create video object
        [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

        %- Update of image names
        fileName = strcat(name,'_Avanco_',num2str(k),'_V(uL)_',num2str((volGotaAdvIni+k*volGotaVarAdv),'%.2f'));
        fullvidName = strcat(fileName,vidDiskFmt);
        fullfilevidPath = fullfile(filePath,fullvidName);
        
        %- Preparation of video object
        diskLogger = VideoWriter(fullfilevidPath,vidComp);
        diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
        diskLogger.Quality = vidQual; %Definition of the video quality
        vidobj.DiskLogger = diskLogger;
        %vidobj.FramesPerTrigger = framerate*tCapt; %Total frames to acquire
        vidobj.FramesPerTrigger = nImagesAdv; %Total frames to acquire
        vidobj.FrameGrabInterval = frameIntervalAdv; %Frame interval

        %- Execution of measurement
        fprintf("------------------------------------------------------------------- \n");
        fprintf("Execution of advancing measurement %d of %d \n",k,nStopsAdv);
        clockPrepVidAdv(k) = toc;
        clockStopAdv(k) = clockStopAdv(k) + clockPrepVidAdv(k); %Update time spent in each stop during advancing
        clockAvanco = clockAvanco + clockPrepVidAdv(k); %Update total advancing time

        %- Additional increase of drop volume
        tic
        fprintf('Please wait, additional increase of drop volume... \n');
        if taxaDepAdv <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
            smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
            wMotorEns = uLs2RPM(taxaDepAdv,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
            tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
                for i = 1:floor(nPassosvolDepAdv) %Additional increase of drop volume
                    move(smSeringa,1); %Executes 1 rotation step on the motor
                    pause(tadpasso);
                end
        elseif (0.32 < taxaDepAdv) && (taxaDepAdv <= 1)
            wRPM = fuLs2RPM(taxaDepAdv); %Calculation of the necessary rotational speed to be sent to the motor in RPM
            smSeringa.RPM = floor(wRPM);
            move(smSeringa,floor(nPassosvolDepAdv)); %Motor movement
        end
        release(smSeringa); %Motor release
        fprintf('Additional increase of drop volume completed! \n');
        erroVolDepAdv = -passo2vol(nPassosvolDepAdv - floor(nPassosvolDepAdv),passoMot,passoSer); %Error in the volume released during advancing (uL)
        erroVolTotalAdv = erroVolTotalAdv + erroVolDepAdv;
        volSerTotal = volSerTotal - volGotaVarAdv; %Update total liquid volume in the syringe
        pause(1);
        clockAdvGota(k) = toc; %Time spent in the increase of drop volume during a stop
        clockStopAdv(k) = clockStopAdv(k) + clockAdvGota(k); %Update time spent each stop during advancing
        clockAvanco = clockAvanco + clockAdvGota(k); %Update total advancing time

        %- Syringe displacement (syringe departure from the drop)
        tic
        fprintf('Please wait, syringe displacement... \n');
        moverSeringa(passoMot,passoFuso,smSeringaZ,"up",dispZSer,"high"); %Syringe movement
        release(smSeringaZ); %Motor release
        fprintf('Syringe displacement completed! \n');
        clockDepSerAdv(k) = toc;
        clockStopAdv(k) = clockStopAdv(k) + clockDepSerAdv(k); %Update time spent each stop during advancing
        clockAvanco = clockAvanco + clockDepSerAdv(k); %Update total advancing time

        %- Waiting time for drop reach equilibrium
        tic
        fprintf('Please wait, waiting time for drop reach equilibrium...\n');
        if clockDispSerAdv1(k) >= tEsperaAdv %Syringe displacement time is longer than waiting time for drop equilibrium
            tEsperaAdvReal(k) = clockDispSerAdv1(k); %Calculation of actual waiting time
        else %Syringe displacement time is less than waiting time for drop equilibrium
            tEsperaAdvReal(k) = tEsperaAdv;
            deltat = tEsperaAdv - clockDispSerAdv1(k);
            pause(deltat)
        end
        %pause(tEsperaAdv);
        fprintf('Waiting time completed! \n');
        clockEspAdv(k) = toc;
        clockStopAdv(k) = clockStopAdv(k) + clockEspAdv(k); %Update time spent each stop during advancing
        clockAvanco = clockAvanco + clockEspAdv(k); %Update total advancing time

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
        clockVidAdv(k) = toc;
        clockStopAdv(k) = clockStopAdv(k) + clockVidAdv(k); %Update time spent each stop during advancing
        clockAvanco = clockAvanco + clockVidAdv(k); %Update total advancing time
        
        %- Showing an example image
        tic
        close all %Close all open windows
        frame = getsnapshot(vidobj); %Capture an image
        imshow(frame); %Show image
        movegui(gcf,'center') %Centralize image window
        set(gcf,'Name',strcat(fileName,'_Example'))
        clockExAdv = toc;
        clockStopAdv(k) = clockStopAdv(k) + clockExAdv; %Update time spent each stop during advancing
        clockAvanco = clockAvanco + clockExAdv; %Update total advancing time

        %- Converting video to images
        tic
        fprintf('Please wait, converting video to images... \n');
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        for m = 1: nImagesAdv %Getting a image from every video frame
            frame = read(vidFile,m);
            framegray = rgb2gray(frame);
            if iEns < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (iEns >= 10) && (iEns < 100)
                prefixo = "0";
            elseif (iEns >= 100) && (iEns < 1000)
                prefixo = "";
            end
            fullimName = strcat(prefixo,num2str(iEns),'_',fileName,imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            iEns = iEns + 1; %Interval of image name
        end
        delete(fullfilevidPath) %Delete video
        fprintf('Conversion of video to images completed! \n');
        clockImagesAdv(k) = toc;
        clockStopAdv(k) = clockStopAdv(k) + clockImagesAdv(k); %Update time spent each stop during advancing
        clockAvanco = clockAvanco + clockImagesAdv(k); %Update total advancing time
   
        %- Syringe displacement (syringe return to the drop)
        tic
        fprintf('Please wait, syringe displacement... \n');
        moverSeringa(passoMot,passoFuso,smSeringaZ,"down",dispZSer,"high"); %Syringe movement
        release(smSeringaZ); %Motor release
        fprintf('Syringe displacement completed! \n');
        clockRetSerAdv(k) = toc; %Time spent for syringe displacement (syringe return to the drop) during each stop
        clockStopAdv(k) = clockStopAdv(k) + clockRetSerAdv(k); %Update time spent each stop during advancing
        clockAvanco = clockAvanco + clockRetSerAdv(k); %Update total advancing time
        
        %- Syringe position monitoring (showing an example image)
        tic
        close all %Close all open windows
        frame = getsnapshot(vidobj); %Capture an image
        imshow(frame); %Show image
        movegui(gcf,'center') %Centralize image window
        set(gcf,'Name',strcat(fileName,'_CheckNeedlePosition'))
        %%{
        hold on;
        vidResValues = vidobj.VideoResolution; %Get resolution of the video object
        stepSizeY = vidResValues(2)/8; % Distance between lines in Y direction
        stepSizeX = vidResValues(1)/4; %Distance between lines in X direction
        for row = 1 : stepSizeY : vidResValues(2)
        yline(row, '-.r', 'LineWidth', 1);
        end
        for col = 1 : stepSizeX : vidResValues(1)
        xline(col, '-.r', 'LineWidth', 1);
        end
        %}
        clockCheckSerAdv = toc;
        clockStopAdv(k) = clockStopAdv(k) + clockCheckSerAdv; %Update time spent each stop during advancing
        clockAvanco = clockAvanco + clockCheckSerAdv; %Update total advancing time
        
        %- Deleting video object
        tic
        delete(vidobj)
        clear vidobj vidsrc
        clockExcVidAdv = toc;
        clockStopAdv(k) = clockStopAdv(k) + clockExcVidAdv; %Update time spent each stop during advancing
        clockAvanco = clockAvanco + clockExcVidAdv; %Update total advancing time
    end

    %2) RECEDING ROUTINE
    %- Begin of experiment
    tic
    fprintf('--------------------- EXECUTION OF DROP RECEDING ROUTINE ----------------------\n');
    %fprintf('Press any key to start the receding routine... \n');
    pause(2);
    close all %Close all open windows 
    %delete(vidobj) %Delete any video objects that may have been created
    clear vidobj vidsrc %Clear video related objects
    clockRecuo1 = toc;
    clockRecuo = clockRecuo + clockRecuo1; %Update total receding time

    %- Increase of drop volume
    tic
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, increasing drop volume...\n');
    move(smSeringa,floor(nPassosvolDepRecPrep)); %Drive motor connected to the syringe
    release(smSeringa); %Motor release
    fprintf('Increase of drop volume completed! \n');
    erroVolDepRecPrep = -passo2vol(nPassosvolDepRecPrep - floor(nPassosvolDepRecPrep),passoMot,passoSer); %Error in the drop preparation volume before receding (uL)
    erroVolTotalRec = erroVolTotalRec + erroVolDepRecPrep;
    volSerTotal = volSerTotal - volGotaVarRecPrep; %Update total liquid volume in the syringe
    pause(1); %Wait to start liquid withdrawal from the drop 
    clockRecGotaPrep = toc; %Time of preparaing drop volume before receding in s
    clockRecuo = clockRecuo + clockRecGotaPrep; %Update total receding time

    %- Initial withdrawal of drop volume
    tic
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, initial withdrawal of drop volume at high flow rate... \n');
    move(smSeringa,-floor(nPassosvolAdqRecIniAlta)); %Drive motor connected to the syringe
    release(smSeringa); %Motor release
    pause(1);
    fprintf('Please wait, initial withdrawal of drop volume at low flow rate... \n');
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    wMotorEns = uLs2RPM(taxaDepRecIni,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
    tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
    for i = 1:floor(nPassosvolAdqRecIniBaixa) %Initial withdrawal of drop volume
        move(smSeringa,-1); %Executes 1 rotation step on the motor
        pause(tadpasso);
    end
    release(smSeringa); %Motor release
    fprintf('Initial withdrawal of drop volume completed! \n');
    erroVolAdqRecIni = passo2vol((nPassosvolAdqRecIniAlta + nPassosvolAdqRecIniBaixa) - (floor(nPassosvolAdqRecIniAlta) + floor(nPassosvolAdqRecIniBaixa)),passoMot,passoSer); %Error in the initial drop volume before receding (uL)
    erroVolTotalRec = erroVolTotalRec + erroVolAdqRecIni;
    volSerTotal = volSerTotal + volGotaVarRecIni; %Update total liquid volume in the syringe
    fprintf('Please wait, waiting time for drop reach equilibrium... \n');
    pause(tEsperaRec); %Waiting time for drop reach equilibrium
    fprintf('Waiting time completed! \n');
    clockRecGotaIni = toc; %Time for initial drop receding in s before receding
    clockRecuo = clockRecuo + clockRecGotaIni; %Update total receding time

    %- Receding measurement series
    %iRec = 1; %Counter for all captured images during receding
    for k = 1:nStopsRec
        tic
        %- Create video object
        [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

        %- Update of image names
        fileName = strcat(name,'_Recuo_',num2str(k),'_V(uL)_',num2str((volGotaRecIni-k*volGotaVarRec),'%.2f'));
        fullvidName = strcat(fileName,vidDiskFmt);
        fullfilevidPath = fullfile(filePath,fullvidName);

        %- Preparation of video object
        diskLogger = VideoWriter(fullfilevidPath,vidComp);
        diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
        diskLogger.Quality = vidQual; %Definition of the video quality
        vidobj.DiskLogger = diskLogger;
        %vidobj.FramesPerTrigger = framerate*tCapt; %Total de frames a adquirir
        vidobj.FramesPerTrigger = nImagesRec; %Total frames to acquire
        vidobj.FrameGrabInterval = frameIntervalRec; %Frame interval

        %- Execution of measurement
        fprintf("------------------------------------------------------------------- \n");
        fprintf("Execution of receding measurement %d of %d \n",k,nStopsRec);
        clockPrepVidRec(k) = toc;
        clockStopRec(k) = clockStopRec(k) + clockPrepVidRec(k); %Update time spent each stop during receding
        clockRecuo = clockRecuo + clockPrepVidRec(k); %Update total receding time

        %- Additional decrease of drop volume
        tic
        fprintf('Please wait, additional decrease of drop volume... \n');
        if taxaDepRec <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
            smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
            wMotorEns = uLs2RPM(taxaDepRec,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
            tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
                for i = 1:floor(nPassosvolAdqRec) %Additional decrease of drop volume
                    move(smSeringa,-1); %Executes 1 rotation step on the motor
                    pause(tadpasso);
                end
        elseif (0.32 < taxaDepRec) && (taxaDepRec <= 1)
            wRPM = fuLs2RPM(taxaDepRec); %Calculation of the necessary rotational speed to be sent to the motor in RPM
            smSeringa.RPM = floor(wRPM);
            move(smSeringa,-floor(nPassosvolAdqRec)); %Motor movement
        end
        release(smSeringa); %Motor release
        fprintf('Additional decrease of drop volume completed! \n');
        erroVolAdqRec = passo2vol(nPassosvolAdqRec - floor(nPassosvolAdqRec),passoMot,passoSer); %Error in the volume acquired during receding (uL)
        erroVolTotalRec = erroVolTotalRec + erroVolAdqRec;
        volSerTotal = volSerTotal + volGotaVarRec; %Update total liquid volume in the syringe
        pause(1);
        clockRecGota(k) = toc;
        clockStopRec(k) = clockStopRec(k) + clockRecGota(k); %Update time spent each stop during receding
        clockRecuo = clockRecuo + clockRecGota(k); %Update total receding time

        %- Syringe displacement (syringe departure from the drop)
        tic
        fprintf('Please wait, syringe displacement... \n');
        moverSeringa(passoMot,passoFuso,smSeringaZ,"up",dispZSer,"high"); %Syringe movement
        release(smSeringaZ); %Motor release
        fprintf('Syringe displacement completed! \n');
        clockDepSerRec(k) = toc;
        clockStopRec(k) = clockStopRec(k) + clockDepSerRec(k); %Update time spent each stop during receding
        clockRecuo = clockRecuo + clockDepSerRec(k); %Update total receding time

        %- Waiting time for drop reach equilibrium
        tic
        fprintf('Please wait, waiting time for drop reach equilibrium...\n');
        if clockDispSerRec1(k) >= tEsperaRec %Syringe displacement time is longer than waiting time for drop equilibrium
            tEsperaRecReal(k) = clockDispSerRec1(k); %Calculation of actual waiting time
        else %Syringe displacement time is less than waiting time for drop equilibrium
            tEsperaRecReal(k) = tEsperaRec;
            deltat = tEsperaRec - clockDispSerRec1(k);
            pause(deltat)
        end
        %pause(tEsperaRec);
        fprintf('Waiting time completed! \n');
        clockEspRec(k) = toc;
        clockStopRec(k) = clockStopRec(k) + clockEspRec(k); %Update time spent each stop during receding
        clockRecuo = clockRecuo + clockEspRec(k); %Update total receding time

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
        clockVidRec(k) = toc;
        clockStopRec(k) = clockStopRec(k) + clockVidRec(k); %Update time spent each stop during receding
        clockRecuo = clockRecuo + clockVidRec(k); %Update total receding time
        
        %- Showing an example image
        tic
        close all %Close all open windows
        frame = getsnapshot(vidobj); %Capture an image
        imshow(frame); %Show image
        movegui(gcf,'center') %Centralize image window
        set(gcf,'Name',strcat(fileName,'_Example'))
        clockExRec = toc;
        clockStopRec(k) = clockStopRec(k) + clockExRec; %Update time spent each stop during receding
        clockRecuo = clockRecuo + clockExRec; %Update total receding time

        %- Converting video to images
        tic
        fprintf('Please wait, converting video to images... \n');
        %vid2frame(filePath,fullvidName,fileName,imFmt,tCaptInterval);
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        %i = 1;
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
            fullimName = strcat(prefixo,num2str(iEns),'_',fileName,imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            iEns = iEns + 1; %Interval of image name
        end
        delete(fullfilevidPath) %Deletar vídeo
        fprintf('Conversion of video to images completed! \n');
        clockImagesRec(k) = toc;
        clockStopRec(k) = clockStopRec(k) + clockImagesRec(k); %Update time spent each stop during receding
        clockRecuo = clockRecuo + clockImagesRec(k); %Update total receding time
        
        %- Syringe displacement (syringe return to the drop)
        tic
        fprintf('Please wait, syringe displacement... \n');
        moverSeringa(passoMot,passoFuso,smSeringaZ,"down",dispZSer,"high"); %Syringe movement
        release(smSeringaZ); %Motor release
        fprintf('Syringe displacement completed! \n');
        clockRetSerRec(k) = toc;
        clockStopRec(k) = clockStopRec(k) + clockRetSerRec(k); %Update time spent each stop during receding
        clockRecuo = clockRecuo + clockRetSerRec(k); %Update total receding time
        
        %- Syringe position monitoring (showing an example image)
        tic
        close all %Close all open windows
        frame = getsnapshot(vidobj); %Capture an image
        imshow(frame); %Show image
        movegui(gcf,'center') %Centralize image window
        set(gcf,'Name',strcat(fileName,'_CheckNeedlePosition'))
        %%{
        hold on;
        vidResValues = vidobj.VideoResolution; %Get resolution of the video object
        stepSizeY = vidResValues(2)/8; % Distance between lines in Y direction
        stepSizeX = vidResValues(1)/4; %Distance between lines in X direction
        for row = 1 : stepSizeY : vidResValues(2)
        yline(row, '-.r', 'LineWidth', 1);
        end
        for col = 1 : stepSizeX : vidResValues(1)
        xline(col, '-.r', 'LineWidth', 1);
        end
        %}
        clockCheckSerRec = toc;
        clockStopRec(k) = clockStopRec(k) + clockCheckSerRec; %Update time spent each stop during receding
        clockRecuo = clockRecuo + clockCheckSerRec; %Update total receding time

        %- Deleting video object
        tic
        delete(vidobj)
        clear vidobj vidsrc
        clockExcVidRec = toc;
        clockStopRec(k) = clockStopRec(k) + clockExcVidRec; %Update time spent each stop during receding
        clockRecuo = clockRecuo + clockExcVidRec; %Update total receding time
    end
    clockEnsaio = clockAvanco + clockRecuo; %Total time spent in the experiment including advancing and receding 
    erroVolTotal = erroVolTotalAdv + erroVolTotalRec; %Total error in the final drop volume including advancing and receding

    %- Export/save experiment description
    filenametxt = "Description.txt";
    fullfilenametxt = fullfile(filePath,filenametxt);
    varnames = {'Datetime','volGotaStart[uL]','volGotaAdvIni[uL]','taxaDepAdvIni[uL]','volGotaAdvMax[uL]','taxaDepAdv[uL/s]','deslZSer[mm]','tEsperaAdv[s]','tEsperaAdvReal[s]','tCaptAdv[s]','nImagesAdv','nStopsAdv','volGotaRecPrep[uL]','volGotaRecIni[uL]','volGotaRecEnd[uL]','taxaDepRecIni[uL/s]','taxaDepRec[uL/s]','tEsperaRec[s]','tEsperaRecReal[s]','tCaptRec[s]','nImagesRec','nStopsRec','volGotaVarAdv[uL]','volGotaVarRec[uL]','totalFrameAdv','totalFrameRec','frameIntervalAdv','frameIntervalRec','erroVolDepAdvIni[uL]','erroVolDepAdv[uL]','erroVolTotalAdv[uL]','erroVolDepRecPrep[uL]','erroVolAdqRecIni[uL]','erroVolAdqRec[uL]','erroVolTotalRec[uL]','erroVolTotal[uL]','clockAdvGotaIni[s]','clockEspAdvGotaIni[s]','clockAvgPrepVidAdv[s]','clockTotPrepVidAdv[s]','clockAvgAdvGota[s]','clockTotAdvGota[s]','clockAvgDepSerAdv[s]','clockTotDepSerAdv[s]','clockAvgEspAdv[s]','clockTotEspAdv[s]','clockAvgVidAdv[s]','clockTotVidAdv[s]','clockAvgImagesAdv[s]','clockTotImagesAdv[s]','clockAvgRetSerAdv[s]','clockTotRetSerAdv[s]','clockAvgStopAdv[s]','clockTotStopAdv[s]','clockTotalAvanco[s]','clockRecGotaPrep[s]','clockRecGotaIni[s]','clockAvgPrepVidRec[s]','clockTotPrepVidRec[s]','clockAvgRecGota[s]','clockTotRecGota[s]','clockAvgDepSerRec[s]','clockTotDepSerRec[s]','clockAvgEspRec[s]','clockTotEspRec[s]','clockAvgVidRec[s]','clockTotVidRec[s]','clockAvgImagesRec[s]','clockTotImagesRec[s]','clockAvgRetSerRec[s]','clockTotRetSerRec[s]','clockAvgStopRec[s]','clockTotStopRec[s]','clockTotalRecuo[s]','clockTotalEnsaio[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
    T = table(datetime,volGotaStart,volGotaAdvIni,taxaDepAdvIni,volGotaAdvMax,taxaDepAdv,dispZSer,tEsperaAdv,round(mean(tEsperaAdvReal),2),tCaptAdv,nImagesAdv,nStopsAdv,volGotaRecPrep,volGotaRecIni,volGotaRecEnd,taxaDepRecIni,taxaDepRec,tEsperaRec,round(mean(tEsperaRecReal),2),tCaptRec,nImagesRec,nStopsRec,volGotaVarAdv,volGotaVarRec,totalFrameAdv,totalFrameRec,frameIntervalAdv,frameIntervalRec,round(erroVolDepAdvIni,4),round(erroVolDepAdv,4),round(erroVolTotalAdv,4),round(erroVolDepRecPrep,4),round(erroVolAdqRecIni,4),round(erroVolAdqRec,4),round(erroVolTotalRec,4),round(erroVolTotal,4),round(clockAdvGotaIni,2),round(clockEspAdvGotaIni,2),round(mean(clockPrepVidAdv),2),round(sum(clockPrepVidAdv),2),round(mean(clockAdvGota),2),round(sum(clockAdvGota),2),round(mean(clockDepSerAdv),2),round(sum(clockDepSerAdv),2),round(mean(clockEspAdv),2),round(sum(clockEspAdv),2),round(mean(clockVidAdv),2),round(sum(clockVidAdv),2),round(mean(clockImagesAdv),2),round(sum(clockImagesAdv),2),round(mean(clockRetSerAdv),2),round(sum(clockRetSerAdv),2),round(mean(clockStopAdv),2),round(sum(clockStopAdv),2),round(clockAvanco,2),round(clockRecGotaPrep,2),round(clockRecGotaIni,2),round(mean(clockPrepVidRec),2),round(sum(clockPrepVidRec),2),round(mean(clockRecGota),2),round(sum(clockRecGota),2),round(mean(clockDepSerRec),2),round(sum(clockDepSerRec),2),round(mean(clockEspRec),2),round(sum(clockEspRec),2),round(mean(clockVidRec),2),round(sum(clockVidRec),2),round(mean(clockImagesRec),2),round(sum(clockImagesRec),2),round(mean(clockRetSerRec),2),round(sum(clockRetSerRec),2),round(mean(clockStopRec),2),round(sum(clockStopRec),2),round(clockRecuo,2),round(clockEnsaio,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
    writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file
    pause(5);
end
end

