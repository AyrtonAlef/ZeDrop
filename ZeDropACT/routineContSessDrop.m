function [volSerTotal] = routineContSessDrop(vidAdap,vidID,vidFmt,imFmt,smSeringa,passoMot,passoSer,wMin,framerate,limTotalFrames,volSerTotal)
%ROUTINECONTSESSDROP Execution of continuous sessile drop experiment 
% (needle in). Includes advancing and receding routines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execução de ensaio contí­nuo (needle in) de gota séssil. Inclui rotinas
%   de avanço e recuo da gota. Imagens são capturadas concomitantemente com
%   o aumento/redução contínuo(a) do volume da gota. Inicialmente, 
%   realiza-se um aumento do volume da gota até que se alcance o volume 
%   indicado (volGotaAdvIni) para a realização da rotina de avanço. Após a 
%   conclusão da rotina de avanço, realiza-se um novo aumento do volume da 
%   gota e subsequente redução do volume da gota, em alta e baixa 
%   velocidade, até que se alcance o volume de lí­quido para a realização 
%   da rotina de recuo (volGotaRecIni). A rotina requer que o usuário 
%   defina o volume inicial da gota (volGotaStar), o volume máximo da gota 
%   após o avanço (volGotaAdvMax), a taxa de deposição da gota durante o 
%   avanço (taxaDepAdv), tempo de espera para equilíbrio da gota durante o 
%   avanço (tEsperaAdv), número de imagens capturadas durante o avanço 
%   (nImagesAdv), volume de preparação da gota para o recuo 
%   (volGotaRecPrep), volume inicial da gota antes do recuo 
%   (volGotaRecIni), taxa de deposição da gota durante o recuo (taxaDepRec)
%   , tempo de espera para equilíbrio da gota durante o recuo (tEsperaRec) 
%   e o número de imagens capturadas durante o recuo (nImagesRec).  A 
%   rotina permite a escolha do  diretório em que as imagens serão salvas, 
%   assim como a definição do nome das imagens. Ao final, um arquivo de 
%   texto é exportado com todas os parâmetros relevantes sobre o ensaio. 
%   Antes da execução do ensaio, requer-se que a agulha esteja inserida na 
%   gota séssil.
%   Execution of continuous sessile drop experiment (needle in). Includes
%   advancing and receding routines. Images are captured concurrently with 
%   a continuous increase/decrease of drop volume. Initially, the drop 
%   volume is increased until the indicated volume is reached 
%   (volGotaAdvIni) to carry out the advancing routine. After completing 
%   the advancing routine, a new increase in the drop volume is carried out
%   and a subsequent reduction of the drop volume, at high and low flow 
%   rates, until the liquid volume is reached to carry out the receding 
%   routine (volGotaRecIni). The routine requires the user to define the 
%   initial drop volume (volGotaStart), the maximum drop volume after 
%   advancing (volGotaAdvMax), the flow rate during advancing (taxaDepAdv),
%   waiting time for drop equilibrium during advancing (tEsperaAdv), number
%   of images captured during advancing (nImagesAdv), preparation drop 
%   volume for receding (volGotaRecPrep), initial drop volume before 
%   receding (volGotaRecIni), flow rate during receding (taxaDepRec),
%   waiting time for drop equilibrium during receding (tEsperaRec) and the 
%   number of images captured during receding (nImagesRec). The routine 
%   allows the user to choose the directory where the images will be 
%   exported, as well as to define the name of the images. At the end, a 
%   text file is exported with all relevant parameters of the experiment. 
%   Prior to start the experiment, the needle is required to be inserted 
%   into the sessile drop.  
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

fprintf('--------- SESSILE DROP (NEEDLE IN) - CONTINUOUS EXPERIMENT (ADVANCING AND RECEDING) -----------\n');
% Experiment parameters (entered by the user)
fprintf('1) Experiment parameters during advancing: \n');
volGotaStart = input("- Initial drop volume (uL): "); %Initial drop volume
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
    taxaDepAdv = input("- - Flow rate during advancing (uL/s) (max. 1 uL/s): "); %Drop flow rate in uL/s during advancing
    if taxaDepAdv > 1 %Check that the flow rate does not exceed the maximum rotation speed 
        fprintf("Unable to proceed with the routine. Please lower the flow rate. \n")
    else
        break;
    end
end
tEsperaAdv = input("- Waiting time for drop equilibrium during advancing (s): "); %Waiting time for the drop to reach equilibrium during advancing protocol
nImagesAdv = input("- Number of images captured during advancing: "); %Number of images captured during advancing

fprintf('2) Experimental parameters during receding: \n');
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
nImagesRec = input("- Number of captured images during receding: "); %Number of captured images during receding

% Experiment parameters (calculated)
volGotaVarAdvIni = volGotaAdvIni - volGotaStart; %Volume variation during intial drop increase before advancing
volGotaVarAdv = volGotaAdvMax - volGotaAdvIni; %Volume variation during advancing
tEnsaioAdv = volGotaVarAdv/taxaDepAdv;  %Total experiment time (drop volume increase) during advancing in s
totalFrameAdv = ceil(framerate*tEnsaioAdv); %Number of total frames acquired during advancing
frameIntervalAdv = round(totalFrameAdv/nImagesAdv); %Frame interval during advancing
tCaptIntervalAdv = frameIntervalAdv/framerate; %Time interval between captured frames during advancing in s 
volCaptIntervalAdv = taxaDepAdv*tCaptIntervalAdv; %Volume interval between captured frames during advancing in uL 

volGotaVarRecPrep = volGotaRecPrep - volGotaAdvMax;  %Volume variation during preparation (drop increase) before receding
volGotaVarRecIni = abs(volGotaRecPrep - volGotaRecIni); %Volume variation during the initial stage of receding (drop decrease)
volGotaVarRecIniAlta = floor(volGotaVarRecIni*0.8); %Volume variation during the initial receding stage (drop decrease) at high flow rate
volGotaVarRecIniBaixa = volGotaVarRecIni - volGotaVarRecIniAlta; %Volume variation during the initial receding stage (drop decrease) at low flow rate
volGotaVarRec = abs(volGotaRecIni - volGotaRecEnd); %Volume change during receding
tEnsaioRec = volGotaVarRec/taxaDepRec; %Total experiment time (drop volume decrease) during receding in s
totalFrameRec = ceil(framerate*tEnsaioRec); %Total number of frames acquired during receding
frameIntervalRec = round(totalFrameRec/nImagesRec); %Frame interval during receding
tCaptIntervalRec = frameIntervalRec/framerate; %Time interval between captured frames during receding in s
volCaptIntervalRec = taxaDepRec*tCaptIntervalRec; %Volume variation between captured frames during receding in uL

erroVolTotalAdv = 0;
erroVolTotalRec = 0;
clockAvanco = 0;
clockRecuo = 0;
clockEnsaio = 0;
%volCaptInterval = taxaDepGota*tCaptInterval; %Volume change between images in uL

if volGotaRecPrep >= volSerTotal  %Check that there is enough volume in the syringe to perform the procedure
    fprintf("Unable to proceed with the routine. Please increase the liquid volume in the syringe. \n")
else
    %- Configuração de movimentação do motor
    nPassosvolDepAdvIni = vol2passo(volGotaVarAdvIni,passoMot,passoSer); %No of steps corresponding to the deposited volume during preparation for advancing
    nPassosvolDepAdv = vol2passo(volGotaVarAdv,passoMot,passoSer); %No of steps corresponding to the deposited volume between each measurement during advancing
    nPassosvolDepRecPrep = vol2passo(volGotaVarRecPrep,passoMot,passoSer); %No of steps corresponding to the volume released during preparation for receding
    nPassosvolAdqRecIniAlta = vol2passo(volGotaVarRecIniAlta,passoMot,passoSer); %No of steps corresponding to the volume withdrawal before the start of receding in high flow rate
    nPassosvolAdqRecIniBaixa = vol2passo(volGotaVarRecIniBaixa,passoMot,passoSer); %No of steps corresponding to the volume withdrawal before the start of receding in low flow rate
    nPassosvolAdqRec = vol2passo(volGotaVarRec,passoMot,passoSer); %No of steps corresponding to the volume withdrawal during receding
    
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
    fprintf('-------------------- EXECUTION OF DROP ADVANCING ROUTINE ---------------------\n');
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
    fprintf('Please wait, inital drop increase... \n');
    wMotorEns = uLs2RPM(taxaDepAdvIni,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
    tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
    for i = 1:floor(nPassosvolDepAdvIni) %Initial drop increase
        move(smSeringa,1); %Executes 1 rotation step on the motor
        pause(tadpasso);
    end
    release(smSeringa); %Motor release
    fprintf('Initial drop increase completed! \n');
    erroVolDepAdvIni = -passo2vol(nPassosvolDepAdvIni - floor(nPassosvolDepAdvIni),passoMot,passoSer); %Error in the initial drop volume formed before advancing (uL)
    erroVolTotalAdv = erroVolTotalAdv + erroVolDepAdvIni;
    volSerTotal = volSerTotal - volGotaVarAdvIni; %Update total liquid volume in the syringe
    clockAdvGotaIni = toc; %Time for initial drop increase in s before advancing
    clockAvanco = clockAvanco + clockAdvGotaIni; %Update total advancing time

    %- Waiting time for drop reach equilibrium
    tic
    fprintf('Please wait, waiting time for drop reach equilibrium...\n');
    pause(tEsperaAdv); %Waiting time for drop reach equilibrium
    fprintf('Waiting time completed!\n');
    clockEspAdvGotaIni = toc; %Time for waiting for drop equilibrium in s
    clockAvanco = clockAvanco + clockEspAdvGotaIni; %Update total advancing time
    
    %- Continuous increase of drop volume
    tic
    %-- Create video object
    [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

    %- Update of image names
    %fileName = strcat(name,'_Avanco_',num2str(k),'_V(uL)_',num2str((volGotaAdvIni+k*volGotaVarAdv),'%.2f'));
    %fileName = name;
    fileName = strcat(name,"_Avanco");
    fullvidName = strcat(fileName,vidDiskFmt);
    fullfilevidPath = fullfile(filePath,fullvidName);

    %- Preparation of video object
    diskLogger = VideoWriter(fullfilevidPath,vidComp);
    diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
    diskLogger.Quality = vidQual; %Definition of the video quality
    vidobj.DiskLogger = diskLogger;
    if totalFrameAdv <= limTotalFrames %Continuous frame capture during video recording
        vidobj.FramesPerTrigger = totalFrameAdv; %Total frames to acquire
        vidobj.FrameGrabInterval = 1; %Frame interval
    else %Capturing spaced frames during video recording
        vidobj.FramesPerTrigger = nImagesAdv; %Total frames to acquire
        vidobj.FrameGrabInterval = frameIntervalAdv; %Frame interval
    end
    clockPrepVidAdv = toc;
    clockAvanco = clockAvanco + clockPrepVidAdv; %Update total advancing time

    %-- Start video recording
    tic
    fprintf('Start video recording. \n');
    start(vidobj) %Start video recording
    %-- Continuous increase in drop volume
    fprintf('Please wait, continuous increase in drop volume... \n');
    if taxaDepAdv <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
        smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
        wMotorEns = uLs2RPM(taxaDepAdv,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
        tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed   
        for i = 1:floor(nPassosvolDepAdv) %Continuous drop volume release
            move(smSeringa,1); %Executes 1 rotation step on the motor
            pause(tadpasso);
        end
    elseif (0.32 < taxaDepAdv) && (taxaDepAdv <= 1)
        wRPM = fuLs2RPM(taxaDepAdv); %Calculation of the necessary rotational speed to be sent to the motor in RPM
        taxaDepAdv = fRPM2uLs(ceil(wRPM)); %Update flow rate
        % Update variables
        tEnsaioAdv = volGotaVarAdv/taxaDepAdv;
        totalFrameAdv = ceil(framerate*tEnsaioAdv);
        frameIntervalAdv = round(totalFrameAdv/nImagesAdv);
        tCaptIntervalAdv = frameIntervalAdv/framerate;
        volCaptIntervalAdv = taxaDepAdv*tCaptIntervalAdv;
        % Syringe plunger movement
        smSeringa.RPM = round(wRPM);
        move(smSeringa,floor(nPassosvolDepAdv)); %Motor movement
    end
    release(smSeringa); %Motor release
    fprintf('Continuous increase in drop volume completed! \n');
    erroVolDepAdv = passo2vol(nPassosvolDepAdv - floor(nPassosvolDepAdv),passoMot,passoSer); %Error in the volume released during the experiment
    erroVolTotalAdv = erroVolTotalAdv + erroVolDepAdv; %Update error on total drop volume during advancing
    volSerTotal = volSerTotal - volGotaVarAdv; %Update total liquid volume in the syringe
    %-- Finishing video recording
    wait(vidobj) %Wait until all frames are captured
    stop(vidobj);
    %--- Wait for all frames to be written to disk
    while vidobj.FramesAcquired ~= vidobj.DiskLoggerFrameCount
        pause(.1);
    end
    fprintf('Video recording completed! \n');
    clockVidAdv = toc;
    clockAvanco = clockAvanco + clockVidAdv; %Update total advancing time
      
    %-- Showing an example image
    tic
    close all %Close all open windows
    frame = getsnapshot(vidobj); %Capture an image
    imshow(frame); %Show image
    movegui(gcf,'center') %Centralize image window
    clockExAdv = toc;
    clockAvanco = clockAvanco + clockExAdv; %Update total advancing time
   
    %- Converting video to images
    tic
    fprintf('Please wait, converting video to images... \n');
    iEns = 1;
    if totalFrameAdv <= limTotalFrames %Continuous frame capture during video recording
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        %iEns = 1;
        iAdv = 1;
        for k = 1: round(frameIntervalAdv): totalFrameAdv %Getting a image from every video second
            frame = read(vidFile,k);
            framegray = rgb2gray(frame);
            if iEns < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (iEns >= 10) && (iEns < 100)
                prefixo = "0";
            elseif (iEns >= 100) && (iEns < 1000)
                prefixo = "";
            end
            %fullimName = strcat(prefixo,num2str(i),'_',fileName,'_V(uL)_',num2str((volGotaIni + k*volCapIntervalframe),'%.2f'),imFmt);
            fullimName = strcat(prefixo,num2str(iEns),'_',fileName,'_V(uL)_',num2str((volGotaAdvIni + iAdv*volCaptIntervalAdv),'%.2f'),imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            iAdv = iAdv + 1;
            iEns = iEns + 1; %Interval of image name
        end
    else %Capturing spaced frames during video recording
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        for k = 1: nImagesAdv %Getting a image from every video frame
            frame = read(vidFile,k);
            framegray = rgb2gray(frame);
            if iEns < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (iEns >= 10) && (iEns < 100)
                prefixo = "0";
            elseif (iEns >= 100) && (iEns < 1000)
                prefixo = "";
            end
            %fullimName = strcat(fileName,'_',num2str(i),'_V(uL)_',num2str((volGotaIni + k*volCaptInterval),'%.2f'),imFmt);
            fullimName = strcat(prefixo,num2str(iEns),'_',fileName,'_V(uL)_',num2str((volGotaAdvIni + k*volCaptIntervalAdv),'%.2f'),imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            iEns = iEns + 1; %Interval of image name
        end
    end
    delete(fullfilevidPath) %Delete video
    fprintf('Conversion of video to images completed! \n');
    clockImagesAdv = toc;
    clockAvanco = clockAvanco + clockImagesAdv; %Update total advancing time

    %- Deleting video object
    tic
    delete(vidobj)
    clear vidobj vidsrc
    clockExcVidAdv = toc;
    clockAvanco = clockAvanco + clockExcVidAdv; %Update total advancing time

    %2) RECEDING ROUTINE
    %- Begin experiment
    tic
    fprintf('-------------------- EXECUTION OF DROP RECEDING ROUTINE ---------------------\n');
    %fprintf('Press any key to proceed with the receding routine... \n');
    %pause;
    close all %Close all open windows
    %delete(vidobj) %Delete any video objects that may have been created
    clear vidobj vidsrc %Clear video related objects
    clockRecuo1= toc;
    clockRecuo = clockRecuo + clockRecuo1; %Update total receding time

    %- Increase of drop volume
    tic
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, increasing drop volume...\n');
    move(smSeringa,floor(nPassosvolDepRecPrep)); %Drive motor connected to syringe
    release(smSeringa); %Motor release
    fprintf('Increase of drop volume completed!\n');
    erroVolDepRecPrep = -passo2vol(nPassosvolDepRecPrep - floor(nPassosvolDepRecPrep),passoMot,passoSer); %Error in the preparation drop volume before receding (uL)
    erroVolTotalRec = erroVolTotalRec + erroVolDepRecPrep;
    volSerTotal = volSerTotal - volGotaVarRecPrep; %Update total liquid volume in the syringe
    pause(2); %Wait to start liquid withdrawal from the drop 
    clockRecGotaPrep = toc; %Time of preparaing drop volume before receding in s
    clockRecuo = clockRecuo + clockRecGotaPrep; %Update total receding time

    %- Initial withdrawal of drop volume
    tic
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, withdrawal of drop volume at high flow rate... \n');
    move(smSeringa,-floor(nPassosvolAdqRecIniAlta)); %Drive motor connected to the syringe
    release(smSeringa); %Motor release
    pause(1);
    fprintf('Please wait, withdrawal of drop volume at low flow rate...\n');
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    wMotorEns = uLs2RPM(taxaDepRecIni,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
    tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
    for i = 1:floor(nPassosvolAdqRecIniBaixa) %Initial withdrawal of drop volume
        move(smSeringa,-1); %Executes 1 rotation step on the motor
        pause(tadpasso);
    end
    release(smSeringa); %Motor release
    fprintf('Withdrawal of drop volume completed! \n');
    erroVolAdqRecIni = passo2vol((nPassosvolAdqRecIniAlta + nPassosvolAdqRecIniBaixa) - (floor(nPassosvolAdqRecIniAlta) + floor(nPassosvolAdqRecIniBaixa)),passoMot,passoSer); %Error in the initial drop volume before receding (uL)
    erroVolTotalRec = erroVolTotalRec + erroVolAdqRecIni;
    volSerTotal = volSerTotal + volGotaVarRecIni; %Update total liquid volume in the syringe
    clockRecGotaIni = toc; %Waiting time for drop reach equilibrium
    clockRecuo = clockRecuo + clockRecGotaIni; %Update total receding time

    %- %Waiting time for drop reach equilibrium
    tic
    fprintf('Please wait, waiting time for drop reach equilibrium... \n');
    pause(tEsperaRec); %Waiting time for drop reach equilibrium
    fprintf('Waiting time completed! \n');
    clockEspRecGotaIni = toc; %Time for initial drop receding in s before receding 
    clockRecuo = clockRecuo + clockEspRecGotaIni; %Update total receding time

    %- Continuous decrease of drop volume
    tic
    %-- Create video object
    [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

    %- Update of image names
    %fileName = strcat(name,'_Avanco_',num2str(k),'_V(uL)_',num2str((volGotaAdvIni+k*volGotaVarAdv),'%.2f'));
    %fileName = name;
    fileName = strcat(name,"_Recuo");
    fullvidName = strcat(fileName,vidDiskFmt);
    fullfilevidPath = fullfile(filePath,fullvidName);

    %- Preparation of video object
    diskLogger = VideoWriter(fullfilevidPath,vidComp);
    diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
    diskLogger.Quality = vidQual; %Definition of the video quality
    vidobj.DiskLogger = diskLogger;
    if totalFrameRec <= limTotalFrames%Continuous frame capture during video recording
        vidobj.FramesPerTrigger = totalFrameRec; %Total frames to acquire
        vidobj.FrameGrabInterval = 1; %Frame interval
    else %Capturing spaced frames during video recording
        vidobj.FramesPerTrigger = nImagesRec; %Total frames to acquire
        vidobj.FrameGrabInterval = frameIntervalRec; %Frame interval
    end
    clockPrepVidRec = toc;
    clockRecuo = clockRecuo + clockPrepVidRec; %Update total receding time

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
        for i = 1:floor(nPassosvolAdqRec) %Continuous drop volume withdrawal
            move(smSeringa,-1); %Executes 1 rotation step on the motor
            pause(tadpasso);
        end
    elseif (0.32 < taxaDepRec) && (taxaDepRec <= 1)
        wRPM = fuLs2RPM(taxaDepRec); %Calculation of the necessary rotational speed to be sent to the motor in RPM
        taxaDepRec = fRPM2uLs(ceil(wRPM)); %Update flow rate
        % Update variables
        tEnsaioRec = volGotaVarRec/taxaDepRec;
        totalFrameRec = ceil(framerate*tEnsaioRec);
        frameIntervalRec = round(totalFrameRec/nImagesRec);
        tCaptIntervalRec = frameIntervalRec/framerate;
        volCaptIntervalRec = taxaDepRec*tCaptIntervalRec;
        % Syringe plunger movement
        smSeringa.RPM = round(wRPM);
        move(smSeringa,-floor(nPassosvolAdqRec)); %Motor movement
    end
    release(smSeringa); %Motor release
    fprintf('Continuous decrease in drop volume completed! \n');
    erroVolAdqRec = passo2vol(nPassosvolAdqRec - floor(nPassosvolAdqRec),passoMot,passoSer); %Error in the volume withdrawal during the experiment
    erroVolTotalRec = erroVolTotalRec + erroVolAdqRec; %Update error on total drop volume
    volSerTotal = volSerTotal + volGotaVarRec; %Update total liquid volume in the syringe
    %-- Finishing video recording
    wait(vidobj) %Wait until all frames are captured
    stop(vidobj);
    %--- Wait for all frames to be written to disk
    while vidobj.FramesAcquired ~= vidobj.DiskLoggerFrameCount
        pause(.1);
    end
    fprintf('Video recording completed! \n');
    clockVidRec = toc;
    clockRecuo = clockRecuo + clockVidRec; %Update total receding time
      
    %-- Showing an example image
    tic
    close all %Close all open windows
    frame = getsnapshot(vidobj); %Capture an image
    imshow(frame); %Show image
    movegui(gcf,'center') %Centralize image window
    clockExRec = toc;
    clockRecuo = clockRecuo + clockExRec; %Update total receding time
   
    %- Converting video to images
    tic
    fprintf('Please wait, converting video to images... \n');
    if totalFrameRec <= limTotalFrames %Continuous frame capture during video recording
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        iRec = 1;
        for k = 1: round(frameIntervalRec): totalFrameRec %Getting a image from every video second
            frame = read(vidFile,k);
            framegray = rgb2gray(frame);
            if iEns < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (iEns >= 10) && (iEns < 100)
                prefixo = "0";
            elseif (iEns >= 100) && (iEns < 1000)
                prefixo = "";
            end
            %fullimName = strcat(prefixo,num2str(i),'_',fileName,'_V(uL)_',num2str((volGotaIni + k*volCapIntervalframe),'%.2f'),imFmt);
            fullimName = strcat(prefixo,num2str(iEns),'_',fileName,'_V(uL)_',num2str((volGotaRecIni - iRec*volCaptIntervalRec),'%.2f'),imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            iRec = iRec + 1; 
            iEns = iEns + 1; %Interval of image name
        end
    else %Capturing spaced frames during video recording
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        for k = 1: nImagesRec %Getting a image from every video frame
            frame = read(vidFile,k);
            framegray = rgb2gray(frame);
            if iEns < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (iEns >= 10) && (iEns < 100)
                prefixo = "0";
            elseif (iEns >= 100) && (iEns < 1000)
                prefixo = "";
            end
            %fullimName = strcat(fileName,'_',num2str(i),'_V(uL)_',num2str((volGotaIni + k*volCaptInterval),'%.2f'),imFmt);
            fullimName = strcat(prefixo,num2str(iEns),'_',fileName,'_V(uL)_',num2str((volGotaRecIni - k*volCaptIntervalRec),'%.2f'),imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            iEns = iEns + 1; %Interval of image name
        end
    end
    delete(fullfilevidPath) %Delete video
    fprintf('Conversion of video to images completed! \n');
    clockImagesRec = toc;
    clockRecuo = clockRecuo + clockImagesRec; %Update total receding time

    %- Deleting video object
    tic
    delete(vidobj)
    clear vidobj vidsrc
    clockExcVidRec = toc;
    clockRecuo = clockRecuo + clockExcVidRec; %Update total receding time

    % Calculation of total time and error
    clockEnsaio = clockAvanco + clockRecuo; %Total time spent on the entire routine including advancing and receding
    erroVolTotal = erroVolTotalAdv + erroVolTotalRec; %Total error in final drop volume including advancing and receding

    %- Export/save experiment description
    filenametxt = "Description.txt";
    fullfilenametxt = fullfile(filePath,filenametxt);
    varnames = {'Datetime','volGotaStart[uL]','volGotaAdvIni[uL]','taxaDepAdvIni[uL]','volGotaAdvMax[uL]','taxaDepAdv[uL/s]','tEsperaAdv[s]','tEnsaioAdv[s]','nImagesAdv','totalFrameAdv','frameIntervalAdv','tCaptIntervalAdv[s]','volCaptIntervalAdv[uL]','volGotaRecPrep[uL]','volGotaRecIni[uL]','volGotaRecEnd[uL]','taxaDepRecIni[uL/s]','taxaDepRec[uL/s]','tEsperaRec[s]','tEnsaioRec[s]','nImagesRec','totalFrameRec','frameIntervalRec','tCaptIntervalRec[s]','volCaptIntervalRec[uL]','volGotaVarAdv[uL]','volGotaVarRec[uL]','erroVolDepAdvIni[uL]','erroVolDepAdv[uL]','erroVolTotalAdv[uL]','erroVolDepRecPrep[uL]','erroVolAdqRecIni[uL]','erroVolAdqRec[uL]','erroVolTotalRec[uL]','erroVolTotal[uL]','clockAdvGotaIni[s]','clockEspAdvGotaIni[s]','clockPrepVidAdv[s]','clockVidAdv[s]','clockImagesAdv[s]','clockTotalAvanco[s]','clockRecGotaPrep[s]','clockRecGotaIni[s]','clockEspRecGotaIni[s]','clockPrepVidRec[s]','clockVidRec[s]','clockImagesRec[s]','clockTotalRecuo[s]','clockTotalEnsaio[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
    T = table(datetime,volGotaStart,volGotaAdvIni,taxaDepAdvIni,volGotaAdvMax,taxaDepAdv,tEsperaAdv,tEnsaioAdv,nImagesAdv,totalFrameAdv,frameIntervalAdv,tCaptIntervalAdv,volCaptIntervalAdv,volGotaRecPrep,volGotaRecIni,volGotaRecEnd,taxaDepRecIni,taxaDepRec,tEsperaRec,tEnsaioRec,nImagesRec,totalFrameRec,frameIntervalRec,tCaptIntervalRec,volCaptIntervalRec,volGotaVarAdv,volGotaVarRec,round(erroVolDepAdvIni,4),round(erroVolDepAdv,4),round(erroVolTotalAdv,4),round(erroVolDepRecPrep,4),round(erroVolAdqRecIni,4),round(erroVolAdqRec,4),round(erroVolTotalRec,4),round(erroVolTotal,4),round(clockAdvGotaIni,2),round(clockEspAdvGotaIni,2),round(clockPrepVidAdv,2),round(clockVidAdv,2),round(clockImagesAdv,2),round(clockAvanco,2),round(clockRecGotaPrep,2),round(clockRecGotaIni,2),round(clockEspRecGotaIni,2),round(clockPrepVidRec,2),round(clockVidRec,2),round(clockImagesRec,2),round(clockRecuo,2),round(clockEnsaio,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
    writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file
    pause(5);
end
end

