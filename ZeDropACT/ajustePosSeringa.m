function ajustePosSeringa(smSeringaZ,passoMot,passoFuso)
%AJUSTEPOSSERINGA Adjustment of the vertical position of the syringe
%   INPUT
% smSerignaZ - Syringe displacement motor control
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% passoSeringa - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Allows the adjustment of the vertical positioning of the syringe. It
%   allows the user to define the value [mm], direction (up/down) and speed
%   [mm/min] of the displacement. Three speeds are available: high,
%   corresponding to a speed of 9.5 mm/min, medium, corresponding to a
%   speed of 6.09 mm/min, and low, corresponding to a speed of 2.864
%   mm/min.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('---------------------- ADJUST SYRINGE VERTICAL POSITION ------------------------\n');
%{
fprintf('Possíveis procedimentos:\n');
fprintf('- Verifique o posicionamento da câmera e do estereomicroscópio; \n');
fprintf('-- A agulha deve estar verticalmente alinhada com o campo de visão da câmera e apontada para baixo; \n')
fprintf('- Verifique o posicionamento da agulha; \n');
fprintf('-- Posicione a agulha de forma que durante a medição a gota esteja centralizada na imagem. \n')
%}

% Preview do vídeo da câmera
%{
f = figure('Name', 'Live video','NumberTitle','off'); %Set live video figure
uicontrol('String', 'Close', 'Callback', 'close(gcf)'); %Create close button
vidResValues = vidobj.VideoResolution; %Get resolution of the video object
vidSizeRatio = vidResValues(2)/vidResValues(1); %Calculate video object size ratio
nBands = vidobj.NumberOfBands;
fsize = 1000; %Width of the video that will appear on the figure
f.Position(3:4) = [fsize fsize*vidSizeRatio]; %Redimension live video figure
movegui('center') %Centralize live video figure in the screen
hImage = image( zeros(vidResValues(2), vidResValues(1), nBands) ); 
preview(vidobj, hImage); %Preview live video object
% -Criar grida sobre o vídeo
rows = vidResValues(2);
columns = vidResValues(1);
hold on;
stepSizeY = vidResValues(2)/2; % Distance between lines in Y direction
stepSizeX = vidResValues(1)/2; %Distance between lines in X direction
for row = 1 : stepSizeY : rows
    yline(row, '-.r', 'LineWidth', 1);
end
for col = 1 : stepSizeX : columns
    xline(col, '-.r', 'LineWidth', 1);
end
%}

%Rotina de ajuste do posisionamento da seringa
%txt = "Y";
while 1
    resp = input("Do you like to proceed with the adjustment of the syringe vertical position (Y/N) [N]? ","s");
    if resp == "Y"
        deslZmm = input("Displacement (mm): ");
        sentidoDeslZ = input("Direction (up/down) [up]: ", "s");
        if isempty(sentidoDeslZ)
            sentidoDeslZ = "up";
        end
        velDeslZ = input("Speed (low/medium/high) [low]: ", "s");
        if isempty(velDeslZ)
            velDeslZ = "low";
        end
        fprintf('Please wait, moving syringe... \n');
        moverSeringa(passoMot,passoFuso,smSeringaZ,sentidoDeslZ,deslZmm,velDeslZ); %Moving the syringe
        release(smSeringaZ); %Motor release
        fprintf('Syringe displacement completed! \n');
    else
        break;
    end
end
end

