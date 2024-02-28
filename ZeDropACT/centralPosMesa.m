function [posMesaX,posMesaY] = centralPosMesa(smMesaX,smMesaY,posMesaX,posMesaY,passoMot,passoFuso)
%CENTRALPOSMESA Positioning stage centering (solid substrate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Centralização automática da mesa de posicionamento. Retorno da mesa de
%   posicionamento (substrato sólido) ao centro de cada eixo. De forma a 
%   prover um rápido retorno, uma alta velocidade de deslocamento é usada 
%   (9,5 mm/min). *Requer que o sistema de posicionamento da mesa esteja 
%   conectado.
%   Automatic centering of the positioning stage. Return of the positioning
%   stage (solid substrate) to the center of each axis. In order to provide
%   a quick return, a high travel speed is used (9.5 mm/min).*Requires the
%   connection of the positioning stage system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% smMesaX - Control of the X axis motor of the positioning stage
% smMesaY - Control of the Y axis motor of the positioning stage
% posMesaX - Current X axis position of stage in mm  
% posMesaY - Current Y axis position of stage in mm
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% passoFuso - Linear displacement of the spindle corresponding to one revolution of the spindle [mm/rev]
%   OUTPUT
% posMesaX - Current X axis position of stage in mm  
% posMesaY - Current Y axis position of stage in mm  

fprintf('---------------- CENTER STAGE (SOLID SUBSTRATE) ---------------- \n');
% Centering stage routine
posXlimPos = 8; %Positive position limit on the X axis
posXlimNeg = -8; %Negative position limit on the X axis
posYlimPos = 24; %Positive position limit on the Y axis
posYlimNeg = -24; %Negative position limit on the Y axis

% Calculation of displacements
deslXmm = posMesaX*(-1);
deslYmm = posMesaY*(-1);

% Stage movement on the X axis
fprintf('Please wait, stage movement along the X axis... \n');
if deslXmm >= 0
    sentido = "positivo";
elseif deslXmm < 0
    sentido = "negativo";
end
velDeslX = "high";
moverMesa(passoMot,passoFuso,smMesaX,sentido,abs(deslXmm),velDeslX); %Move stage
release(smMesaX); %Release motor
posMesaX = posMesaX + deslXmm; %Update of the current position of the stage on the X axis
fprintf('Stage movement along the X axis completed! \n');
% Stage movement on the Y axis
fprintf('Please wait, stage movement along the Y axis... \n');
if deslYmm >= 0
    sentido = "positivo";
elseif deslYmm < 0
    sentido = "negativo";
end
velDeslY = "high";
moverMesa(passoMot,passoFuso,smMesaY,sentido,abs(deslYmm),velDeslY); %Move stage
release(smMesaY); %Release motor
posMesaY = posMesaY + deslYmm; %Update of the current position of the stage on the Y axis
fprintf('Stage movement along the X axis completed! \n');
pause(2);
end

