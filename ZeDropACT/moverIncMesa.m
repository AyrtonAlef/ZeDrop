function  moverIncMesa(smInc,passoMotor,nEntradasSemFim,nDentesCoroa,alpha,sentido,velocidade)
%MOVERINCMESA Tilt stage according to a predetermined angle (angle), direction (sentido) and speed (velocidade)
%   INPUT
% smInc - Control of the motor connected to the tilting stage
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% nEntradasSemFim - Worm number of starts
% nDentesCoroa - Number of worm gear teeth
% angle - tilt angle [°]
% sentido - Stage tilt direction (clockwise - CW/ counterclockwise - CCW)
% velocidade - tilt rate (low/medium/high)

nPassos = inc2passo(alpha,passoMotor,nEntradasSemFim,nDentesCoroa); %Conversion of tilt in ° to the necessary number of motor steps
if strcmp(sentido,'CCW') %Determination of tilt direction
    i = 1;
elseif strcmp(sentido,'CW')
    i = -1;
end
if strcmp(velocidade,'low') %Determination of tilt rate
    smInc.RPM = 2; %RPM = 2 -> wreal = 1,432 RPM -> vbaixa = 0,0826 °/s
elseif strcmp(velocidade,'medium') 
    smInc.RPM = 8; %RPM = 8 -> wreal = 3,045 RPM -> vmedia = 0,1760 °/s
elseif strcmp(velocidade,'high')
    smInc.RPM = 158; %RPM = 125 -> wreal = 4,767 RPM -> vmedia = 0,2750 °/s
end
move(smInc,i*round(nPassos)); %Move stepper motor
end



