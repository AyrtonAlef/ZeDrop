function [angle] = passo2inc(nPassos,passoMotor,nEntradasSemFim,nDentesCoroa)
%PASSO2INC Returns how much the stage is tilted given a number of motor 
% steps.
%   INPUT
% nPassos - Number of motor steps
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% nEntradasSemFim - Worm number of starts
% nDentesCoroa - Number of worm gear teeth
%   OUPUT
% angle - tilt to be performed on the stage [°]

passoModInc = 360*nEntradasSemFim/nDentesCoroa; %tilt stage pitch [°/rev]
angle = nPassos * passoModInc/passoMotor;
end

