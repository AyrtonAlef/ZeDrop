function [nPassos] = inc2passo(angle,passoMotor,nEntradasSemFim,nDentesCoroa)
%INC2PASSO Returns how many motor steps are needed (nPassos) to perform a 
% certain tilt on the stage [°] (angleInc)
%   INPUT
% angle - Tilt to be perfomed on the stage [°]
% passoMotor - Number of steps required for one motor revolution [steps/rev]
% nEntradasSemFim - Worm number of starts
% nDentesCoroa - Number of worm gear teeth
%   OUPUT
% nPassos - Number of motor steps

passoModInc = 360*nEntradasSemFim/nDentesCoroa; %Tilting stage pitch [°/rev]
nPassos = angle * passoMotor/passoModInc;
end

