function [nPassos] = mm2passo(passoMot,passoFuso,deslmm)
%PASSOS2VOL Returns the motor steps required for a given linear spindle 
% displacement [mm] 
%   INPUT
% passoMotor - Number of steps required for one motor revolution [steps/rev]
% passoFuso - Linear displacement of the spindle corresponding to one revolution of the spindle [mm/rev]
% deslmm - Displacement [mm]
%   OUTPUT
% nPassos - Number of motor steps

nPassos = (passoMot/passoFuso)*deslmm;
end

