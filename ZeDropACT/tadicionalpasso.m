function [tadpasso] = tadicionalpasso(passoMotor,wMotor,wMin)
%TADICIONALPASSO Returns the delay time required per step to reach the 
% desired rotational speed.
%   INPUT
% passoMotor - Number of motor steps
% wMotor - Motor rotation speed to be achieved [RPM]
% wMin - Minimum stable rotational speed that the motor can execute [RPM]
%   OUTPUT
% tadpasso - Time that must be added between each step to reach the desired rotational speed [s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Returns the necessary delay time per step (tadpasso) to reach the
%   necessary rotational speed in the motor (wMotor) corresponding to the
%   desired flow rate. It is only relevant when the minimum motor rotation
%   speed (wMin) is greater than the speed corresponding to the desired
%   flow rate (wMotor).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tadpasso = (60/passoMotor)*((1/wMotor)-(1/wMin));
end

