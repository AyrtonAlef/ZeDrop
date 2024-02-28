function [wMotor] = incs2RPM(taxaInc,nEntradasSemFim,nDentesCoroa)
%INCS2RPM Returns the motor rotation speed [RPM] corresponding to a given 
% tilt rate (taxaInc) [°/s].
%   INPUT
% taxaInv - Tilt rate [°/s]
% nEntradasSemFim - Worm number of starts
% nDentesCoroa - Number of worm gear teeth
%   OUTPUT
% wMotor - motor rotation speed [RPM]

wMotor = (taxaInc*nDentesCoroa*60)/(360*nEntradasSemFim);
end