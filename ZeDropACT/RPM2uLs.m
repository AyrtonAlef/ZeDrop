function [taxaDep] = RPM2uLs(wMotor,passoSeringa)
%RPM2ULS Conversion of RPM to uL/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Returns the flow rate (uL/s) corresponding to a motor rotation speed
%   (RPM).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taxaDep = wMotor*passoSeringa/60;
end

