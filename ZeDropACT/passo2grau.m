function [grau] = passo2grau(nPassos,nDentesCoroa,passoMotor)
%PASSO2GRAU Conversion of motor steps to degrees (°)
grau = (nPassos*360)/(nDentesCoroa*passoMotor);
end

