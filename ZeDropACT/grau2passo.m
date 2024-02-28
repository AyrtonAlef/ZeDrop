function [nPassos] = grau2passo(grau,nDentesCoroa,passoMotor)
%GRAU2PASSO Conversion of degrees (°) to motor steps.
nPassos = (nDentesCoroa*passoMotor*grau)/360;
end