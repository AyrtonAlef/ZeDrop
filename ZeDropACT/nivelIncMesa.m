function [incMesaTot] = nivelIncMesa(smInc,incMesaTot,passoMot,nEntradasSemFim,nDentesCoroa)
%INCMESATOT Nivelamento da mesa (substrato sólido),
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Nivelamento automático da mesa. Retorno da mesa (substrato sólido) ao 
%   nível horizontal. De forma a prover um rápido retorno, uma alta taxa de
%   inclinação é usada (0,275 °/s). *Requer que o módulo de inclinação 
%   esteja conectado.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% smInc - Controle do motor de posicionamento conectado ao módulo de inclinação
% incMesaTot - Inclinação atual da mesa [°]
% passoMotor - Número de passos necessário para uma revolução no motor [passos/rev]
% nEntradasSemFim -Número de entradas do parafuso sem-fim
% nDentesCoroa - Número de dentes da coroa
%   OUTPUT
% incMesaTot - Inclinação atual da mesa[°]

fprintf('-------------- NIVELAMENTO DA MESA (SUBSTRATO SÓLIDO) -------------- \n');
% Rotina de nivelamento da mesa
% - Cálculo da inclinação
angle = incMesaTot*(-1);

%- Inclinação da mesa
fprintf('Aguarde, inclinação da mesa... \n');
if angle >= 0
    sentido = "AH";
elseif angle < 0
    sentido = "H";
end
velocidade = "alta";
moverIncMesa(smInc,passoMot,nEntradasSemFim,nDentesCoroa,angle,sentido,velocidade); %Inclinando mesa
release(smInc); %Liberação do motor
incMesaTot = incMesaTot + angle; %Atualização da posição atual da mesa no eixo X;
fprintf('Fim da inclinação da mesa! \n');
pause(2);
end

