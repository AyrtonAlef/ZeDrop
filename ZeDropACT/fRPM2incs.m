function [RPM] = fRPM2incs(taxaInc,nEntradasSemFim,nDentesCoroa)
%FRPM2INCS Returns the stage tilt rate (taxaInc) corresponding to a given
% motor RPM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Returns the stage tilt rate (°/s) corresponding to a given motor RPM.
%   Use of cubic interpolation for curve fitting to experimental data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- Experimental values collected by preliminary tests with the stepper motors
RPMExp = [2 5 10 15 20 30 40 60 80 100 120 150 200 300];
RPMReal = [1.4322 2.4375 3.3165 3.7328 3.9589 4.2068 4.3584 4.5505 4.6605 4.7240 4.7607 4.7889 4.8053 4.8106];
%taxaDepExp = [0.315795 0.537468 0.731277 0.823086 0.872947 0.927602 0.961031 1.003375 1.027636 1.041647 1.049740 1.055950 1.059577 1.060731];
taxaIncExp = RPMReal.*(nEntradasSemFim*6/nDentesCoroa);

%- Fitting a fucntion to the experimental values collected
[f,gof] = fit(RPMExp',taxaIncExp', 'cubicinterp'); %Fit a pieciwise cubic interpolation
%{
%- Check fitting
gof
plot(f,RPMExp,taxaIncExp,'ro')
ylabel("taxaInc (°/s)");
xlabel("RPM")
%}

%- Conversion of motor RPM to tilt rate (°/s)
taxaInc = f(RPM);
end

