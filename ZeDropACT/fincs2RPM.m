function [RPM] = fincs2RPM(taxaInc,nEntradasSemFim,nDentesCoroa)
%FINCS2RPM Returns the required RPM to be sent to the motor to achieve a 
% desired stage tilt rate (taxaInc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Returns the required RPM to be sent to the motor to achieve a desired 
%   stage tilt rate (°/s). Use of a cubic interpolation to fit the curve to
%   the experimental data. Valid only for tilt rates between the maximum
%   and minimum motor rotation speed limits (wmin = 1.4171 RPM and wmax = 
%   4.8526 RPM).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- Experimental values collected by preliminary tests with the stepper motors
RPMExp = [2 5 10 15 20 30 40 60 80 100 120 150 200 300];
RPMReal = [1.4322 2.4375 3.3165 3.7328 3.9589 4.2068 4.3584 4.5505 4.6605 4.7240 4.7607 4.7889 4.8053 4.8106];
%taxaDepExp = [0.315795 0.537468 0.731277 0.823086 0.872947 0.927602 0.961031 1.003375 1.027636 1.041647 1.049740 1.055950 1.059577 1.060731];
taxaIncExp = RPMReal.*(nEntradasSemFim*6/nDentesCoroa);

%- Fitting a fucntion to the experimental values collected
[f,gof] = fit(taxaIncExp', RPMExp','cubicinterp'); %Fit a pieciwise cubic interpolation
%{
%- Check fitting
gof
plot(f,taxaIncExp,RPMExp,'ro')
xlabel("taxaInc (°/s)");
ylabel("RPM")
%}

%- Conversion of the tilt rate (°/s) into input for the motor (RPM)
RPM = f(taxaInc);
end

