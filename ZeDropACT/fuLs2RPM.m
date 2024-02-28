function [RPM] = fuLs2RPM(taxaDepGota)
%FULS2RPM Returns the required RPM to be given to the motor to achieve a 
% desired drop flow rate (taxaDepGota).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Returns the required RPM to be given to the motor to achieve a desired
%   flow rate (uL/s). Use of cubic interpolation for curve fitting to
%   experimental data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- Experimental values collected by preliminary tests with the stepper motors
RPMExp = [2 5 10 15 20 30 40 60 80 100 120 150 200 300];
taxaDepExp = [0.315795 0.537468 0.731277 0.823086 0.872947 0.927602 0.961031 1.003375 1.027636 1.041647 1.049740 1.055950 1.059577 1.060731];

%- Fitting a fucntion to the experimental values collected
[f,gof] = fit(taxaDepExp', RPMExp','cubicinterp'); %Fit a pieciwise cubic interpolation
%{
%- Check fitting
gof
plot(f,taxaDepExp,RPMExp,'ro')
xlabel("taxaDep (uL/s)");
ylabel("RPM")
%}

%- Conversion of flow rate (uL/s) into motor input (RPM)
RPM = f(taxaDepGota);
end

