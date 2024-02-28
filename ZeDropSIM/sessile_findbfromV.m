function [x_ly_sim_cm_otm,z_ly_sim_cm_otm,theta_otm,b_otm,r0_otm,r_w_otm,r_eq_otm,h_otm,A_otm,V_otm] = sessile_findbfromV(c,theta,V,bmin,bmax,interval_b,tolV,S_span,initialCo)
%SESSILE_FINDBFROMV Determines the apex curvature (b) and drop properties from drop volume (V)
%   INPUT:
% c - capillary constant for water in cm^-2
% theta - apparent CA in deg
% V - drop volume in uL
% bmin - minimum boundary of b (cm-1) that is gonna be searched
% bmax - maximum boundary of b (cm-1) that is gonna be searched
% interval_b - number of b values created between bmin and bmax for volume search
% tolV - tolerance for volume determination in uL 
% S_span - step variable for ode45 solver
% initialCo - initial conditions for integration

%VARIABLES INITIALIZATION 
x_ly_sim_cm = cell(1,interval_b); %x coordinates of drop profile in cm
z_ly_sim_cm = cell(1,interval_b); %z coordinates of drop profile in cm
theta_sim = zeros(1,interval_b); %theta effectively used for drop creation
r_w_cm = zeros(1,interval_b); %wetted radius in cm
%r_w_pxs = zeros(1,interval_b); %wetted radius in pxs
r_eq_cm = zeros(1,interval_b); %drop equatorial (maximum) radius in cm
%r_eq_pxs = zeros(1,interval_b); %drop equatorial (maximum) radius in cm
h_cm = zeros(1,interval_b); %drop height in cm
%h_pxs = zeros(1,interval_b); %drop height in pxs
A_cm2 = zeros(1,interval_b); %Drop surface area in cm²
V_uL = zeros(1,interval_b); %Drop volume in uL calculated concomitantly with the Laplace curve generation
%{
b_otm = length(V);
theta_otm = length(V);
r_w_otm = length(V);
r_eq_otm = length(V);
h_otm =  length(V);
A_otm = length(V);
V_otm = length(V);
%}
diffV = 1;
iter = 1;

%MAIN CODE
while diffV > tolV && iter <= 5
    fprintf('Searching for V = %.2f uL, between bmin = %.4f and bmax = %.4f... \n',V,bmin,bmax);
    b = linspace(bmin,bmax,interval_b);
    r0 = 1./b;
    for i=1:length(b)
        %Determination of the sessile drop profile by numerical integration of the Young-Laplace equation
        [x_ly_sim_cm{i},z_ly_sim_cm{i},theta_sim(i),r_w_cm(i),r_eq_cm(i),h_cm(i),A_cm2(i),V_uL(i)] = profile_sessile_drop_ly(b(i),c,theta,S_span,initialCo);
    end
    diffV = min(abs(V_uL-V));
    fprintf('diffV = %f\n',diffV);
    indx_V = find(abs(V_uL-V) == diffV,1,'last');
    if indx_V-1<=0
        bmin = b(indx_V);
        bmax = b(indx_V+1);
    elseif indx_V+1>=length(b)
        bmin = b(indx_V-1);
        bmax = b(indx_V);
    else
        bmin = b(indx_V-1);
        bmax = b(indx_V+1);
    end
    iter = iter +1;
end
x_ly_sim_cm_otm = x_ly_sim_cm{indx_V};
z_ly_sim_cm_otm = z_ly_sim_cm{indx_V};
theta_otm = theta_sim(indx_V);
b_otm = b(indx_V);
r0_otm = r0(indx_V);
r_w_otm = r_w_cm(indx_V);
r_eq_otm = r_eq_cm(indx_V);
h_otm =  h_cm(indx_V);
A_otm = A_cm2(indx_V);
V_otm = V_uL(indx_V);
end

