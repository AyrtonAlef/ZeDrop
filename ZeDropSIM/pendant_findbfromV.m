function [x_ly_sim_cm_otm,z_ly_sim_cm_otm,b_otm,r0_otm,d_eq_otm,h_eq_otm,h_otm,d_s_otm,A_otm,V_otm] = pendant_findbfromV(c,r_h_cm,tol,V,bmin,bmax,interval_b,tolV,S_span,initialCo)
%PENDANT_FINDBFROMV Determines the apex curvature (b) and drop properties of pendant drop from drop volume (V)
%   INPUT:
% c - capillary constant for water in cm^-2
% r_h_cm - holder outer radius in cm
% tol - tolerance to find the intersection points of the drop with the holder radius
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
d_eq_cm = zeros(1,interval_b); %Drop maximum diameter in cm
h_eq_cm = zeros(1,interval_b); %Distance from the drop apex to the horizontal position of drop maximum diameter in cm
d_s_cm = zeros(1,interval_b); %Drop diameter at a distance of d_eq from the drop apex in cm
h_cm = zeros(1,interval_b); %Distance from the drop apex to the needle tip in cm
A_cm2 = zeros(1,interval_b); %Drop surface area in cm²
V_uL = zeros(1,interval_b); %Drop volume in uL calculated concomitantly with the Laplace curve generation

diffV = 1;
iter = 1;
%MAIN CODE
while diffV > tolV && iter <= 5
    fprintf('Searching for V = %.2f uL, between bmin = %.4f and bmax = %.4f... \n',V,bmin,bmax);
    b = linspace(bmin,bmax,interval_b);
    r0 = 1./b;
    for i=1:length(b)
        %Determination of the pendant drop profile by numerical integration of the Young-Laplace equation
        %[x_ly_sim_cm{i},z_ly_sim_cm{i},theta_sim(i),r_w_cm(i),r_eq_cm(i),h_cm(i),A_cm2(i),V_uL(i)] = profile_sessile_drop_ly(b(i),c,theta,S_span,initialCo);
        [x_ly_cm,z_ly_cm,x_ly_sim_cm{i},z_ly_sim_cm{i},indx_last,indx_second,indx_deq,indx_ds,d_eq_cm(i),h_eq_cm(i),h_cm(i),d_s_cm(i),A_cm2(i),V_uL(i)] = profile_pendant_drop_ly(b(i),c,r_h_cm,S_span,initialCo,tol);
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
b_otm = b(indx_V);
r0_otm = r0(indx_V);
d_eq_otm = d_eq_cm(indx_V);
d_s_otm = d_s_cm(indx_V);
h_eq_otm = h_eq_cm(indx_V);
h_otm =  h_cm(indx_V);
A_otm = A_cm2(indx_V);
V_otm = V_uL(indx_V);
end

