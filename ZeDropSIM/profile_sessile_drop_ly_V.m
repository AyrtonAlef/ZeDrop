function [x_ly_sim_cm_otm,z_ly_sim_cm_otm,theta_otm,b_otm,r_w_otm,r_eq_otm,h_otm,A_otm,V_otm] = profile_sessile_drop_ly_V(V,c,theta,S_span,initialCo,bmin_ini,bmax_ini,interval_b,tolV,itermax)
%SESSDROPFINDBRELATEDTOV Calculation of sessile drop profile by the 
% numerical integration of the Young-Laplace equation given curvature at 
% the apex (b), liquid capillary constant (c) and contact angle (theta).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Return sessile drop important properties given liquid capillary
%   constant [cm-2](c), curvatures at apex [cm-1](b) and contact angle
%   [°](CA).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
% c - liquid capillary constant in cm-2
% V - array of drop volumes in cm-2
% theta - contact angle in degree
% S_span - step variable for ode45 solver
% initialCo - initial conditions for integration
% bmin_ini - minimum curvature at the apex that is going to be searched in cm-1
% bmax_ini - maximum curvature at the apex that is going to be searched in cm-1
% interval_b - number of b values created between bmin and bmax for volume search
% tolV - tolerance for volume determination in uL
% itermax - maximum number of iterations for volume search

%   OUPUT:
% x_ly_sim_cm - x coordinates of sessile drop profile in cm
% z_ly_sim_cm - z coordinates of sessile drop profile in cm
% theta_otm - contact angle effectively used for drop simulation in degree
% b_otm - curvatures at the apex effectuvely used for drop simulation in cm-1
% r0_otm - radius of curvature at the apex in cm
% r_w_otm - wetting radius in cm
% r_eq_otm - equatorial (maximum) radius of the drop in cm
% h_otm - drop height in cm
% A_otm - drop surface area in cm2
% V_otm - drop volume in uL

%{
%- Configurations for searching b values
bmin_ini = 0.05; %Minimum curvature at the apex that is going to be searched in cm-1
bmax_ini = 10; %Maximum curvature at the apex that is going to be searched in cm-1
interval_b = 50; %Number of b values created between bmin and bmax for volume search
tolV = 1e-3; %Tolerance for volume determination in uL 
itermax = 5; %Maximum number of iterations for volume search
%}

%VARIABLES INITIALIZATION
%x_ly_sim_cm = NaN(1e6,interval_b); 
%z_ly_sim_cm = NaN(1e6,interval_b);
x_ly_sim_cm = cell(1,interval_b);
z_ly_sim_cm = cell(1,interval_b);
theta_sim = zeros(1,interval_b); %Contact angle (theta) effectively used for drop simulation
r_w_cm = zeros(1,interval_b);  %Wetting radius in cm
r_eq_cm = zeros(1,interval_b); %Equatorial (maximum) radius of the drop in cm
drop_height_cm = zeros(1,interval_b); %Drop height in cm
A_cm2 = zeros(1,interval_b); %Drop surface area in cm²
V_uL = zeros(1,interval_b); %Drop volume in uL calculated concomitantly with the Laplace curve generation
%V_otm = length(V);
%x_ly_sim_cm_otm = cell(1,length(V));
%z_ly_sim_cm_otm = cell(1,length(V));
%theta_otm = length(V);
%b_otm = length(V);
%r0_otm = length(V);
%r_w_otm = length(V);
%r_eq_otm = length(V);
%h_otm =  length(V);
%A_otm = length(V);

% Find b values related to V
%for j = 1:length(V)
    diffV = 1; %reset volume difference
    bmin = bmin_ini; %reset bmin
    bmax = bmax_ini; %reset bmax
    iter = 1;
    while diffV > tolV && iter <= itermax
        %fprintf('Searching for V = %.2f uL, between bmin = %.4f and bmax = %.4f (iteration %d)... \n',V(j),bmin,bmax,iter);
        fprintf('Searching for V = %.2f uL, between bmin = %.4f and bmax = %.4f (iteration %d)... \n',V,bmin,bmax,iter);
        b = linspace(bmin,bmax,interval_b);
        r0 = 1./b;
        for i = 1:length(b)
            %Determination of the sessile drop profile by numerical integration of the Young-Laplace equation
            [x_ly_sim_cm{i},z_ly_sim_cm{i},theta_sim(i),r_w_cm(i),r_eq_cm(i),drop_height_cm(i),A_cm2(i),V_uL(i)] = profile_sessile_drop_ly(b(i),c,theta,S_span,initialCo);
            %[x_ly_sim_cm_b,z_ly_sim_cm_b,theta_sim(i),r_w_cm(i),r_eq_cm(i),drop_height_cm(i),A_cm2(i),V_uL(i)] = profile_sessile_drop_ly(b(i),c,theta,S_span,initialCo);
            %x_ly_sim_cm{i} = x_ly_sim_cm_b; 
        end
        %diffV = min(abs(V_uL-V(j)));
        diffV = min(abs(V_uL-V));
        fprintf('diffV = %f\n',diffV);
        %indx_V = find(abs(V_uL-V(j)) == diffV,1,'last');
        indx_V = find(abs(V_uL-V) == diffV,1,'last');
        if (indx_V-1) <= 0
            bmin = b(indx_V);
            bmax = b(indx_V+1);
        elseif (indx_V+1) >= length(b)
            bmin = b(indx_V-1);
            bmax = b(indx_V);
        else
            bmin = b(indx_V-1);
            bmax = b(indx_V+1);
        end
        iter = iter +1;
    end
    %{
    V_otm(j) = V_uL(indx_V);
    x_ly_sim_cm_otm{j} = x_ly_sim_cm{indx_V};
    z_ly_sim_cm_otm{j} = z_ly_sim_cm{indx_V};
    %indx_Profile = ~isnan(x_ly_sim_cm(:,indx_V));
    %x_ly_sim_cm_otm(j) = x_ly_sim_cm(indx_Profile,indx_V);
    %indx_z_ly = ~isnan(z_ly_sim_cm(:,indx_V));
    %z_ly_sim_cm_otm(j) = z_ly_sim_cm(indx_Profile,indx_V);
    theta_otm(j) = theta_sim(indx_V);
    b_otm(j) = b(indx_V);
    r0_otm(j) = r0(indx_V);
    r_w_otm(j) = r_w_cm(indx_V);
    r_eq_otm(j) = r_eq_cm(indx_V);
    h_otm(j) =  drop_height_cm(indx_V);
    A_otm(j) = A_cm2(indx_V);
    %}
    V_otm = V_uL(indx_V);
    x_ly_sim_cm_otm = x_ly_sim_cm{indx_V};
    z_ly_sim_cm_otm = z_ly_sim_cm{indx_V};
    %indx_Profile = ~isnan(x_ly_sim_cm(:,indx_V));
    %x_ly_sim_cm_otm(j) = x_ly_sim_cm(indx_Profile,indx_V);
    %indx_z_ly = ~isnan(z_ly_sim_cm(:,indx_V));
    %z_ly_sim_cm_otm(j) = z_ly_sim_cm(indx_Profile,indx_V);
    theta_otm = theta_sim(indx_V);
    b_otm = b(indx_V);
    r0_otm = r0(indx_V);
    r_w_otm = r_w_cm(indx_V);
    r_eq_otm = r_eq_cm(indx_V);
    h_otm =  drop_height_cm(indx_V);
    A_otm = A_cm2(indx_V);
%end
end