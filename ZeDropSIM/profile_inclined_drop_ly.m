function [X_data_cm,Z_data_cm,W,U,theta0,CA_calc_end,CA_adv_end,CA_rec_end,d_w_cm,drop_h_cm,area_cm2,volume_uL] = profile_inclined_drop_ly(b,c,alpha_deg,CA,N)
%PROFILE_INCLINED_DROP_LY Calculation of inclined drop profile by solving
%PDE derived from Young-Laplace equation using a numerical solution

%   INPUT
% b - curvature at the apex in cm^-1
% c - capillary constant in cm^-2
% alpha_deg - inclination angle in degrees
% N - number of elements along theta/ Suggestion: give N odd numbers, so the middle element is given at theta = 90°
%   OUTPUT
% theta0 - theta position where phi = CA at a W distance from the apex of the drop in radians
% CA_calc - intrinsic CA calculated in degrees
% CA_adv - advancing CA in degrees3
% CA_rec - receding CA in degrees
% area - dimensionless surface area of the inclined drop profile
% volume - dimensionless volume of the inclined drop profile
% W_data - dimensionless axial distance of the inclined drop profile
% U_data - dimensionless radius of the inclined drop profile

%INITIALIZATION
%- Liquid properties
R0 = 1/b; %Radius of curvature at the apex in cm;
beta = (R0^2)*c; %Bond Number beta = R0^2*g*delta_ro/gamma or beta = (R0^2)*c

%- Configuration parameters
alpha = deg2rad(alpha_deg); %Inclination angle in radians

%- Parameters
W_min = 0;
W_max = 6;
%delta_W = 0.0001;
delta_W_ini = 1e-6; %Highest precision near the apex
W_ini = 1e-2; %Axial distance from the apex until the drop profile is calculated with higher precision
delta_W_aft = 1e-4; %Precision for all the rest of the drop
%W = W_min:delta_W:W_max; %Dimensionless axial distance
delta_W = delta_W_ini;
W = zeros(1,round(W_max/delta_W)+1);
theta_min = 0;
theta_max = pi;
theta = linspace(theta_min,theta_max,N); %Angle between U and the X axis in radians
%theta_deg = rad2deg(theta); %theta vector in degrees
delta_theta = theta(end)/(length(theta)-1);
tol_CA = 1e-3; %tolerance for the precision in determining the inclination drop profile according to the apparent CA

%- Variables
U = zeros(length(W),length(theta)); %Dimensionless radius
V = zeros(length(W),length(theta)); %Dimensionless V = 𝜕U/𝜕W
Q = zeros(length(W),length(theta)); %Dimensionless Q = 𝜕U/𝜕theta
h = zeros(length(W),length(theta)); %Parameter to chase for gravity perpendicularity
phi = zeros(length(W),length(theta)); % Angle of intersection between the fluid interface and the solid in radians
theta0 = zeros(length(W),1); %Theta position where phi = CA at a W distance from the apex of the drop in radians
phi_theta0 = zeros(length(W),1); %phi at theta0 position in degrees
CA_calc = zeros(length(W),1); %Intrinsic apparent CA for a W distance from the apex of the drop in degrees
CA_adv = zeros(length(W),1); %Advancing CA in degrees
CA_rec = zeros(length(W),1); %Receding CA in degrees
%s = zeros(1,length(theta)); %Dimensionless arc length between W step for each theta position
area = zeros(length(W),1); %Dimensionless area;
volume = zeros(length(W),1); %Dimensionless volume

%Initial and boundary conditions
%- "Derivatives of both U and V with respect to theta equal zero at 0 and pi" 
%   (Larkin, 1967)
Q(:,1) = 0;
Q(:,end) = 0;
%- "Computations were started at k = 2, not k = 1. For the region between 
%   zero and one (the first delta_W increment) the surface was assumed to 
%   be spherical. Starting values along the line k = 2, as obtained from 
%   the spherical approximation"(Larkin,1967)
%-- Filling values for k = 2
W(2) = W(1) + delta_W; 
U(2,:) = (2*delta_W)^.5;
V(2,:) = (1-delta_W)./U(2,:);
Q(2,:) = 0;
h(2,:) = cos(theta) + (sin(theta)./U(2,:)).*Q(2,:);
phi(2,:) = (pi/2) + atan(V(2,:).*(1 + (Q(2,:)./U(2,:)).^2).^-.5);
indx_theta0 = find(abs(h(2,:)) == min(abs(h(2,:))),1,'first'); %Find indx where h = 0
theta0(2) = theta(indx_theta0); %Determining theta0
phi_theta0(2) = rad2deg(phi(2,indx_theta0)); %Determining phi at theta0 in degrees
CA_calc(2) = 180 - phi_theta0(2); %Determination of the apparent CA 
CA_adv(2) = CA_calc(2);
CA_rec(2) = CA_calc(2);
area(2) = pi*(mean(U(2,:))^2 + W(2)^2); %Dimensionless area assuming sphericity
volume(2) = (1/6)*pi*W(2)*(3*mean(U(2,:))^2+W(2)^2); %Dimensionless volume assuming sphericity

%- Determining inclined drop profile 
k = 2;
diff_CA = 1;
flag = 0; %flag for changing delta_W after reaching an axial distance W_ini from the apex 
while abs(diff_CA) > tol_CA
%}
%for k = 2:(length(W)-1)
    R = zeros((length(theta)-2),1); %Results matrix
    %R1 = zeros((length(theta)-2),1); %Results matrix
    M = zeros(length(theta)-2); %Coefficients matrix
    for j = 2:(length(theta)-1)
        %Nonlinear coefficients of PDE
        A = 1 + ((Q(k,j)^2)/(U(k,j)^2));
        B = -(V(k,j)*Q(k,j))/(U(k,j)^2);
        C = (1 + V(k,j)^2)/(U(k,j)^2); %(Larkin,1967)
        %C = 1 + ((V(k,j)^2)/(U(k,j)^2)); %(Xu and Wang, 2015)
        D = ((2 + beta*(U(k,j)*sin(alpha)*cos(theta(j))- W(k)*cos(alpha))) * (1+(V(k,j)^2)+(Q(k,j)^2)/(U(k,j)^2))^(3/2)) - (1/U(k,j))*(1+(V(k,j)^2)+2*((Q(k,j)^2)/(U(k,j)^2)));
        %{
        D11 = 2 + beta*(U(k,j)*sin(alpha)*cos(theta(j)) - W(k)*cos(alpha));
        D12 = (1 + V(k,j)^2 + ((Q(k,j)^2)/(U(k,j)^2)))^(3/2);
        D13 = (-1/U(k,j)) * (1 + V(k,j)^2 + 2*((Q(k,j)^2)/(U(k,j)^2))); 
        D1(k,j) = D11 * D12 + D13;
        %}
        
        %Creation of results matrix R
        %R(j-1,1) = (V(k,j)*(A/delta_W)) - D - (C/(2*delta_theta))*(Q(k,j+1)-Q(k,j-1)); 
        R(j-1,1) = (A/delta_W)*V(k,j) - (C/(2*delta_theta))*(Q(k,j+1)-Q(k,j-1)) - D;
        %{
        R11 = (A/delta_W)*V(k,j);
        R12 = (C/(2*delta_theta))*(Q(k,j+1)-Q(k,j-1));
        R13 = D;
        R1(j-1,1) = R11 - R12 - R13; 
        %}
        %R(k,j-1) = (A(k,j)/delta_W)*V(k,j) - (C(k,j)/(2*delta_theta))*(Q(k,j+1)-Q(k,j-1)) - D(k,j);
        %{
        R11 = (A(k,j)/delta_W)*V(k,j);
        R12 = (C(k,j)/(2*delta_theta))*(Q(k,j+1)-Q(k,j-1));
        R13 = D1(k,j);
        R1(k,j-1) = R11 - R12 - R13; 
        %}
        
        %Creation of coefficients matrix M
        %- 1st method
        %{
        if j == 2
            M(j-1,j-1) = (A/delta_W)-(B/delta_theta); %Main diagonal
            M(j-1,j) = (B/delta_theta); %Diagonal above main diagonal
        elseif j == (length(theta)-1)
            M(j-1,j-1) = (A/delta_W)+(B/delta_theta); %Main diagonal
            M(j-1,j-2) = -(B/delta_theta); %Diagonal below main diagonal
        else
            M(j-1,j-1) = (A/delta_W); %Main diagonal
            M(j-1,j) = (B/delta_theta); %Diagonal above main diagonal
            M(j-1,j-2) = -(B/delta_theta); %Diagonal below main diagonal
        end
        %}
        %{
        if j == 2
            M(j-1,j-1) = (A(k,j)/delta_W)-(B(k,j)/delta_theta); %Main diagonal
            M(j-1,j) = (B(k,j)/delta_theta); %Diagonal above main diagonal
        elseif j == (length(theta)-1)
            M(j-1,j-1) = (A(k,j)/delta_W)+(B(k,j)/delta_theta); %Main diagonal
            M(j-1,j-2) = -(B(k,j)/delta_theta); %Diagonal below main diagonal
        else
            M(j-1,j-1) = (A(k,j)/delta_W); %Main diagonal
            M(j-1,j) = (B(k,j)/delta_theta); %Diagonal above main diagonal
            M(j-1,j-2) = -(B(k,j)/delta_theta); %Diagonal below main diagonal
        end
        %}
        %- 2nd method
        %%{
        if j == 2 %first line
            M(j-1,j-1) = (A/delta_W) - (4/3)*(B/delta_theta); %Main diagonal
            M(j-1,j) = (4/3)*(B/delta_theta); %Diagonal above main diagonal
        elseif j == (length(theta)-1) %last line
            M(j-1,j-1) = (A/delta_W) + (4/3)*(B/delta_theta); %Main diagonal
            M(j-1,j-2) = -(4/3)*(B/delta_theta); %Diagonal below main diagonal
        else
            M(j-1,j-1) = (A/delta_W); %Main diagonal
            M(j-1,j) = (B/delta_theta); %Diagonal above main diagonal
            M(j-1,j-2) = -(B/delta_theta); %Diagonal below main diagonal
        end
        %}
        %{
        if j == 2 %first line
            M(j-1,j-1) = (A(k,j)/delta_W) - (4/3)*(B(k,j)/delta_theta); %Main diagonal
            M(j-1,j) = (4/3)*(B(k,j)/delta_theta); %Diagonal above main diagonal
        elseif j == (length(theta)-1) %last line
            M(j-1,j-1) = (A(k,j)/delta_W) + (4/3)*(B(k,j)/delta_theta); %Main diagonal
            M(j-1,j-2) = -(4/3)*(B(k,j)/delta_theta); %Diagonal below main diagonal
        else
            M(j-1,j-1) = (A(k,j)/delta_W); %Main diagonal
            M(j-1,j) = (B(k,j)/delta_theta); %Diagonal above main diagonal
            M(j-1,j-2) = -(B(k,j)/delta_theta); %Diagonal below main diagonal
        end
        %}
    end
    %{
    % Compute M-inverse
    M_inv = inv(M);
    % Left-multiply R by M-inverse. S is the solution vector
    S = M_inv*R;
    %}
    %Resolution of the system of equations. S is the solution vector
    S = M\R;
    %Prova real
    %S = M\R(k,:)';
    V(k+1,2:(length(theta)-1)) = S';
    
    % Boundary conditions/ Filling V values for theta = 0 (j=1) and theta = pi (j=N)
    %- 1st method
    %{
    V(k+1,1) = V(k+1,2);
    V(k+1,end) = V(k+1,end-1);
    %}
    %- 2nd method
    %%{
    V(k+1,1) = (4*V(k+1,2)-V(k+1,3))/3;
    V(k+1,end) = (4*V(k+1,end-1)-V(k+1,end-2))/3;
    %}
    
    %Calculating U and Q values
    for j = 2:length(theta)-1
        Q(k+1,j) = (delta_W/(2*delta_theta))*(V(k+1,j+1)-V(k+1,j-1)) + Q(k,j);
        U(k+1,j) = (delta_W/2)*(V(k+1,j)+V(k,j)) + U(k,j);
    end
    %- 1st method
    %{
    U(k+1,1) = U(k+1,2);
    U(k+1,end) = U(k+1,end-1);
    %}
    %- 2nd method
    %%{
    U(k+1,1) = (4*U(k+1,2)-U(k+1,3))/3;
    U(k+1,end) = (4*U(k+1,end-1)-U(k+1,end-2))/3;
    %}
    
    %Calculating h and phi values
    h(k+1,:) = cos(theta) + (sin(theta)./U(k+1,:)).*Q(k+1,:);
    phi(k+1,:) = (pi/2) + atan(V(k+1,:).*(1 + (Q(k+1,:)./U(k+1,:)).^2).^-.5);
    %{
    for j=1:length(theta)
        h1(k+1,j) = cos(theta(j)) + (sin(theta(j))/U(k+1,j))*Q(k+1,j);
        phi1(k+1,j) = (pi/2) + atan(V(k+1,j)*(1 + (Q(k+1,j)/U(k+1,j))^2)^(-1/2));
    end
    %}
    
    %Calculation of surface area and volume
    s = ((delta_W)^2 + (U(k+1,:)-U(k,:)).^2).^.5;
    sum_A = 0;
    sum_V = 0;
    for j=1:(length(theta)-1)
        sum_A = sum_A + ((s(j)+s(j+1)) * (U(k+1,j)+U(k+1,j+1)));
        sum_V = sum_V + ((U(k+1,j)+U(k+1,j+1))^2);
    end
    delta_A = (delta_theta/2) * sum_A; %Increment of area
    area(k+1) = area(k) + delta_A;
    delta_V = ((delta_W*delta_theta)/4) * sum_V;
    volume(k+1) = volume(k) + delta_V;
    
    %Determining theta position and phi where the line of solid-liquid-vapor contact moves perpendicular to g
    indx_theta0 = find(abs(h(k+1,:)) == min(abs(h(k+1,:))),1,'first'); %Find indx where h = 0
    theta0(k+1) = theta(indx_theta0); %Determining theta0
    phi_theta0(k+1) = rad2deg(phi(k+1,indx_theta0)); %Determining phi at theta0 in degrees
    CA_calc(k+1) = 180 - phi_theta0(k+1); %Determination of the apparent CA 
    CA_adv(k+1) = 180 - rad2deg(phi(k+1,1)); %Determination of the advancing CA
    CA_rec(k+1) = 180 - rad2deg(phi(k+1,end)); %Determination of the receding CA
    
    %%{
    diff_CA = CA - CA_calc(k+1);
    if diff_CA < 0 %CA_calc > CA
        W(k+1) = W(k) + delta_W;
        delta_W = delta_W/2; %Reduce of W step
        diff_CA = 1;
    else
        W(k+1) = W(k) + delta_W;
        k = k +1;
    end
    
    if W(k) >= W_ini && flag == 0 %Changing delta_W after W_ini is reached only one time
        delta_W = delta_W_aft;
        flag = 1;
    end
    %}
end
%%{
%Trim variables according to CA
W = W(1,1:k);
volume = volume(1:k,1);
area = area(1:k,1);
U = U(1:k,:);
phi = phi(1:k,:);
phi_theta0 = phi_theta0(1:k,:);
CA_calc = CA_calc(1:k,1);
CA_calc_end = CA_calc(end,1);
CA_adv = CA_adv(1:k,1);
CA_adv_end = CA_adv(end,1);
CA_rec = CA_rec(1:k,1);
CA_rec_end = CA_rec(end,1);
%}
%{
%Plotting inclined drop profile
plot(U(:,1),-W')
hold on
plot(-U(:,end),-W')
axis equal
%}
%Trim inclined drop profile for theta=0 and theta=pi
U_adv = U(:,end);
U_rec = U(:,1);

%Transformation to cartesian coordinates and conversion to dimensions
X_adv_cm = R0*U_adv*cos(theta(end)); %x coordinates of the advancing side of the drop in cm
X_rec_cm = R0*U_rec*cos(theta(1)); %x coordinates of the receding side of the drop in cm
Z_cm = R0*W'; %z coordinates in cm
area_cm2 = (R0^2)*area(end); %inclined drop surface in cm² 
volume_uL = ((R0^3)*volume(end))*1000; %inclined drop volume in mm³ or uL

%Build of the inclined drop profile
X_data_cm = [wrev(X_adv_cm); X_rec_cm(2:end)]; %x coordinates of the inclined drop profile in cm
Z_data_cm = [wrev(Z_cm); Z_cm(2:end)]; %z coordinates of the inclined drop profile in cm

%Drop properties
d_w_cm = X_data_cm(end) - X_data_cm(1); %distance between contact points in cm
drop_h_cm = max(abs(Z_data_cm)); %drop height in cm
end

