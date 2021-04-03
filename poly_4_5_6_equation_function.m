function [S, V, A, J, s, v, a, j, theta, C_r, C_f, x_rf, y_rf, Rp, phi, x_fff, y_fff, Rb_fff]= poly_4_5_6_equation_function(beta_r,beta_r_deg,...
    beta_f,beta_f_deg,beta_d, beta_d_deg, h,Rf, Rb, omega)
%% Calculating SVAJ for the Rise

coeff_r = [1, 0, 0, 0, 0, 0, 0, 0; %s start
           1, 1, 1, 1, 1, 1, 1, 1; %s apex
           0, (1/beta_r_deg), 0, 0, 0, 0, 0, 0; %v start
           0, (1/beta_r_deg), (2/beta_r_deg), (3/beta_r_deg), (4/beta_r_deg), (5/beta_r_deg), 6/beta_r_deg, 7/beta_r_deg; %v apex
           0, 0, 2/(beta_r_deg)^2, 0, 0, 0, 0, 0; %a start
           0, 0, 2/(beta_r_deg)^2, 6/(beta_r_deg)^2, 12/(beta_r_deg)^2, 20/(beta_r_deg)^2, 30/(beta_r_deg)^2, 42/(beta_r_deg)^2; %a apex
           0, 0, 0, 6/(beta_r_deg)^3, 0, 0, 0, 0; %j start
           0, 0, 0, 6/(beta_r_deg)^3, 24/(beta_r_deg)^3, 60/(beta_r_deg)^3, 120/(beta_r_deg)^3, 210/(beta_r_deg)^3;]; %j apex

bc_r = [0, h, 0, 0, 0, 0, 0, 0]';

%solve for coeffecients
C_r = linsolve(coeff_r, bc_r);

%define theta for the rise
theta_r = linspace(0, beta_r_deg,beta_r_deg);

x_r = theta_r./beta_r_deg;

s_r = C_r(1)+C_r(2)*(x_r)+C_r(3)*(x_r).^2+C_r(4)*(x_r).^3+C_r(5)*(x_r).^4 ...
    +C_r(6)*(x_r).^5+C_r(7)*(x_r).^6+C_r(8)*(x_r).^7;

v_r = (1/beta_r).*(C_r(2)+2*C_r(3)*x_r+3*C_r(4)*(x_r).^2+4*C_r(5)*(x_r).^3 ...
    +5*C_r(6)*(x_r).^4+6*C_r(7)*(x_r).^5+7*C_r(8)*(x_r).^6);

a_r = (1/beta_r)^2.*(2*C_r(3)+6*C_r(4)*x_r+12*C_r(5)*(x_r).^2+20*C_r(6)*(x_r).^3 ...
    +30*C_r(7)*(x_r).^4+42*C_r(8)*(x_r).^5);

j_r = (6*C_r(4)+24*C_r(5)*(x_r)+60*C_r(6)*(x_r).^2+120*C_r(7)*(x_r).^3 ...
    +210*C_r(8)*(x_r).^4);

%% Calculating SVAJ for the Fall

coeff_f = [1, 0, 0, 0, 0, 0, 0, 0; %s apex
           1, 1, 1, 1, 1, 1, 1, 1; %s fall end
           0, 1/beta_f_deg, 2/beta_f_deg, 3/beta_f_deg, 4/beta_f_deg, 5/beta_f_deg, 6/beta_f_deg, 7/beta_f_deg; %v apex
           0, 1/beta_f_deg, 0, 0, 0, 0, 0, 0; %v fall end
           0, 0, 2/(beta_f_deg)^2, 0, 0, 0, 0, 0; %a start
           0, 0, 2/(beta_f_deg)^2, 6/(beta_f_deg)^2, 12/(beta_f_deg)^2, 20/(beta_f_deg)^2, 30/(beta_f_deg)^2, 42/(beta_f_deg)^2; %a fall
           0, 0, 0, 6/(beta_f_deg)^3, 0, 0, 0, 0; %j apex
           0, 0, 0, 6/(beta_f_deg)^3, 24/(beta_f_deg)^3, 60/(beta_f_deg)^3, 120/(beta_f_deg)^3, 210/(beta_f_deg)^3]; %j fall end

%Defining boundary conditions
bc_f = [h, 0, 0, 0, 0, 0, 0, 0]';

%solving for coeffecient values
C_f = linsolve(coeff_f,bc_f);

%defining theta for the fall
theta_f = linspace(0,beta_f_deg,beta_f_deg);

x_f = theta_f./beta_f_deg;

s_f = C_f(1)+C_f(2)*(x_f)+C_f(3)*(x_f).^2+C_f(4)*(x_f).^3+C_f(5)*(x_f).^4 ...
    +C_f(6)*(x_f).^5+C_f(7)*(x_f).^6+C_f(8)*(x_f).^7;

v_f = (1/beta_f).*(C_f(2)+2*C_f(3)*x_f+3*C_f(4)*(x_f).^2+4*C_f(5)*(x_f).^3 ...
    +5*C_f(6)*(x_f).^4+6*C_f(7)*(x_f).^5+7*C_f(8)*(x_f).^6);

a_f = (1/beta_f)^2.*(2*C_f(3)+6*C_f(4)*x_f+12*C_f(5)*(x_f).^2+20*C_f(6)*(x_f).^3 ...
    +30*C_f(7)*(x_f).^4+42*C_f(8)*(x_f).^5);

j_f = (6*C_f(4)+24*C_f(5)*(x_f)+60*C_f(6)*(x_f).^2+120*C_f(7)*(x_f).^3 ...
    +210*C_f(8)*(x_f).^4);

%% Calculating SVAJ for the dwell

theta_d = linspace(0,beta_d_deg,beta_d_deg);
s_d = theta_d .* 0;
v_d = theta_d .* 0;
a_d = theta_d .* 0;
j_d = theta_d .* 0;

offset_f_deg = beta_r_deg;
offset_d_deg = offset_f_deg + beta_f_deg;

%sticking matrixes together
theta = [theta_r, theta_f + offset_f_deg, theta_d + offset_d_deg];
s = [s_r, s_f, s_d];
v = [v_r, v_f, v_d];
a = [a_r, a_f, a_d];
j = [j_r, j_f, j_d];

% Calculating For the Actual SVAJ values
S = s;
V = v.*omega;
A = a.*(omega)^2;
J = j.*(omega)^3;

%% Calculating pressure angle and generating the Cam geometry
% Pressure angle

Rb = linspace(Rb,Rb+h,length(theta));
Rp = Rf+Rb;

for i = 1:length(Rp)
    phi = atan(v./(s+Rp(i)))*180/pi;
    max_phi(i) = max(phi);
    clear phi
end
k = find(max_phi<30,1);

phi = max_phi(k); %final pressure angle for the given cam

% Roller Follower Geometry

Rp = Rf + Rb(k);

% Coordinates for the Roller Follower

x_rf = cosd(theta).*(Rp+s);
y_rf = sind(theta).*(Rp+s);


%Radius of Curvature

Rb_fff = 0:.01:2;
i2 = [];
while isempty(i2)
rho_min = Rb_fff + min(s+a);
i2 = find(rho_min>0,1);
    if isempty(i2)
        Rb_fff = Rb_fff + 2;
    end
end

% Flat Faced Follower Geometry

Rb_fff = Rb_fff(i2);

x_fff = (Rb_fff+s).*sind(theta)+v.*cosd(theta);
y_fff = (Rb_fff+s).*cosd(theta)+v.*sind(theta);

end


