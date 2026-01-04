clc; clear; close all;

%% Load data
load('Exercise1.mat')

h = Channel(1).h; 
nu= Channel(1).nu; 

base = 0.30; % base length [m]
U_0 = 0.30;

[u_mean, v_mean, uu_rms, vv_rms, uv_rms, y] = compute_velocities_rms(Channel, U_0, base);

%{
plot(u_mean, y, 'LineWidth', 2);
grid on;
ylabel(' y [m]')
xlabel(' u̅ [m/s]')
%}

V = (1/h) * trapz(y, u_mean); %eqn: 10.3

fprintf('The depth averaged velocity V = %.4f\n', V);



Area = h*base; % area
P = 2*h+base; %perimiter
r_h = Area/P; %given in question

Re = r_h * V / nu; %Reynolds

f = 0.0557/(Re^(0.25)); % eqn 10.5

U_f_estimate = sqrt(f/2)*V; % eqn 10.4
 
fprintf('Estimated friction velocity U_f = %.4f\n', U_f_estimate);


% Part 1
y_plus_estimate = U_f_estimate*y/nu;
Re_tau_estimate = h*U_f_estimate/nu;

%Logarimic layer bounds
Upper_bound_y_plus = 0.1*Re_tau_estimate;
Lower_bound_y_plus = 30;

Upper_bound_y = Upper_bound_y_plus*nu/U_f_estimate;
Lower_bound_y = Lower_bound_y_plus*nu/U_f_estimate;

% data in bounds
valid_indices = y >= Lower_bound_y & y <= Upper_bound_y;

% Log of y data and normal u
log_y = log(y(valid_indices));          
u_mean_data = u_mean(valid_indices);   

% fit
p = polyfit(log_y, u_mean_data, 1);   
y_fit = linspace(Lower_bound_y, Upper_bound_y, 100);
log_y_fit = log(y_fit);                
u_fit = polyval(p, log_y_fit); 



figure(1);
semilogy(u_mean, y, 'LineWidth', 2); 
hold on;  
% plot bounds y
yline(Upper_bound_y, 'k--', 'LineWidth', 1); 
yline(Lower_bound_y, 'k--', 'LineWidth', 1);  
%plot fit
semilogy(u_fit, y_fit, 'r-', 'LineWidth', 1.5);  % Fitted line in red


% text
A = p(1);
B = p(2);  

equation_str = sprintf(' u̅  = %.3f * ln(y) + %.3f', A, B);

text(0.1, max(y) * 0.9, equation_str, ...
    'FontSize', 12, 'Color', 'black', 'BackgroundColor', 'white', 'EdgeColor', 'black');
grid on;
ylabel('y [m]')
xlabel('u̅ [m/s]')
hold off; 

% finding U_f
y_specific = mean(y(valid_indices));  
%U_f_calc = (mean(u_mean(valid_indices)) - B) / (A * log(y_specific)); % eqn 3.42 u_mean = A*U_f *ln(y)+B
U_f_calc = 0.033/2.5;
fprintf('Calculated U_f = %.4f\n', U_f_calc);



% Plot 1

y_plus_calc = U_f_calc*y/nu;
U_plus_calc = u_mean/U_f_calc;


%region bouds
%y_plus_bounds = [0, 5, 5, 30, 30, Upper_bound_y_plus, Upper_bound_y_plus, 1000];
%U_plus_bounds = [0, 25, 25, 0, 0, 25, 25, 0]; % Fixed bounds for U^+

figure(2);  
semilogy(U_plus_calc, y_plus_calc, 'LineWidth', 2);
hold on;

% Viscous sublayer
fill([0, 25, 25, 0], ...  % x-coordinates for viscous sublayer
     [1, 1, 5, 5],  [0.4, 0.8, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Buffer layer
fill([0, 25, 25, 0], ...  % x-coordinates for buffer layer
     [5, 5, 30, 30], [1, 0.6, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Logarithmic layer
fill([0, 25, 25, 0], ...  % x-coordinates for logarithmic layer
     [30, 30, Upper_bound_y_plus, Upper_bound_y_plus], [0.6, 1, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Outer layer
fill([0, 25, 25, 0], ...  % x-coordinates for outer layer
     [Upper_bound_y_plus, Upper_bound_y_plus, 1000, 1000], [0.0, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Plot the data
grid on;
ylim([0 1001]);  % Set y-axis limits from 0 to 100
xlim([0 25]);   % Set x-axis limits from 0 to 25
ylabel('y^+');
xlabel('u̅/U_f');

% Add legend if necessary
legend({'Data','Viscous Sublayer', 'Buffer Layer', 'Logarithmic Layer', 'Outer Layer'}, 'Location', 'Best');
hold off;

figure(3); 
grid on;
semilogy(U_plus_calc, y / h, 'LineWidth', 2);
hold on;

% Viscous sublayer
fill([0, 25, 25, 0], ...  % x-coordinates for viscous sublayer
     [1*nu/(U_f_calc * h), 1*nu/(U_f_calc * h), 5*nu/(U_f_calc * h), 5*nu/(U_f_calc * h)], [0.4, 0.8, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Buffer layer
fill([0, 25, 25, 0], ...  % x-coordinates for buffer layer
     [5*nu/(U_f_calc * h), 5*nu/(U_f_calc * h), 30*nu/(U_f_calc * h), 30*nu/(U_f_calc * h)], [1, 0.6, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Logarithmic layer
fill([0, 25, 25, 0], ...  % x-coordinates for logarithmic layer
     [30*nu/(U_f_calc * h), 30*nu/(U_f_calc * h), Upper_bound_y_plus*nu/(U_f_calc * h), Upper_bound_y_plus*nu/(U_f_calc * h)], [0.6, 1, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Outer layer
fill([0, 25, 25, 0], ...  % x-coordinates for outer layer
     [Upper_bound_y_plus*nu/(U_f_calc * h), Upper_bound_y_plus*nu/(U_f_calc * h), 1000*nu/(U_f_calc * h), 1000*nu/(U_f_calc * h)],[0.0, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Plot the data
grid on;
ylim([0 1000*nu/(U_f_calc * h)]);  % Set y-axis limits from 0 to 100
xlim([0 25]);   % Set x-axis limits from 0 to 25
ylabel('y/h');
xlabel('u̅/U_f');

% Add legend if necessary
legend({'Data','Viscous Sublayer', 'Buffer Layer', 'Logarithmic Layer', 'Outer Layer'}, 'Location', 'Best');
hold off;

%eqn 3.108
kappa = 0.4;
A_d = 25;
u_mean_van_driest = 2* U_f_calc * cumtrapz(y_plus_calc, 1 ./ ( 1 + ( 1 + 4 .* kappa^2 .* y_plus_calc.^2 .* ( 1 -exp(- y_plus_calc/A_d)).^2) .^(0.5)));
u_mean_van_driest_Plus = u_mean_van_driest/U_f_calc;


figure(4);  
semilogy(U_plus_calc, y_plus_calc, 'LineWidth', 2);
hold on;
semilogy(u_mean_van_driest_Plus, y_plus_calc, 'LineWidth', 2);

% Viscous sublayer
fill([0, 25, 25, 0], ...  % x-coordinates for viscous sublayer
     [1, 1, 5, 5], [0.4, 0.8, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Buffer layer
fill([0, 25, 25, 0], ...  % x-coordinates for buffer layer
     [5, 5, 30, 30], [1, 0.6, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Logarithmic layer
fill([0, 25, 25, 0], ...  % x-coordinates for logarithmic layer
     [30, 30, Upper_bound_y_plus, Upper_bound_y_plus], [0.6, 1, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Outer layer
fill([0, 25, 25, 0], ...  % x-coordinates for outer layer
     [Upper_bound_y_plus, Upper_bound_y_plus, 1000, 1000], [0.0, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  

% Plot the data
grid on;
ylim([0 1001]);  % Set y-axis limits from 0 to 100
xlim([0 25]);   % Set x-axis limits from 0 to 25
ylabel('y^+');
xlabel('u̅/U_f');

% Add legend if necessary
legend({'Data','Van Driest velocity distribution','Viscous Sublayer', 'Buffer Layer', 'Logarithmic Layer', 'Outer Layer'}, 'Location', 'Best');
hold off;

% comparison with relative eror
relative_error = abs((U_plus_calc - u_mean_van_driest_Plus) ./ U_plus_calc);

figure(5);
semilogy(relative_error, y_plus_calc, 'k-', 'LineWidth', 2);  % Error plot
title('Relative Error Between Data and Van Driest Distribution');
xlabel('Relative Error');
ylabel('y^+');
grid on;

% for y^+ 0 - 100

% Filter for 0 < y+ < 100
idx = y_plus_calc > 0 & y_plus_calc < 100;
y_plus_valid = y_plus_calc(idx);
uu_rms_valid = uu_rms(idx);
vv_rms_valid = vv_rms(idx);
uv_rms_valid = uv_rms(idx);

uu_rms_valid = [0; uu_rms_valid];
y_plus_valid = [0; y_plus_valid];
vv_rms_valid = [0; vv_rms_valid];
uv_rms_valid = [0; uv_rms_valid];


figure(6);
plot(y_plus_valid, uu_rms_valid / U_f_calc, 'LineWidth', 2);
%ylabel('y^+');
%xlabel('$\sqrt{\overline{u^{\prime 2}}}/U_f$','Interpreter','latex');
grid on;
%{xlim([0 100]); 
%ylim([0 3]); %}

hold on;
plot(y_plus_valid, vv_rms_valid / U_f_calc,'LineWidth', 2);
%ylabel('y^+');
%xlabel('$\sqrt{\overline{v^{\prime 2}}}/U_f$','Interpreter','latex');
grid on;

 

hold on;
plot(y_plus_valid, uv_rms_valid / U_f_calc, 'LineWidth', 2);
xlabel('y^+');
%ylabel('RMS of fluctuating velocity/U_f');
grid on;

legend("$\sqrt{\overline{u^{\prime 2}}}/U_f$", "$\sqrt{\overline{v^{\prime 2}}}/U_f$", "$\sqrt{-\overline{u^{\prime}v^{\prime}}}/U_f$", "Interpreter", "latex")
%xlim([0 100]); 
% maybe combine? 

y_plus_calc= [0; y_plus_calc];
uu_rms = [0; uu_rms];
vv_rms = [0; vv_rms];
uv_rms = [0; uv_rms];

y(end) = [];

figure(9); 
plot(y / h, uu_rms / U_f_calc, 'LineWidth', 2);

grid on;

hold on;
plot(y / h, vv_rms / U_f_calc, 'LineWidth', 2);

grid on;

hold on;
plot(y / h,uv_rms / U_f_calc, 'LineWidth', 2);
xlabel('y/h');
grid on;
xlim([0 1]);
%ylabel('RMS of fluctuating velocity/U_f');

legend("$\sqrt{\overline{u^{\prime 2}}}/U_f$", "$\sqrt{\overline{v^{\prime 2}}}/U_f$", "$\sqrt{-\overline{u^{\prime}v^{\prime}}}/U_f$", "Interpreter", "latex")



k = 0.5 * (uu_rms.^2 + vv_rms.^2 + 1.8.*(vv_rms).^2);


figure(12);
plot(y / h, k / U_f_calc^2, 'LineWidth', 2);
xlabel('y/h');
ylabel('$k/U_f^2$','Interpreter','latex');
grid on;

 %%%% ASK NOT GIVEN %%%% !!!!
 % maybe can be found from dynamic viscocity? mu/rho?
 % Assumed it is water as stated for now 
rho = 997;
mu = nu*rho;

Tau = rho * U_f_calc^2 *(1-y/h);

u_mean(end) = [];



Reynolds_stress = Tau - mu* gradient(u_mean, y);
uv = rho*uv_rms.^2;

figure(15);
plot(y/h, Reynolds_stress / U_f_calc^2,'LineWidth', 2);
hold on;
plot(y/h, uv/U_f_calc^2,'x','MarkerSize', 8,'LineWidth', 2);
xlabel('y/h');
ylabel('$ - \rho \overline{u^{\prime}v^{\prime}} /U_f^2$','Interpreter','latex');
grid on;
legend("Theoretical", "Experimental")
ylim([0 inf])
hold off;



EE = Reynolds_stress.* gradient(u_mean, y);
Tau_plot = Tau.* gradient(u_mean, y);

uv_EE = uv.* gradient(u_mean, y);

% Main figure
figure(16);
plot(EE, y/h, 'LineWidth', 2); % Replace this with the data you want to plot on the logarithmic scale
hold on;
grid on;
ylabel('y/h'); % Set the label for the right y-axis

% Main plot settings
xlim([-5 20]);
%ylim([0 0.1]);

% Plot additional data
plot(uv_EE, y/h, 'x', 'MarkerSize', 8, 'LineWidth', 2);
plot(Tau_plot, y/h, 'LineWidth', 2);
legend('$ Theoretical: - \rho \overline{u^{\prime}v^{\prime}} \frac{\partial \overline{u}}{\partial y}$', ...
       '$ Experimental: - \rho \overline{u^{\prime}v^{\prime}} \frac{\partial \overline{u}}{\partial y}$', ...
       '$\overline{\tau} \frac{\partial \overline{u}}{\partial y}$', 'Interpreter', 'latex');

% Create a zoomed-in inset axes
axes('Position', [0.535 0.35 0.35 0.35]); % Adjust the position and size as needed
box on; % Add box around the inset
hold on;

% Plot the zoomed-in data
plot(EE, y/h, 'LineWidth', 2); % Same data for zoomed-in view
plot(uv_EE, y/h, 'x', 'MarkerSize', 8, 'LineWidth', 2);
plot(Tau_plot, y/h, 'LineWidth', 2);

% Adjust the zoomed-in view limits
xlim([0 10]); % Set x-limits for zoomed-in section
ylim([0 0.1]); % Set y-limits for zoomed-in section
grid on; % Add grid to zoomed-in section
ylabel('y/h'); % Set label for zoomed-in section
title('Zoomed In'); % Optional: Title for the zoomed-in section
hold off;




velocity_RANS = readmatrix('velocity.csv');
energy_k_RANS = readmatrix('energy_k.csv');
yparameter = readmatrix("yparameter.csv");

figure(20);
plot(k / U_f_calc^2,y / h, 'LineWidth', 2);
hold on;
plot(energy_k_RANS, yparameter, 'LineWidth', 2);
ylabel('y/h');
xlabel('$k/U_f^2$','Interpreter','latex');
legend("Experimental", "RANS");
grid on;
hold off;
xlim([0 1]);

figure(21);
plot(u_mean / U_f_calc ,y / h, 'LineWidth', 2);
hold on;
plot(velocity_RANS, yparameter, 'LineWidth', 2);
ylabel('y/h');
xlabel('$\overline{u}/U_f$','Interpreter','latex');
legend("Experimental", "RANS");
grid on;
hold off;
turb_energy_production = rho * uv_rms.^2 .* gradient(u_mean, y);



figure(22);
semilogy(turb_energy_production, y / h, 'LineWidth', 2); % Switch axes
xlabel('$  -\rho \overline{u^{\prime}v^{\prime}} \frac{\partial \overline{u}}{\partial y}$', 'Interpreter', 'latex');
ylabel('y/h');
grid on;
xlim([0 1]);
%}



