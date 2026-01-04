clc; clear; close all;

%% Load Data
nu_visc = 1.51 * 1e-5;
f_sample = 50000;
load("Exercise3.mat");
data = Jet(12);

%% Question 1
fluc_u = data.u - mean(data.u);  % u' = u-u_mean
u_mean =  mean(data.u);

figure(1);
plot(data.t, fluc_u);
grid on;
xlabel("Time (s)");
ylabel("u' (m/s)");
hold off;

figure(2);
plot(data.t, fluc_u);
hold on;
yline(0, 'k--');
grid on;
xlabel("Time (s)");
xlim([4.00 4.02])
ylim([-6 6])
ylabel("u' (m/s)");
hold off;


%% Question 2
%Compute statistics
u_bar = mean(data.u);
sigma_u = std(data.u);
var_u = var(data.u);
S_u = skewness(data.u);
K_u = kurtosis(data.u);
turb_intesity = sigma_u /u_bar;


%% Question 3
% Compare with normal distribution
[pdf_u,x] = ksdensity(fluc_u);
pn = 1/sigma_u/sqrt(2*pi)*exp(-fluc_u.^2/(2*var_u));

mean_fluc = mean(fluc_u);

figure(3);
plot(fluc_u, pn,'.');
hold on;
plot(x, pdf_u, '.','MarkerSize',8);
xlim([-15 15])
grid on;
xlabel("u' (m/s)");
ylabel("p(u')");
legend('Normal Distribution', 'Turbulent velocity signal', 'Mean','\sigma_u');

%% Question 4
% Time correlation
n=1; 
R_E(n) = mean(fluc_u.^2)/var_u;
while R_E(n)>-0.01
    n=n+1;
    R_E(n) = mean(fluc_u(1:end-n+1).*fluc_u(n:end))/var_u;
end

dt = data.t(2) - data.t(1);


T_E = trapz(data.t(1:n),R_E);
%These two calculations are equivalent according to Furhmann.
R_tt = 2*(R_E(2) - R_E(1))/dt^2;
tau_E2 = 1/sqrt(-0.5*R_tt);
tau_E = sqrt(2*var_u/mean((diff(fluc_u)/dt).^2));

figure(4);
subplot(2,1,1)
plot(data.t(1:n), R_E, 'LineWidth', 1.5, 'DisplayName', 'Computed');
hold on;
plot(data.t(1:10:n), exp(-data.t(1:10:n)/T_E), '.', 'MarkerSize', 10, 'DisplayName', 'e^{-\tau/T_E}');
ylim([-0.05, 1.05]);
xlim([0 data.t(n)]);

% Find where the curve crosses zero
cross_index = find(diff(sign(R_E)) ~= 0, 1); % Find the first zero-crossing point
if ~isempty(cross_index)
    tau_cross = data.t(cross_index); % Corresponding tau value
    y_cross = R_E(cross_index); % Corresponding R_E value, should be near zero
    plot([tau_cross, tau_cross], [-0.05, 1.05], '--k', 'LineWidth', 1, 'DisplayName', 'R_E(\tau)=0'); % Vertical line at crossing
end

grid on;
xlabel("\tau (s)");
ylabel("R_E(\tau)");
legend('Location', 'best');

subplot(2,1,2)
plot(data.t(1:n), R_E, 'LineWidth', 1.5, 'DisplayName', 'R_E(\tau)');
hold on;
plot(data.t(1:n), 1 - data.t(1:n).^2 / tau_E^2,'r', 'LineWidth', 1.5, 'DisplayName', '1-\tau^2/\tau_E^2');
plot([0 T_E T_E], [1 1 0], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Transition Line');
ylim([0, 1.05]);
xlim([0 3*T_E]);
grid on;
xlabel("\tau (s)");
ylabel("R_E(\tau)");
text(tau_E + 0.0001, 0.07, '\tau_E');
text(T_E + 0.0001, 0.07, 'T_E');
%legend('Location', 'best');


figure(9);
plot(data.t, fluc_u);
hold on;
yline(0, 'k--');
x_fill = [4.00742-tau_E/2, 4.00742+tau_E/2, 4.00742+tau_E/2, 4.00742-tau_E/2];
y_fill = [-10, -10, 8, 8]; 
fill(x_fill, y_fill, 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
grid on;
x_fill = [4.01200-T_E/2, 4.012+T_E/2, 4.012+T_E/2, 4.012-T_E/2];
y_fill = [-10, -10, 8, 8];
fill(x_fill, y_fill, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
xlabel("Time (s)");
xlim([4.00 4.02])
ylim([-6 8])
ylabel("u' (m/s)");
legend('Velocity fluctuations','','\tau_E Micro scale example','T_E Macro scale example','Location','best')
hold off;

%% Question 5
% Taylor macro length scale
Lambda_f = u_bar * T_E;
lambda_f = u_bar * tau_E;

Re = lambda_f*u_mean/nu_visc;

%% Question 6
% Kolmogorov scale
epsilon = 30 * nu_visc * var_u / (lambda_f^2);
tau_K = sqrt(nu_visc / epsilon);
eta_K = (nu_visc^3/epsilon)^(1/4);

%% Question 7
N = length(fluc_u);         
f_sample = 1/dt;           
f_N = f_sample / 2;         
df = f_sample / N;          

%% FFT Analysis
fft_values = fft(fluc_u);    
Pxx_fft = (1/(f_sample * N)) * abs(fft_values).^2; % Normalize
Pxx_fft(2:end-1) = 2 * Pxx_fft(2:end-1); 
fft_freq = (0:(N/2)) * df; 

% Welch's Method
% Use pwelch to calculate PSD
[welch_values, welch_freq] = pwelch(fluc_u, [], [], [], f_sample);

integral_signal = trapz(fft_freq(1:(N/2)),Pxx_fft(1:(N/2)));

integral_welch = trapz(welch_freq, welch_values);

figure(5);
loglog(fft_freq, Pxx_fft(1:N/2+1), '.'); % Plot FFT power spectrum
hold on;
grid on;
loglog(welch_freq, welch_values, '.');   % Plot Welch's PSD
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density S(f) (m^2/s)'); %I am not sure whether this should power or energy we can ask about it
legend('Computed - FFT', 'Computed - pwelch','Location','best');
hold off;
xlim([0 25000])


%% Question 8 - 10
% Wavenumber spectrum 
S_f = Pxx_fft(1:(N/2));
F_k = u_bar /(4*pi) .*S_f;

k = 2*pi*fft_freq(1:(N/2)) / u_bar;
integral_f = 2*trapz(k, F_k); % Multiplying by 2 = variance

%Welch
S_f_Welch = welch_values;
F_k_Welch = u_bar /(4*pi) .*S_f_Welch;

k_Welch = 2*pi*welch_freq / u_bar;
integral_f_Welch = 2*trapz(k_Welch, F_k_Welch);

%von Karman
F_k_karman = (Lambda_f * var_u / pi) ./ (1 + 70.78 .* (k .* Lambda_f / (2 * pi)).^2).^(5/6);

% Plot wave number spectrum
figure(6);
loglog(k, F_k, '.', 'MarkerSize', 8);
hold on;
grid on;
loglog(k_Welch, F_k_Welch, '.', 'MarkerSize', 8);
loglog(k, F_k_karman, 'LineWidth', 1.5);

% Find index where von Karman wavenumber is 10^4
[~, end_idx] = min(abs(k - 1e4)); 
k_end = k(end_idx); 
F_k_end = F_k_karman(end_idx); 

% Define the inertial subrange
k_min = 1e1;
inertial_idx = k >= k_min & k <= k_end; 
k_inertial = k(inertial_idx); 

% Compute the -5/3 slope line
F_k_ref = F_k_end * (k_inertial / k_end).^(-5/3); 
loglog(k_inertial, F_k_ref, 'k-.', 'LineWidth', 1.5); 

xline(1e2, 'r--', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(4500, 'r--', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');

% Characteristic wave numbers - question 10
k_macro = 1 / Lambda_f;
k_micro = 1 / lambda_f;
k_kolmogorov = 1 / eta_K;

xline(k_macro, 'm--', 'k_{macro} \approx \Lambda_f^{-1}', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(k_micro, 'm--', 'k_{micro} \approx \lambda_f^{-1}', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(k_kolmogorov, 'm--', 'k_{kolmogorov} \approx \eta_f^{-1}', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
x_fill = [1e-2, 1e2, 1e2, 1e-2]; 
y_fill = [1e-14, 1e-14, 1e0, 1e0]; 
fill(x_fill, y_fill, 'c', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); 
x_fill = [k_kolmogorov, 1e5, 1e5, k_kolmogorov]; 
y_fill = [1e-14, 1e-14, 1e0, 1e0]; 
fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); 
x_fill = [4500, k_kolmogorov, k_kolmogorov, 4500]; 
y_fill = [1e-14, 1e-14, 1e0, 1e0]; 
fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); 
x_fill = [1e2, 4500, 4500, 1e2];
y_fill = [1e-14, 1e-14, 1e0, 1e0]; 
fill(x_fill, y_fill, 'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); 
line([1e2, 1e5], [0.5, 0.5], 'Color', 'g', 'LineWidth', 1.5); 

legend('Computed - FFT', 'Computed - pwelch', 'von Karman equation', '-5/3 slope','','' ...
    ,'k_{macro} = 36.421','k_{micro} = 319.832','k_{Kolmogorov} = 1.898 \times 10^4','Energy Containing subrange','Dissipation subrange','','Inertial subrange','Universal equilibrium range','Location','southwest', 'FontSize', 14);
ylabel('Wave number spectrum F(k) (m^3/s^2)', 'FontSize', 16);
xlabel('Wavenumber k (m^{-1})', 'FontSize', 16);
set(gca, 'FontSize', 14); 
hold off;


%% Question 11
% Constants
nu_visc = 1.51e-5; 
D = 0.03;
num_tests = 12; 


Re = zeros(1, num_tests);
Lambda_f = zeros(1, num_tests);
lambda_f = zeros(1, num_tests);
eta_K = zeros(1, num_tests);
L = zeros(1, num_tests);
var_u = zeros(1, num_tests);

% Load data 
load("Exercise3.mat"); 
for i = 1:num_tests
    data = Jet(i);
    
    % Extract fluctuations
    fluc_u = data.u - mean(data.u);
    u_bar = mean(data.u); % Mean velocity
    var_u(i) = var(data.u);

    % Time parameters
    dt = data.t(2) - data.t(1);
    
    % Turbulent integral scale
    n = 1; 
    R_E(n) = mean(fluc_u.^2) / var_u(i);
    while R_E(n) > 0
        n = n + 1;
        R_E(n) = mean(fluc_u(1:end-n+1) .* fluc_u(n:end)) / var_u(i);
    end
    
    % Ensure same length
    R_E = R_E(:);
    T_E = trapz(data.t(1:length(R_E)), R_E); 
    
    % Taylor micro length scale
    tau_E = sqrt(2 * var(fluc_u) / mean((diff(fluc_u) / dt).^2));
    lambda_f(i) = u_bar * tau_E;
    
    % Taylor macro length scale
    Lambda_f(i) = u_bar * T_E;
    
    % Kolmogorov scales
    epsilon = 30 * nu_visc * var(fluc_u) / (lambda_f(i)^2); % Energy dissipation rate
    eta_K(i) = (nu_visc^3 / epsilon)^(1/4); % Kolmogorov length scale
    
    % Reynolds number
    Re(i) = u_bar * D / nu_visc; % Reynolds number

    L(i) = var_u(i)^(3/2) / epsilon;
end

% Compute scale ratios
ratio_Lambda_eta = Lambda_f ./ eta_K; 
ratio_Lambda_lambda = Lambda_f ./ lambda_f;

figure;
loglog(Re, ratio_Lambda_eta, 'o-', 'LineWidth', 1.5, 'DisplayName', '\Lambda_f / \eta_K');
hold on;

C_1_eta = ratio_Lambda_eta(1); 
alpha_eta = 3/4; 

Re_trend = logspace(log10(min(Re)), log10(max(Re)), 100); 
trend_eta = C_1_eta * (Re_trend / Re(1)).^alpha_eta; 

loglog(Re_trend, trend_eta, '--', 'LineWidth', 1.5, 'DisplayName', ['Trendline ~' sprintf('%.2f', C_1_eta) '  Re^{3/4}']);

grid on;
xlabel('Reynolds number (Re)');
ylabel('Length scale ratio \Lambda_f / \eta_K');
legend('show', 'Location', 'northwest');
hold off;

figure;
loglog(Re, ratio_Lambda_lambda, 'o-', 'LineWidth', 1.5, 'DisplayName', '\Lambda_f / \lambda_f');
hold on;

C_2_lambda = ratio_Lambda_lambda(1); 
beta_lambda = 1/2; 

trend_lambda = C_2_lambda * (Re_trend / Re(1)).^beta_lambda; 
loglog(Re_trend, trend_lambda, '--', 'LineWidth', 1.5, 'DisplayName', ['Trendline ~' sprintf('%.2f', C_2_lambda) '  Re^{1/2}']);


grid on;
xlabel('Reynolds number (Re)');
ylabel('Length scale ratio \Lambda_f / \lambda_f');
legend('show', 'Location', 'northwest');
hold off;