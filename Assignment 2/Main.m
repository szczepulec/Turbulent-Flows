clc; clear; close all;

%% Info
Cases = getAllWaveConditions();
nu = 1.14e-6;

%Cases
% 1 - Laminar 
% 2 - Transition 
% 3 - Turbulent 
% 4 - Turbulent (rough)

%% Question 1
for i = 1:4
    % Rough turbulent Case 4
    if Cases(i).ks > 0
        f_w = exp(5.5 * (Cases(i).a_over_ks^(-0.16)) - 6.7); % eqn 5.69
    % Laminar flow
    elseif Cases(i).Re <= 1.5e5
        f_w = 2 / sqrt(Cases(i).Re); % eqn 5.59
    % Turbulent flow
    elseif Cases(i).Re >= 5e5
        f_w = 0.035 / (Cases(i).Re^(0.16)); % eqn 5.60
    end

    %% maximum velocity
    U_fm = sqrt(f_w / 2) * Cases(i).U0m; % eqn 5.57

    % update struct
    Cases(i).f_w = f_w;
    Cases(i).U_fm = U_fm;
end

%% Question 2
for i = 1:4
    if Cases(i).caseNumber == 4
        delta_y = Cases(i).ks*0.02; % eqn 9.54
        ks = Cases(i).ks;
        ks_plus = ks*Cases(i).U_fm/nu; % eqn 9.24
    else
        delta_y = nu/Cases(i).U_fm; % eqn 9.49
        ks = 0.1* nu/Cases(i).U_fm; % eqn 9.24
        ks_plus = ks*Cases(i).U_fm/nu; % eqn 9.24
    end

    % update struct
    Cases(i).delta_y = delta_y;
    Cases(i).ks = ks;
    Cases(i).ks_plus = ks_plus;
end

CaseInfo = Cases;
save('CaseInfo.mat', 'CaseInfo');

%% Case 1 %%
load('CaseInfo.mat')
T = CaseInfo(1).T;
nu = 1.14e-6;
rho = 1000;
mu = rho*nu;
U0m =  CaseInfo(1).U0m;
hm = 0.145;

%% Laminar theory solution 
t = 0:0.1:T; 
y = 0:0.0001:hm;
wt = [0 45 90 135]; 
figure; set(gcf, 'Position', [100, 100, 1200, 400]); hold on;

%% MatRANS
load('out_MatRANS_case1.mat') 

num_phases = length(wt);
for i = 1:num_phases
    omega = 2*pi/T;
    a = U0m/omega; % eqn 5.16
    delta_1 = sqrt(2*nu/omega); % eqn 5.13
    delta = 3*pi/4* delta_1;

    u_theory = U0m * sin(deg2rad(wt(i))) - U0m * exp(-y / delta_1) .* sin(deg2rad(wt(i)) - y / delta_1); % eqn 5.12
    
    time_RANS = deg2rad(wt(i)) / omega;
    index = find(MatRANS.t == time_RANS);

    u_norm_RANS = MatRANS.u(index, :) ./ U0m;
    y_norm_RANS = MatRANS.y ./ a;

    subplot(1, 4, i);
    plot(u_norm_RANS, y_norm_RANS, 'DisplayName', 'RANS k-$\omega$', 'LineWidth', 1.5); hold on;
    plot(-30, -30, 'DisplayName', '', 'HandleVisibility', 'off');
    plot(u_theory / U0m, y / a, 'DisplayName', 'Laminar Theory', 'LineWidth', 1.5);

    xlabel("$\overline{U} / U_{0m}$", 'Interpreter', 'latex', 'FontSize', 12);
    ylabel("$y / a $", 'Interpreter', 'latex', 'FontSize', 12);
    title(['$\omega t = ' num2str(wt(i)) '^\circ$'], 'FontSize', 12, 'Interpreter', 'latex');
    ylim([0 0.04]); xlim([min(u_theory / U0m) max(u_theory / U0m)+0.2]); grid on;
end
hold off; legend('show', 'Interpreter', 'latex', 'Location', 'north');

wt = -30:1:360; 
tau_0 = mu * U0m / delta_1 * (sin(deg2rad(wt)) + cos(deg2rad(wt)));

wt_RANS_rad = MatRANS.t*omega;
wt_RANS_deg = rad2deg(wt_RANS_rad);
tau_0_norm_RANS = MatRANS.tau0 / (rho * U0m^2);

figure(2);
plot(wt_RANS_deg, tau_0_norm_RANS, 'DisplayName','RANS k-$\omega$', 'LineWidth', 1.5); hold on;
plot(-30, -30, 'DisplayName', '', 'HandleVisibility', 'off');
plot(wt, tau_0 / (rho * U0m^2), 'DisplayName', 'Laminar Theory', 'LineWidth', 1.5);
xlim([0 360]); ylim([-4e-3 4e-3]);
ylabel("$\tau_0 / \rho U_{0m}^2$ [m]", 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$\omega t\ [^\circ]$', 'Interpreter', 'latex', 'FontSize', 12);
legend('show', 'Interpreter', 'latex'); grid on; hold off;

%% Case 2 %%
load('CaseInfo.mat')
T = CaseInfo(2).T;
nu = 1.14e-6; rho = 1000; mu = rho*nu; U0m =  CaseInfo(2).U0m; hm = 0.145;

omega = 2*pi/T; a = U0m/omega; delta_1 = sqrt(2*nu/omega); delta = 3*pi/4* delta_1;

%% Theory
wt_theory = -30:1:360; 
tau_0_Theory = mu * U0m / delta_1 * (sin(deg2rad(wt_theory)) + cos(deg2rad(wt_theory)));
tau_0_norm_Theory = tau_0_Theory / (rho * U0m^2);

%% RANS
load('out_MatRANS_case2.mat') 
wt_RANS_deg = rad2deg(MatRANS.t*omega);
indices_one_period = wt_RANS_deg >= 360*4 & wt_RANS_deg < 360*5;
wt_RANS_deg_one_period = wt_RANS_deg(indices_one_period)-(360*4);

tau_0_RANs = MatRANS.tau0;
tau_0_RANs_5phase = tau_0_RANs(wt_RANS_deg > 360*4 & wt_RANS_deg < 5*360);
tau_0_norm_RANS_5phase = tau_0_RANs_5phase / (rho * U0m^2);

%% Experiment
load('Exercise4.mat')
c = 2;
wt_EXP = WBL(c).omegat_tau0;
tau_0_exp = WBL(c).tau0;
tau_0_norm_EXP = tau_0_exp / (rho * U0m^2);

%% Plot
figure(3);
plot(wt_RANS_deg_one_period, tau_0_norm_RANS_5phase, 'DisplayName', 'RANS k-$\omega$', 'LineWidth', 1.5); hold on;
plot(wt_EXP, tau_0_norm_EXP, '.', 'DisplayName', 'Experimental', 'MarkerSize', 18);
plot(wt_theory, tau_0_norm_Theory, 'DisplayName', 'Laminar Theory', 'LineWidth', 1.5);
xlim([0 360]);
ylabel("$\tau_0 / \rho U_{0m}^2$ [m]", 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$\omega t\ [^\circ]$', 'Interpreter', 'latex', 'FontSize', 12);
legend('show','Interpreter','latex','FontSize',12); grid on; hold off;

%% Theory Cf
Cf_norm_theory = (2*tau_0_Theory/rho) ./ (U0m^2 * sin(deg2rad(wt_theory)+pi/4));
Cf_theory_laminar = 2/(CaseInfo(2).Re)^(0.5);
fprintf('Cf Theory Laminar = %.4f \n', Cf_theory_laminar)
Cf_theory_turbulent = 0.035/(CaseInfo(2).Re)^(0.16);
fprintf('Cf Theory turbulent = %.4f \n', Cf_theory_turbulent)

%% RANS Cf
Cf_norm_RANS = (2*tau_0_RANs_5phase/rho) ./ (U0m^2 * sin(deg2rad(wt_RANS_deg_one_period)+pi/4));
Cf_RANS = (2*max(tau_0_RANs_5phase)/rho) ./ (U0m^2);
fprintf('Cf RANS = %.4f \n', Cf_RANS)

%% Experiment Cf
Cf_norm_EXP = (2*tau_0_exp/rho) ./ (U0m^2 * sin(deg2rad(wt_EXP)+pi/4));

figure(4);
plot(wt_RANS_deg_one_period, Cf_norm_RANS, 'DisplayName', 'RANS k-$\omega$', 'LineWidth', 1.5); hold on;
plot(wt_EXP, Cf_norm_EXP, '.', 'DisplayName', 'Experimental', 'MarkerSize', 18);
plot(wt_theory, Cf_norm_theory, 'DisplayName', 'Laminar Theory', 'LineWidth', 1.5);
xlim([0 130]);
ylabel("$f_w^*$", 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$\omega t\ [^\circ]$', 'Interpreter', 'latex', 'FontSize', 12);
legend('show','Interpreter','latex','FontSize',12); grid on; hold off;

%% Cases 3 & 4 
cases_to_plot = [3 4];
n = [1 4 7 10];
phase_angle_deg = [0 45 90 135];

for c = cases_to_plot
    load(sprintf("out_MatRANS_case%d.mat", c));
    load("Exercise4.mat");

    U0m = MatRANS.U0m;
    T = CaseInfo(c).T;
    omega = 2*pi/T;
    a = U0m / omega;

    % phase times
    phase_angle_rad = deg2rad(phase_angle_deg) + 8*pi; % fifth period
    time_list = phase_angle_rad ./ omega;
    tolerance = 1e-5;
    index_list = zeros(size(time_list));

    for i = 1:length(time_list)
        idx = find(abs(MatRANS.t - time_list(i)) < tolerance, 1);
        index_list(i) = ~isempty(idx) * idx;
    end

    % set limits for gprahs
    if c == 3
        u_xlim = [-0.02 1.1]; u_ylim = [0 0.05];
        k_xlim = [0 0.007]; uv_xlim = [-1e-3 1.7e-3];
    elseif c == 4
        u_xlim = [-0.02 1.1]; u_ylim = [0 0.05];
        k_xlim = [0 0.011]; uv_xlim = [-1e-3 4e-3];
    end

    %% velocity plots
    figure(10*c-27); set(gcf, 'Position', [100,100,1200,400]);
    for i = 1:length(phase_angle_rad)
        idx = index_list(i);
        u_plot = MatRANS.u(idx,:) / U0m;
        y_plot = MatRANS.y / a;
        u_exp = WBL(c).u(:, n(i)) / WBL(c).U0m;
        y_exp = WBL(c).y_u / a;

        subplot(1,length(phase_angle_rad),i);
        plot(u_plot, y_plot, 'LineWidth',1); hold on;
        plot(u_exp, y_exp, 'o', 'MarkerFaceColor',[1 0.5 0], 'MarkerSize',3);
        xlabel("$\overline{U} / U_{0m}$",'Interpreter','latex','FontSize',12);
        ylabel("$y / a$",'Interpreter','latex','FontSize',12);
        title(sprintf('$\\omega t = %d^\\circ$', phase_angle_deg(i)),'Interpreter','latex','FontSize',12);
        xlim(u_xlim); ylim(u_ylim); grid on;
        if i==length(phase_angle_rad), legend('RANS k-$\omega$','Experimental','Interpreter','latex','FontSize',12,'Location','best'); end
    end

    %% k / U0m^2 plots
    figure(10*c-26); set(gcf, 'Position', [100,100,1200,400]);
    for i = 1:length(phase_angle_rad)
        idx = index_list(i);
        k_plot = MatRANS.k(idx,:) / (U0m^2);
        k_exp = 0.65*(WBL(c).uu(:, n(i)) + WBL(c).vv(:, n(i))) / (WBL(c).U0m^2);
        y_k_exp = WBL(c).y_uuvv / a;

        subplot(1,length(phase_angle_rad),i);
        plot(k_plot, y_plot, 'LineWidth',1); hold on;
        plot(k_exp, y_k_exp, 'o','MarkerFaceColor',[1 0.5 0],'MarkerSize',3);
        xlabel("$k / U_{0m}^2$",'Interpreter','latex','FontSize',12);
        ylabel("$y / a$",'Interpreter','latex','FontSize',12);
        title(sprintf('$\\omega t = %d^\\circ$', phase_angle_deg(i)),'Interpreter','latex','FontSize',12);
        xlim(k_xlim); grid on;
        if i==length(phase_angle_rad), legend('RANS k-$\omega$','Experimental','Interpreter','latex','FontSize',12,'Location','best'); end
    end

    %% uv plots
    figure(10*c-25); set(gcf, 'Position', [100,100,1200,400]);
    for i = 1:length(phase_angle_rad)
        idx = index_list(i);
        uv_plot = MatRANS.nu_t(idx,:) .* gradient(MatRANS.u(idx,:), MatRANS.y) / (U0m^2);
        uv_exp = -WBL(c).uv(:, n(i)) / (WBL(c).U0m^2);
        y_uv_exp = WBL(c).y_uv / a;

        subplot(1,length(phase_angle_rad),i);
        plot(uv_plot, y_plot, 'LineWidth',1); hold on;
        plot(uv_exp, y_uv_exp,'o','MarkerFaceColor',[1 0.5 0],'MarkerSize',3);
        xlabel("$-\overline{u'v'} / U_{0m}^2$",'Interpreter','latex','FontSize',12);
        ylabel("$y / a$",'Interpreter','latex','FontSize',12);
        title(sprintf('$\\omega t = %d^\\circ$', phase_angle_deg(i)),'Interpreter','latex','FontSize',12);
        xlim(uv_xlim); grid on;
        if i==length(phase_angle_rad), legend('RANS k-$\omega$','Experimental','Interpreter','latex','FontSize',12,'Location','best'); end
    end

    %% tau0 comparison
    time_tau = MatRANS.t(index_list(1):end);
    omegat = (omega*time_tau - 8*pi)*180/pi;
    tau0_plot = MatRANS.tau0(index_list(1):end) / (MatRANS.rho * U0m^2);
    tau0_exp = WBL(c).tau0 / (MatRANS.rho * U0m^2);
    omegat_exp = WBL(c).omegat_tau0;

    delta_1 = sqrt(2*nu/omega);
    wt_full = -30:1:360;
    tau0_theory = mu*U0m/delta_1*(sin(deg2rad(wt_full)) + cos(deg2rad(wt_full)));

    figure(10*c-24);
    plot(omegat, tau0_plot, 'LineWidth',1.75); hold on;
    plot(omegat_exp, tau0_exp,'o','MarkerFaceColor',[1 0.5 0],'MarkerSize',3);
    plot(wt_full, tau0_theory/(rho*U0m^2),'LineWidth',1.5);
    xlabel('$\omega t\ [^\circ]$','Interpreter','latex','FontSize',12);
    ylabel("$\tau_0 / \rho U_{0m}^2$ [m]",'Interpreter','latex','FontSize',12);
    legend('RANS k-$\omega$','Experimental','Laminar Theory','Interpreter','latex','FontSize',12);
    xlim([0 360]); grid on; hold off;
end

%% Comparison case 3 and 4
phase = 90;
phase_rad = deg2rad(phase) + 8*pi;
time = phase_rad / omega;
idx = 307;

load("out_MatRANS3.mat"); u3 = MatRANS.u(idx,:)/U0m; k3 = MatRANS.k(idx,:)/U0m^2; y3 = MatRANS.y / a;
load("out_MatRANS4.mat"); u4 = MatRANS.u(idx,:)/U0m; k4 = MatRANS.k(idx,:)/U0m^2; y4 = MatRANS.y / a;

figure(99);
subplot(1,2,1);
plot(u3, y3,'b','LineWidth',1); hold on; grid on;
plot(u4, y4,'r','LineWidth',1);
legend('Case 3: Smooth bed','Case 4: Rough bed','Location','best');
xlabel("$\overline{U} / U_{0m}$",'Interpreter','latex','FontSize',14);
ylabel("$y / a$",'Interpreter','latex','FontSize',14); xlim([0 1.1]);

subplot(1,2,2);
plot(k3, y3,'b','LineWidth',1); hold on; grid on;
plot(k4, y4,'r','LineWidth',1);
xlabel("$k / U_{0m}^2$",'Interpreter','latex','FontSize',14);
ylabel("$y / a$",'Interpreter','latex','FontSize',14);
