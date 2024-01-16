%% Computational Modelling of Turbulent Air Flow in the High Speed Wind Leg of Virgina Tech Stability Wind Tunnel
%% Turbulence Modelling Assignment
% Starting file

clear 
close all
clc
% Calculating $U_{ref}$, $\rho_{ref}$, and $\nu$

%Air Properties
gamma = 1.4;
R_uni = 8.314; %Universal Gas Constant
mol_mass = 0.0289645;

%Geometric Properties
y_ts = 1.8304; %y-length of the test section
z_ts = 1.8391; %z-length of the test section
y_sc = 5.4625; %y-length of the settling chamber
z_sc = 5.4686; %z-length of the settling chamber

%Input
u_ref_nu = 3.9e6; %Reference Velo-Visc Ratio
M_ref = 0.193;
Tt_ref = 300.4; %Local (Total) Stagnation Temperature
Pt_ref = 94242; %Local (Total) Stagnation Pressure
Ps_ref = 91832; %Local Static Pressure
Pgauge = Pt_ref-Ps_ref


%Process and Output
Ts_Tt = (1 + (gamma-1)/(2)*M_ref^2); %Tstatic/Ttotal
Ts = Ts_Tt*Tt_ref; %Static Temperature [K]
a_ref = sqrt(gamma*R_uni*Ts/mol_mass); %Speed of Sound [m/s]
rho_ref = Ps_ref/(R_uni*Ts/mol_mass); %Assuming there is the flow is incompressible

%Continuity-Equation
u_ts = M_ref*a_ref;
A_ts = y_ts*z_ts;
A_sc = y_sc*z_sc;
u_sc = (A_ts*u_ts)/(A_sc);
nu = u_ts/(u_ref_nu);


% Resulting Data from CFD

u_inp_cfd = [6 7.4186 7.5847 8];
u_ref_cfd = [54.479 67.3067 68.8082 72.5623];

p = polyfit(u_ref_cfd, u_inp_cfd, 1)

% Generate x values for the plot
x_plot = linspace(50, 75, 75);

% Calculate y values for the plot
y_plot = polyval(p, x_plot);

% Plot the original data and the fitted line
figure(1);
plot(x_plot, y_plot, '-'); % plot the fitted line
hold on
scatter(u_ref_cfd, u_inp_cfd, 'red','filled','o')
ylabel('Input Velocity (m/s)', 'Interpreter','latex', 'FontSize',14)
xlabel('Pitot Static Tube Velocity (m/s)', 'Interpreter','latex', 'FontSize',14)
set(gca,'TickLabelInterpreter','latex','FontSize',14)
hold off;

%Find the point
u_sc_actual = polyval(p, u_ts);
percent_error = abs(u_ts-67.3067)
%%
disp('AIR PROPERTIES AND INPUT')
fprintf('Inlet Velocity (m/s): %.4f ', u_sc_actual);
fprintf('Initial Calculation Input Velocity (m/s): %.4f ', u_sc);
fprintf('Difference (percent): %.4f ', percent_error);
fprintf('Pitot Static Velocity (m/s): %.4f ', u_ts);
fprintf('Reference density (kg/m3): %.4f ', rho_ref);
fprintf('Kinematic viscosity (m2/s): %.4e\n', nu);
%%
%Numerical Results
n_mesh = [1423729 2814664 5613201];
m_imbal = [-0.01039 3.2924872e-05 1e-6];
v_pitot = [67.371689 67.224817 67.472389];
ps_ref_cfd = [91956.266 91828.233 91927.336];
pt_ref_cfd = [94358.844 94273.18 94282.938];

%Analytical Results
n_mesh_analytical = [0 6000000]
v_ts_array = repmat(u_ts, 1, size(n_mesh_analytical, 2))
ps_array = repmat(Ps_ref, 1, size(n_mesh_analytical, 2))
pt_array = repmat(Pt_ref, 1, size(n_mesh_analytical, 2))
%%
figure(2);
ax1 = subplot(2,2,1);
plot(n_mesh, m_imbal, '-o');
ylim([-1, 1])
ylabel('Mass Imbalance (kg/s)','Interpreter', 'latex', 'FontSize', 14);
xlabel('Number of Mesh','Interpreter', 'latex', 'FontSize', 14);
set(gca,'TickLabelInterpreter','latex','FontSize',14)

ax2 = subplot(2,2,2);
plot(n_mesh, v_pitot/, '-o');
hold on
plot(n_mesh_analytical, v_ts_array, '--', 'Color','r')
hold off
ylim([0, 80])
ylabel('Velocity (m/s)','Interpreter', 'latex', 'FontSize', 14);
xlabel('Number of Mesh','Interpreter', 'latex', 'FontSize', 14);
legend('Numerical', 'Analytical','Interpreter', 'latex', 'FontSize', 14, 'Location', 'east');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

ax3 = subplot(2,2,3);
plot(n_mesh, ps_ref_cfd, '-o');
hold on
plot(n_mesh_analytical, ps_array, '--', 'Color','r')
hold off
ylim([0, 120000])
ylabel('Static Pressure (Pa)','Interpreter', 'latex', 'FontSize', 14);
xlabel('Number of Mesh','Interpreter', 'latex', 'FontSize', 14);
legend('Numerical', 'Analytical','Interpreter', 'latex', 'FontSize', 14, 'Location', 'east');
set(gca,'TickLabelInterpreter','latex','FontSize',14)

ax4 = subplot(2,2,4);
plot(n_mesh, pt_ref_cfd, '-o');
hold on
plot(n_mesh_analytical, pt_array, '--', 'Color','r')
hold off
ylim([0,120000])
ylabel('Total Pressure (Pa)','Interpreter', 'latex', 'FontSize', 14);
xlabel('Number of Mesh','Interpreter', 'latex', 'FontSize', 14);
legend('Numerical', 'Analytical','Interpreter', 'latex', 'FontSize', 14, 'Location', 'east');
set(gca,'TickLabelInterpreter','latex','FontSize',14)