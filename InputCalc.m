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
figure;
plot(x_plot, y_plot, '-'); % plot the fitted line
hold on
scatter(u_ref_cfd, u_inp_cfd, 'red','filled','o')
ylabel('Input Velocity (m/s)', 'Interpreter','latex', 'FontSize',14)
xlabel('Pitot Static Tube Velocity (m/s)', 'Interpreter','latex', 'FontSize',14)
set(gca,'TickLabelInterpreter','latex','FontSize',14)
hold off;

%Find the point
u_sc_actual = polyval(p, u_ts);
percent_error = (abs(u_ts-67.3067)/(u_sc_actual))*100;
%%
disp('Air Properties in Pitot Static Flow')
fprintf('Inlet Velocity (m/s): %.4f ', u_sc_actual);
fprintf('Initial Calculation Input Velocity (m/s): %.4f ', u_sc);
fprintf('Difference (percent): %.4f ', percent_error);
fprintf('Pitot Static Velocity (m/s): %.4f ', u_ts);
fprintf('Reference density (kg/m3): %.4f ', rho_ref);
fprintf('Kinematic viscosity (m2/s): %.4e\n', nu);