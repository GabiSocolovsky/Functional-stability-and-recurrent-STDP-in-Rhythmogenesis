%% Figure2c Simulation %%
% This script generates **Figure 2c** from the paper:
% "Functional stability and recurrent STDP in rhythmogenesis"
%
% Description:

%       Computes and plots the dynamics of the order parameters m_E and m_I
%       in a specific point in the R region

% Dependencies:

%       - Two_populations_full_rate_model_history.m
%
% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Parameters %%%%%%%%%%%%%%%
dt=0.01; % time bin
tf=200; % final time
m_e_history=0.2; % value of m_e when -D<t<0
m_i_history=0.8; % value of m_i when -D<t<0
J_ei=0.6; % J_EI order parameter synapse
J_ie=10; % J_IE order parameter synapse
J_ee=0.6; % J_EE order parameter synapse
J_ii=0.4; % J_II order parameter synapse

%%%%%%%%%%%% Dyanmics %%%%%%%%%%%%%%%%%%
[m_e,m_i,T_mean_m_e,T_mean_m_i,time]=proj.common.Two_populations_full_rate_model_history(m_e_history,m_i_history,J_ee,J_ei,J_ie,J_ii,dt,tf); % computes the dynamics of m_E and m_I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%
figure;
plot(5*time,m_i,'LineWidth',2.5) % multiplying by 5 since the units are 5ms, tau_m=1 (5ms), and we want to display units of ms
hold on
plot(5*time,m_e,'LineWidth',2.5) % same here
xlabel('$t \ [\mathrm{ms}]$','interpreter','latex','FontSize',18)
ylabel('$m_X(t)$','interpreter','latex','FontSize',18)
xlim([0 250])
lgd=legend({'$\mathrm{I}$','$\mathrm{E}$'},'Interpreter','latex','Location','NorthEast');
legendJiiTitle=sprintf('$X$','Interpreter','latex');
lgd.Title.String = legendJiiTitle;
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
grid on