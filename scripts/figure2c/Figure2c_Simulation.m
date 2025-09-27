dt=0.01;
tf=200;
m_e_history=0.2;
m_i_history=0.8;
J_ei=0.6;
J_ie=10;
J_ee=0.6;
J_ii=0.4;


[m_e,m_i,T_mean_m_e,T_mean_m_i,time]=proj.common.Two_populations_full_rate_model_history(m_e_history,m_i_history,J_ee,J_ei,J_ie,J_ii,dt,tf);

figure;
plot(5*time,m_i,'LineWidth',2.5)
hold on
plot(5*time,m_e,'LineWidth',2.5)
xlabel('$t \ [\mathrm{ms}]$','interpreter','latex','FontSize',18)
ylabel('$m_X(t)$','interpreter','latex','FontSize',18)
xlim([0 250])
lgd=legend({'$\mathrm{I}$','$\mathrm{E}$'},'Interpreter','latex','Location','NorthEast');
legendJiiTitle=sprintf('$X$','Interpreter','latex');
lgd.Title.String = legendJiiTitle;
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
grid on