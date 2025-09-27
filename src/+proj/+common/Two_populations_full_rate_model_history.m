function [m_e,m_i,T_mean_m_e,T_mean_m_i,time] = Two_populations_full_rate_model_history(m_e_history,m_i_history,J_ee,J_ei,J_ie,J_ii,dt,tf)
%% Two populations rate model %%

I_0=1; % Constant external input
tau=1; % Time scale for the firing rate dynamics - 5ms
D=0.4; % Delay

t0=0; % Initial time
time=t0:dt:tf; % Time vector

N_e=size(J_ei,1); % size of excitatory poulation
N_i=size(J_ei,2); % size of inhibitory population


m_e=zeros(N_e,length(time));
m_i=zeros(N_i,length(time));


if (size(m_e_history,2)==1) && length(size(m_i_history,2))==1  % checks if the history function is a constant value in time

    if (size(m_e_history,1)==1) && length(size(m_i_history,1))==1 % checks if the history is one value that represents the firing rates of all neurons
        m_e_history=m_e_history*ones(1,N_e);
        m_i_history=m_i_history*ones(1,N_i);
    else
        m_e_history=m_e_history';
        m_i_history=m_i_history';
    end
m_e(:,1)=m_e_history + (dt/tau)*(-m_e_history+Th_li_full(I_0 + (J_ee*m_e_history')'/N_e - (J_ei*m_i_history')'/N_i));
m_i(:,1)=m_i_history + (dt/tau)*(-m_i_history+Th_li_full(I_0 + (J_ie*m_e_history')'/N_e - (J_ii*m_i_history')'/N_i));
I_before_D_e=Th_li_full(I_0 + (J_ee*m_e_history')/N_e - (J_ei*m_i_history')/N_i);
I_before_D_i=Th_li_full(I_0 + (J_ie*m_e_history')/N_e - (J_ii*m_i_history')/N_i);

i=1;

for t=(t0+dt):dt:tf

    if t<D
        
        m_e(:,i+1)= m_e(:,i) + (dt/tau)*(-m_e(:,i)+I_before_D_e);    
        m_i(:,i+1)= m_i(:,i) + (dt/tau)*(-m_i(:,i)+I_before_D_i);    
        Inpm_e(:,i+1)=I_before_D_e;
        Inpm_i(:,i+1)=I_before_D_i;        
    else
       
       m_e(:,i+1)= m_e(:,i) + (dt/tau)*(-m_e(:,i)+Th_li_full(I_0 + J_ee*m_e(:,round((t-D)/dt)+1)/N_e-J_ei*m_i(:,round((t-D)/dt)+1)/N_i));                     
       m_i(:,i+1)= m_i(:,i) + (dt/tau)*(-m_i(:,i)+Th_li_full(I_0 + J_ie*m_e(:,round((t-D)/dt)+1)/N_e-J_ii*m_i(:,round((t-D)/dt)+1)/N_i));
       Inpm_e(:,i+1)=I_0 + J_ee*m_e(:,round((t-D)/dt)+1)/max([N_e-1 1])-J_ei*m_i(:,round((t-D)/dt)+1)/N_i;
       Inpm_i(:,i+1)=I_0 + J_ie*m_e(:,round((t-D)/dt)+1)/N_e-J_ii*m_i(:,round((t-D)/dt)+1)/max([N_i-1 1]); 
    end
            
         i=i+1;

end




else   % if m_e_history and m_i_history are not constants
    
  
    m_e(:,1)=m_e_history(:,end) + (dt/tau)*(-m_e_history(:,end)+Th_li_full(I_0 + J_ee*m_e_history(:,end-round(D/dt))/N_e - J_ei*m_i_history(:,end-round(D/dt))/N_i));
    m_i(:,1)=m_i_history(:,end) + (dt/tau)*(-m_i_history(:,end)+Th_li_full(I_0 + J_ie*m_e_history(:,end-round(D/dt))/N_e - J_ii*m_i_history(:,end-round(D/dt))/N_i));
    
    i=1;
    

for t=(t0+dt):dt:tf
    
    if t<D
               
       m_e(:,i+1)= m_e(:,i) + (dt/tau)*(-m_e(:,i)+Th_li_full(I_0 + J_ee*m_e_history(:,round(end-round(D/dt))+i)/N_e-J_ei*m_i_history(:,round(end-round(D/dt))+i)/N_i));                     
       m_i(:,i+1)= m_i(:,i) + (dt/tau)*(-m_i(:,i)+Th_li_full(I_0 + J_ie*m_e_history(:,round(end-round(D/dt))+i)/N_e-J_ii*m_i_history(:,round(end-round(D/dt))+i)/N_i));
       Inpm_e(:,i+1)=I_0 + J_ee*m_e_history(:,round(end-round(D/dt))+i)/max([N_e-1 1])-J_ei*m_i_history(:,round(end-round(D/dt))+i)/N_i;        
       Inpm_i(:,i+1)=I_0 + J_ie*m_e_history(:,round(end-round(D/dt))+i)/N_e-J_ii*m_i_history(:,round(end-round(D/dt))+i)/max([N_i-1 1]); 
    else
       
       m_e(:,i+1)= m_e(:,i) + (dt/tau)*(-m_e(:,i)+Th_li_full(I_0 + J_ee*m_e(:,round((t-D)/dt)+1)/N_e-J_ei*m_i(:,round((t-D)/dt)+1)/N_i));                     
       m_i(:,i+1)= m_i(:,i) + (dt/tau)*(-m_i(:,i)+Th_li_full(I_0 + J_ie*m_e(:,round((t-D)/dt)+1)/N_e-J_ii*m_i(:,round((t-D)/dt)+1)/N_i));
       Inpm_e(:,i+1)=I_0 + J_ee*m_e(:,round((t-D)/dt)+1)/max([N_e-1 1])-J_ei*m_i(:,round((t-D)/dt)+1)/N_i;
       Inpm_i(:,i+1)=I_0 + J_ie*m_e(:,round((t-D)/dt)+1)/N_e-J_ii*m_i(:,round((t-D)/dt)+1)/max([N_i-1 1]);
    end       
         i=i+1;

end
    


end

[T_mean_m_e,T_std_m_e]=proj.common.Period_Peaks_full(m_e,dt); 
[T_mean_m_i,T_std_m_i]=proj.common.Period_Peaks_full(m_i,dt); 

%m_e_steady=I_0*(1-J_ei)/(1+J_ei*J_ie);
%m_i_steady=I_0*(1+J_ie)/(1+J_ei*J_ie);

m_e(m_e(:,:,end)<10^-7)=0;
m_i(m_i(:,:,end)<10^-7)=0;

%T1=Period_current(dt,I_0,J_0,m_i.');



%%% The condition for input of m_i above
%figure;
%for k=1:size(m_e,1)
%plot(time,m_e(k,:))
%hold on
%end
%title('m_E','FontSize',20)

%figure;

%for k=1:size(m_i,1)
%plot(time,m_i(k,:))
%hold on
%end
%title('m_I','FontSize',20)


%x0=0;
%y0=0;
%width=18;
%height=40;

%figure;
%subplot(2,1,1); %% of m_e
%plot(time*5,m_e,'LineWidth',2);
%hold on
%plot(J_i_arr,J_bar_D^2./J_i_arr,'color','black','Linewidth',6); % the limit between oscillatory phase to lower than
%title('J_e nullcline','FontSize',30);
%xlabel('$t $','fontsize',30,'interpreter','latex')
%ylabel('$m_E$','fontsize',30,'interpreter','latex')
%axis([1000 3000 0 0.07])

%set(gca,'FontSize',20)
%set(gcf,'Units','Centimeters','position',[x0,y0,width,height])
%set(gca,'xtick',[])




%1/T_mean_m_e*100
%subplot(2,1,2);
%plot(time*5,m_i,'LineWidth',2);
%w=2*pi/T_mean_m_e;
%f=w/(2*pi); %in KHz
%hold on
%plot(J_i_arr,J_bar_D^2./J_i_arr,'color','black','Linewidth',6); % the limit between oscillatory phase to lower than
%title('J_e nullcline','FontSize',30);
%xlabel('$t \ [\textrm{msec}]$','fontsize',30,'interpreter','latex')
%ylabel('$m_I$','fontsize',30,'interpreter','latex')
%axis([~ ~ 1 1.6])

%set(gca,'FontSize',20)
%set(gcf,'Units','Centimeters','position',[x0,y0,width,height])

% figure;
% plot(5*time,m_e,'LineWidth',3)
% hold on
% plot(5*time,m_i,'LineWidth',3)
% xlabel('$t \ [\mathrm{ms}]$','interpreter','latex','FontSize',18)
% ylabel('$m_X(t)$','interpreter','latex','FontSize',18)
% lgd=legend({'$\mathrm{E}$','$\mathrm{I}$'},'Interpreter','latex','Location','NorthEast');
% legendJiiTitle=sprintf('$X$','Interpreter','latex');
% lgd.Title.String = legendJiiTitle;
% set(gca,'FontSize',14)
% set(gca,'TickLabelInterpreter','latex')
% grid on

%plot(Inpm_e)
%hold on
%plot(Inpm_i)
end