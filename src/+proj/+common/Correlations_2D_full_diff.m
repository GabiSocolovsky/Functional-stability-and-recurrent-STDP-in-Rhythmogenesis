%% Correlations_2D
function [Corr_ie_final,Corr_ei_final,Corr_ee_final,Corr_ii_final,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended] = Correlations_2D_full_diff(J_ee,J_ei,J_ie,J_ii,dt,tf)
%% History of m , -D<t<0 %%
m_e_history=0.1;
m_i_history=0.5;
N_e=size(J_ei,1); % size of excitatory poulation
N_i=size(J_ei,2); % size of inhibitory population
%m_history=0.0000001; % The value of m in the time -D<=t<0
%time_hist=-D:dt:(0-dt);
%m_history=0.5+0.4*cos(100*time_hist);
%m_history=rand(1,length(time_hist));
%% Calculation of the Correlation of m
global Delta_min Delta_max

[m_e,m_i,T_mean_m_e,T_mean_m_i,time]=proj.common.Two_populations_full_rate_model_history(m_e_history,m_i_history,J_ee,J_ei,J_ie,J_ii,dt,tf);
Delta_min=-200;
Delta_max=200;

T_mean_m_e=round(T_mean_m_e,-log10(dt)); %leaving only valuable digits in accordance with dt!
T_mean_m_i=round(T_mean_m_e,-log10(dt)); %leaving only valuable digits in accordance with dt!
Corr_ie_final=zeros(N_i,N_e,length(Delta_min:dt:Delta_max));
Corr_ei_final=zeros(N_e,N_i,length(Delta_min:dt:Delta_max));

t_f_ind=length(time); % The initial time at which the corrleation is calculated

if T_mean_m_e>0
    m_e_T=m_e(:,(t_f_ind-round(T_mean_m_e/dt)+1):t_f_ind); %
    m_i_T=m_i(:,(t_f_ind-round(T_mean_m_e/dt)+1):t_f_ind); %
else
    m_e_T=nan; %
    m_i_T=nan; %    
end

for n=1:N_e
    for k=1:N_i
        
        if isnan(T_mean_m_e(n))  %  When m_e is non-oscillatory (than m_i is not oscillatory too), T_mean_m_e=T_mean_m_i=NaN
            if isnan(T_mean_m_i(k))
                corr=m_e(n,end)*m_i(k,end); % in the case of non oscillatory m, I  choose an arbitrary slot, since every slot has the mean value (from the "My_delayed...")
                Corr_ie_final(k,n,:)=corr;
                Corr_ei_final(n,k,:)=Corr_ie_final(k,n,:);
            else
                Corr_ie_final(k,n,:)=0;
                Corr_ei_final(n,k,:)=Corr_ie_final(k,n,:);
            end
        elseif isnan(T_mean_m_i(k))
            Corr_ie_final(k,n,:)=0;
            Corr_ei_final(n,k,:)=Corr_ie_final(k,n,:);
        else
            
            m_e_t=m_e(n,(t_f_ind-round(T_mean_m_e(n)/dt)+1):t_f_ind); %
            m_i_t=m_i(k,(t_f_ind-round(T_mean_m_e(n)/dt)+1):t_f_ind); %
            
            
            corr_ei=fliplr(((toeplitz([m_i_t(1) fliplr(m_i_t(2:end))],m_i_t)*m_e_t')./length(m_e_t))');
            
            
            corr_ei=repmat(corr_ei,1,round(length(0:dt:Delta_max)/round(T_mean_m_e(n)/dt))+2);
            corr_ei=[(corr_ei(2:end)) corr_ei];
            
            corr_ei=corr_ei(round((length(corr_ei)-length(Delta_min:dt:Delta_max))/2)+1:(length(corr_ei)-(length(corr_ei)-length(Delta_min:dt:Delta_max))/2));
            
            Corr_ei_final(n,k,:)=corr_ei;
            Corr_ie_final(k,n,:)=fliplr(corr_ei);
        end
    end
end

if N_e==1 && N_i==1 % autocorrelations only for ORDER PARAMETERS (when Ne=1 and Ni=1)
    if isnan(T_mean_m_e)  %  When m_e is non-oscillatory (than m_i is not oscillatory too), T_mean_m_e=T_mean_m_i=NaN
        if isnan(T_mean_m_i)
            Corr_ee_final(:)=m_e(end)*m_e(end)*ones(1,length(Delta_min:dt:Delta_max)); % in the case of non oscillatory m, I  choose an arbitrary slot, since every slot has the mean value (from the "My_delayed...")
            Corr_ii_final(:)=m_i(end)*m_i(end)*ones(1,length(Delta_min:dt:Delta_max)); % in the case of non oscillatory m, I  choose an arbitrary slot, since every slot has the mean value (from the "My_delayed...")
        else
            Corr_ee_final(:)=0*ones(1,length(Delta_min:dt:Delta_max));
            Corr_ii_final(:)=0*ones(1,length(Delta_min:dt:Delta_max));
        end
    elseif isnan(T_mean_m_i(k))
            Corr_ee_final(:)=0*ones(1,length(Delta_min:dt:Delta_max));
            Corr_ii_final(:)=0*ones(1,length(Delta_min:dt:Delta_max));
    else
        
        m_e_t=m_e((t_f_ind-round(T_mean_m_e(n)/dt)+1):t_f_ind); %
        m_i_t=m_i((t_f_ind-round(T_mean_m_e(n)/dt)+1):t_f_ind); %
        
        
        corr_ee=fliplr(((toeplitz([m_e_t(1) fliplr(m_e_t(2:end))],m_e_t)*m_e_t')./length(m_e_t))');
        corr_ee=repmat(corr_ee,1,round(Delta_max/T_mean_m_e)+1);
        corr_ee=corr_ee(1:round(Delta_max/dt)+1);
        Corr_ee_final=[fliplr(corr_ee(2:end)) corr_ee];
        
        corr_ii=fliplr(((toeplitz([m_i_t(1) fliplr(m_i_t(2:end))],m_i_t)*m_i_t')./length(m_i_t))');
        corr_ii=repmat(corr_ii,1,round(Delta_max/T_mean_m_i)+1);
        corr_ii=corr_ii(1:round(Delta_max/dt)+1);
        Corr_ii_final=[fliplr(corr_ii(2:end)) corr_ii];
    end
else
    Corr_ee_final=0; % we assume no STDP dynanmics of individual intra-synapses
    Corr_ii_final=0; % we assume no STDP dynanmics of individual intra-synapses
end

%figure;
Delta_extended=-((length(Corr_ei_final(n,k,:))-1)*dt/2):dt:((length(Corr_ei_final(n,k,:))-1)*dt/2);
%plot(Delta_extended,Corr_ei,'LineWidth',1.5);
%figure;
%plot(Delta_extended,Corr_ie,'LineWidth',1.5);
%xlabel('$\Delta$','Interpreter','latex','FontSize',18)
%ylabel('$\Gamma$','Interpreter','latex','FontSize',19)

end
