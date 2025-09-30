%% Two populations rate model %%%
function [m_e,m_i,T_mean_m_e,T_mean_m_i,time] = Two_populations_full_rate_model_history(m_e_history,m_i_history,J_ee,J_ei,J_ie,J_ii,dt,tf)
% This function is the simulation of the neural dynamics in the paper:
% "Functional stability and recurrent STDP in Rhythmogenesis"

%   Description:

%       Computes the dynamics of the firing rates of the neural network by using
%       Euler's method for ODE. In addition it computes the time period in the
%       case of Rhythmic activity

%   Inputs:

%       m_e_history   -   m_e at time -d<t<0
%       m_i_history   -   m_i at time -d<t<0
%       J_ee          -   E-to-E connectivity matrix
%       J_ei          -   I-to-E connectivity matrix
%       J_ie          -   E-to-I connectivity matrix
%       J_ii          -   I-to-I connectivity matrix
%       dt            -   time bin
%       tf            -   final time of simulation

%   Outputs:

%       m_e                      -   Dynamics of Excitatory neurons firing rates
%       m_i                      -   Dynamics of Inhibitory neurons firing rates
%       T_mean_m_e,T_mean_m_i    -   Time period of m_e and m_i (sometimes
%       different when in only the inhibitory neurons are rhythmic)
%       time                     -   Vector of time

%   Dependencies:

%   Period_Peaks_full.m          -   Computes Time period when in the R region
%
% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Network parameters %%%%%%%
I_0=1; % Constant external input
tau=1; % Time scale for the firing rate dynamics - 5ms
D=0.4; % Delay

t0=0; % Initial time
time=t0:dt:tf; % Time vector

N_e=size(J_ei,1); % Size of excitatory poulation
N_i=size(J_ei,2); % Size of inhibitory population


m_e=zeros(N_e,length(time)); % Empty array of excitatory firing rates dynamics
m_i=zeros(N_i,length(time)); % Empty array of Inhibitory firing rates dynamics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Network Dynamics %%%%%%%%%%%%%%%%%%
if (size(m_e_history,2)==1) && length(size(m_i_history,2))==1  % checks if the history function is a constant value in time
    
    if (size(m_e_history,1)==1) && length(size(m_i_history,1))==1 % checks if the history is one value that represents the firing rates of all neurons
        m_e_history=m_e_history*ones(1,N_e);
        m_i_history=m_i_history*ones(1,N_i);
    else
        m_e_history=m_e_history';
        m_i_history=m_i_history';
    end
    m_e(:,1)=m_e_history + (dt/tau)*(-m_e_history+Th_li_full(I_0 + (J_ee*m_e_history')'/N_e - (J_ei*m_i_history')'/N_i)); % first step
    m_i(:,1)=m_i_history + (dt/tau)*(-m_i_history+Th_li_full(I_0 + (J_ie*m_e_history')'/N_e - (J_ii*m_i_history')'/N_i)); % first step
    I_before_D_e=Th_li_full(I_0 + (J_ee*m_e_history')/N_e - (J_ei*m_i_history')/N_i); % input to m_e when -D<t<0
    I_before_D_i=Th_li_full(I_0 + (J_ie*m_e_history')/N_e - (J_ii*m_i_history')/N_i); % input to m_i when -D<t<0
    
    i=1;
    
    for t=(t0+dt):dt:tf % The dynamics loop
        
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

[T_mean_m_e,T_std_m_e]=proj.common.Period_Peaks_full(m_e,dt);  % compute the time period of m_e
[T_mean_m_i,T_std_m_i]=proj.common.Period_Peaks_full(m_i,dt);  % compute the time period of m_i

m_e(m_e(:,:,end)<10^-7)=0; % Don't work with small numbers
m_i(m_i(:,:,end)<10^-7)=0; % Don't work with small numbers

end