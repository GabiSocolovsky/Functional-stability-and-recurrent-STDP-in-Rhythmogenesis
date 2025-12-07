%% Figure 9 Simulation %%
% This script generates **Figure 9** from the paper:
% "Functional stability and recurrent STDP in rhythmogenesis"
%
% Description:

%       - Computes the continuous STDP order parameter dynamics for sparse
%       connectivity

% Dependencies:

%       - Two_populations_full_rate_model_sparse.m - Continuous neural
%       dynamics with sparse connectivity

% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Basic Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Definitions %%%%%%%%%%%%%
filename='data\SparseDynamicsp0.3.mat'; % 
oldparam = sympref("HeavisideAtOrigin",1/2); % 0.5 at the origin for heaviside
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=1; % time constant tau (5msec)
dt=0.1; % timebin
Tunits=5*10^-3; % 5ms for T=1
tf=200; % final time of neural dynamics
D=0.4; % delay in msec
Jiimean=0.4; % Jii order parameter
Jeemean=0.6; % Jee order parameter
Jhat=(Jeemean+Jiimean)/2;
%%%%%%% Bifurcation properties %%%%%%%
syms wD JbarD
if Jeemean>=Jiimean
    range=[0.1 5 ;0.01 pi/(2*D)];
    Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), T*wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
elseif Jeemean<Jiimean
    range=[0.1 5 ;0.01 pi/D];
    Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), T*wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
end
wD=double(Y.wD);
fD=wD/(2*pi);
JbarD=double(Y.JbarD);
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% STDP and grid Parameters %%%%%%%%%
N_e=60;
N_i=15;
alpha=0.98; % relative depression
mu=0.01;% measure of linearity
Jiemax=20; % J_ie_max
tau_pE=2; % typical potentiation time of excitatory synapses
tau_pI=2; % typical potentiation time of inhibitory synapses
tau_mE=5; % typical depression time of excitatory synapses
tau_mI=3; % typical depression time of inhibitory synapses
thetapI=acos((1+(wD*tau_pI)^2)^-0.5);
thetamI=acos((1+(wD*tau_mI)^2)^-0.5);
thetapE=acos((1+(wD*tau_pE)^2)^-0.5);
thetamE=acos((1+(wD*tau_mE)^2)^-0.5);
lambda_e=8*10^-1; % excitatory synapses learning rate
lambda_i=8*10^-2; % inhibitory synapses learning rate
H_E=-1;
H_I=1;

colo=[0.2, 0.6, 0.2]; % Green

p=0.3; % sparsness of connectivity
Jeimean=0.1; % inital Jei order parameter
Jiemean=0.3; % initial Jie order parameter
f=(1-Jiemean/Jiemax)^mu; % synaptic dependent function
%%%%%%%%% STDP analytical expression %%%%%%%%%
K_Ibar=1-alpha;
K_Iptilphi=cos(thetapI)*cos(thetapI+phi);
K_Imtilphi=cos(thetamI)*cos(thetamI-phi);
K_Itilphi=K_Iptilphi-alpha*K_Imtilphi;
K_Iptilmphi=cos(thetapI)*cos(thetapI-phi);
K_Imtilmphi=cos(thetamI)*cos(thetamI+phi);
K_Itilmphi=K_Iptilmphi-alpha*K_Imtilmphi;

K_Iptilpsi=cos(thetapI)*cos(thetapI+psi);
K_Imtilpsi=cos(thetamI)*cos(thetamI-psi);
K_Itilpsi=K_Iptilpsi-alpha*K_Imtilpsi;

K_Ebar=f-alpha;
K_Eptilphi=cos(thetapE)*cos(thetapE-phi);
K_Emtilphi=cos(thetamE)*cos(thetamE+phi);
K_Etilphi=f*K_Eptilphi-alpha*K_Emtilphi;
K_Eptilmphi=cos(thetapE)*cos(thetapE+phi);
K_Emtilmphi=cos(thetamE)*cos(thetamE-phi);
K_Etilmphi=f*K_Eptilmphi-alpha*K_Emtilmphi;

K_Eptilmpsi=cos(thetapE)*cos(thetapE+psi);
K_Emtilmpsi=cos(thetamE)*cos(thetamE-psi);
K_Etilmpsi=f*K_Eptilmpsi-alpha*K_Emtilmpsi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gabs=(JbarD^2-Jeemean*Jiimean)^0.5; % absolute value of g

%%%%%%%%%%%%%%% Make the inital connectivity matrix %%%%%%%%%%%%%
J_ei_std=Jeimean/10; % SD of Jei
J_ei_noise=normrnd(0,J_ei_std,N_e,N_i); % noise part of Jei
Jei=(Jeimean*ones(N_e,N_i))-mean(J_ei_noise(:))+J_ei_noise; % Jei connectivity matrix
Jei_bin=double(rand(N_e,N_i)>p); % Matrix of living (1) and dead (0) Jei synapses - Bernoulli distributed with probability p 
Jei=Jei.*Jei_bin; % Final matrix of Jei

J_ie_std=Jiemean/10; % SD of Jie
J_ie_noise=normrnd(0,J_ie_std,N_i,N_e); % noise part of Jie
Jie=Jiemean*ones(N_i,N_e)-mean(J_ie_noise(:))+J_ie_noise; % Jie connectivity matrix
Jie_bin=double(rand(N_i,N_e)>p);  % Matrix of living (1) and dead (0) Jie synapses - Bernoulli distributed with probability p 
Jie=Jie.*Jie_bin; % Final matrix Jie

Jee_bin=double(rand(N_e,N_e)>p); % Matrix of living... for Jee
Jee=Jeemean*ones(N_e,N_e).*Jee_bin; % final matrix of Jee
Jii_bin=double(rand(N_i,N_i)>p); % Matrix of living... for Jee
Jii=Jiimean*ones(N_i,N_i).*Jii_bin; % final matrix of Jii

Jeiinitlive=sum(Jei_bin(:)); % Number of inital living Jei
Jieinitlive=sum(Jie_bin(:)); % ... for Jie

tflearn=76000; % total STDP running time
t=1:tflearn; % time for STDP dynamics
J_ei_dot=zeros(N_e,N_i);
J_ie_dot=zeros(N_i,N_e);
dtlearn=1; % timebin for STDP dynamics

Jeimeandyn=[]; % mean Jei through time
Jiemeandyn=[]; % ... Jie
Jiedyn=[]; % Jei through time
Jeidyn=[]; % ... same for Jie

m_e_history=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean); % first history of m_e
m_i_history=(1+Jiemean-Jeemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean); % first history of m_i

m_e_full=[]; % create null m_e for full dynamics
m_i_full=[]; % create null m_i for full dynamics
m_e_full_std=[]; % SD of m_e through time
m_i_full_std=[]; % SD of m_i through time
time_full=[]; % full time

%%%%%% Create matrices of temporal kernels for integration in real time %%%
tmaxSTDP=10*max([tau_mE tau_mI tau_pI tau_pE]); % history time for m for calculating STDP
tSTDP=dt:dt:(tf + tmaxSTDP); % time for STDP
[TSTDP,~]=meshgrid(tSTDP,tSTDP); % make a matrix of these times

Delta=TSTDP-TSTDP'; % Make Delta
KpI=1/tau_pI*exp(-Delta*H_I/tau_pI).*heaviside(H_I*Delta); % temporal kernel, p - plus , m - minus
KmI=1/tau_mI*exp(Delta*H_I/tau_mI).*heaviside(-H_I*Delta);
KpE=1/tau_pE*exp(-Delta*H_E/tau_pE).*heaviside(H_E*Delta);
KmE=1/tau_mE*exp(Delta*H_E/tau_mE).*heaviside(-H_E*Delta);
KI=KpI-alpha*KmI; % we define it outside the dynamics because it does not has a dependence in the synaptic dynamics

Transtime=[]; % The time before each synaptic update
timebinsave=round((tf*tflearn/dt)/50); % When to save


for i=2:length(t) % dynamics
    
    %%% construcing the eigenvectors on each time step %%%
    Jeimean=sum(Jei(:))/Jeiinitlive; % mean synapse Jei (order parameter)
    Jiemean=sum(Jie(:))/Jieinitlive; % mean synapse Jie (order parameter)
    
    Jeimeandyn(i)=Jeimean;
    Jiemeandyn(i)=Jiemean;
    Jeidyn(:,:,i)=Jei;
    Jiedyn(:,:,i)=Jie;
    
            
    [m_e_now,m_i_now,T_mean_m_e,~,tnetwork]=Two_populations_full_rate_model_sparse(m_e_history,m_i_history,Jee,Jei,Jie,Jii,dt,tf,p); % compute the neuronal dynamics m_e and m_i
    
    m_e_now=m_e_now(:,1:end-1); % m_e of this current loop
    m_i_now=m_i_now(:,1:end-1); % m_i of this current loop
    m_e_full_std=[m_e_full_std  std(m_e_now)];
    m_i_full_std=[m_i_full_std  std(m_i_now)];
    
    
    m_e_history=m_e_now(:,end-round(D/dt):end); % that is the history for the m_e and m_i for the calculation of the network dynamics
    m_i_history=m_i_now(:,end-round(D/dt):end); % that is the history for the m_e and m_i for the calculation of the network dynamics
    
    
    %%%%%%% Take the m_e and m_i from previous steps in order to take them into
    % Account in the STDP dynamics %%%%%%%%%%%%%%%%%%%%%%%
    if i==2
        time_full=tnetwork;
        m_e_bef=m_e_now(:,1).*ones(N_e,length(0:dt:tmaxSTDP-dt));
        m_i_bef=m_i_now(:,1).*ones(N_i,length(0:dt:tmaxSTDP-dt));
    else
        time_full=[time_full time_full(end)+tnetwork];
        m_e_bef=m_e_full(:,(end-tmaxSTDP/dt+1):end);
        m_i_bef=m_i_full(:,(end-tmaxSTDP/dt+1):end);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Transtime(i)=time_full(end); % the times between each synaptic update
    
    Jeidot=lambda_i*dt*( [m_i_bef m_i_now]*KI*[m_e_bef m_e_now]' - [m_i_bef 0*m_i_now]*KI*[m_e_bef 0*m_e_now]' )'/tf; % compute the integral in STDP dynamics (see appendix of numerics)
    Jiedotpot=lambda_e*dt*( (proj.common.Th_li_full(1-Jie'/Jiemax)).^mu.*([m_e_bef m_e_now]*KpE*[m_i_bef m_i_now]') - (proj.common.Th_li_full(1-Jie'/Jiemax)).^mu.*([m_e_bef 0*m_e_now]*KpE*[m_i_bef 0*m_i_now]') )'/tf; % compute the integral in STDP dynamics (see appendix of numerics)
    Jiedotdep=lambda_e*dt*( [m_e_bef m_e_now]*KmE*[m_i_bef m_i_now]' - [m_e_bef 0*m_e_now]*KmE*[m_i_bef 0*m_i_now]' )'/tf; % compute the integral in STDP dynamics (see appendix of numerics)
    Jiedot=Jiedotpot-alpha*Jiedotdep; % same 
    
    m_e_full=[m_e_full m_e_now]; % Add the current dynamics to the previous dynamics
    m_i_full=[m_i_full m_i_now]; % Add the current dynamics to the previous dynamics
    
    Jei=(Jei+dtlearn*Jeidot).*Jei_bin.*heaviside((Jei+dtlearn*Jeidot).*Jei_bin); % update synapse
    Jie=(Jie+dtlearn*Jiedot).*Jie_bin.*heaviside((Jie+dtlearn*Jiedot).*Jie_bin); % update synapse
    
    if size(m_e_full,2)>timebinsave
        
        m_e_full=m_e_full(:,end-timebinsave+1:end);
        m_i_full=m_i_full(:,end-timebinsave+1:end);
        
        if ~mod(i,1000)
            
            save(filename,"Jeimeandyn","Jiemeandyn","JbarD","Jiimean","Jiemax");
            i % a flag for last save
        end
        
        
    end
    
    
end
%% Plot figure 9 %%
figure(2)
plot(0:0.01:2,(JbarD^2)./(0:0.01:2),'LineWidth',2.5,'Color','Black')
hold on
plot(Jeimeandyn(2:end),Jiemeandyn(2:end),'.','Color',[0 0.4470 0.7410])
plot((1+Jiimean)*ones(1,length(0:1:Jiemax)),0:1:Jiemax,'LineWidth',2.5,'Color','Black')
plot(Jeimeandyn(2),Jiemeandyn(2),'+','Color',[0.8500, 0.3250, 0.0980],'MarkerSize', 12,'LineWidth',2,'HandleVisibility', 'off'); % No legend entry
plot(Jeimeandyn(end),Jiemeandyn(end),'x','Color',[0, 1, 0],'MarkerSize', 12,'LineWidth',2,'HandleVisibility', 'off'); % No legend entry
xlabel('$J_{EI}$','interpreter','latex','FontSize',18)
ylabel('$J_{IE}$','interpreter','latex','FontSize',18)
xlim([0 1.5])
ylim([0 Jiemax])
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
grid on