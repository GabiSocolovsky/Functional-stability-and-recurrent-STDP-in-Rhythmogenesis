%% The MEAN synaptic dynamics - figure 3b lambda->0 %%
T=1; % time constant tau in 5msec
Tunits=5*10^-3; % 5ms for T=1
dt=0.01;
tf=200;
D=0.4; % delay in msec
Jii=0.4;
Jee=0.6;
Jhat=(Jee+Jii)/2;
syms wD JbarD
if Jee>=Jii
range=[0.1 5 ;0.01 pi/(2*D)];
Y=vpasolve([(JbarD^2-Jee*Jii)^0.5==1/cos(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5))), T*wD==-tan(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
elseif Jee<Jii
range=[0.1 5 ;0.01 pi/D];    
Y=vpasolve([(JbarD^2-Jee*Jii)^0.5==1/cos(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5))), T*wD==-tan(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
end    
wD=double(Y.wD);
JbarD=double(Y.JbarD);
phi=acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5));
psi=acos(Jhat/JbarD);
%%%%%%%%% STDP Parameters %%%%%%%%%
alpha_arr=[0.95:0.01:0.98 0.98]; % relative depression %
Jei_arr=[0.1 0.2 0.5 0.1 0.75];
Jie_arr=[10 9 10 5 6.5];
mu=0.01;%0.015;%0.15;%0.015; % measure of linearity
Jiemax=20; % J_ie_max
tau_pE=2; % typical potentiation time of excitatory synapses
tau_pI=2; % typical potentiation time of inhibitory synapses
tau_mE=5; % typical depression time of excitatory synapses
tau_mI=3; % typical depression time of inhibitory synapses
thetapI=acos((1+(wD*tau_pI)^2)^-0.5);
thetamI=acos((1+(wD*tau_mI)^2)^-0.5);
thetapE=acos((1+(wD*tau_pE)^2)^-0.5);
thetamE=acos((1+(wD*tau_mE)^2)^-0.5);
lambda_e=4;%0*200*15/1000; % excitatory synapses learning rate
lambda_i=0.4; % inhibitory synapses learning rate
lamf=1; % This is a scaling factor for the lambdas - in the FP regiion its 1, in the R region its lamf*some factor (0<some factor<1), it decreases tf as well
H_E=-1;
H_I=1;
%%%%%%%%% Generating the synapses matrix %%%%%%%%%
N_e=1; % Size of excitatory population
N_i=1; % Size of inhibitory population
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% zero order synapses %%%
Jei=0.1;%JbarD^2/J_ie_final;
Jie=10;%J_ie_final;%JbarD^2/Jeimean+0.1;%10^-5;
[Corr_ie_final,Corr_ei_final,~,~,~,~,~,~,Delta_extended]=proj.common.Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf);
K_pE=1/tau_pE*exp(-Delta_extended*H_E/tau_pE).*heaviside(H_E*Delta_extended);
K_mE=1/tau_mE*exp(Delta_extended*H_E/tau_mE).*heaviside(-H_E*Delta_extended);
K_pI=1/tau_pI*exp(-Delta_extended*H_I/tau_pI).*heaviside(H_I*Delta_extended);
K_mI=1/tau_mI*exp(Delta_extended*H_I/tau_mI).*heaviside(-H_I*Delta_extended);


% Define improved yellow and orange
orange = [1.0, 0.5, 0.0];  % Rich orange
green  = [0.0, 0.7, 0.3];  % Natural green

plot(JbarD^2./(0:0.1:20),0:0.1:20,'Color','Black','LineWidth',3)
hold on
plot((1+Jii)*ones(1,length(0:20)),0:20,'Color','Black','LineWidth',3)
xlabel('$J_{EI}$','interpreter','latex','FontSize',18)
ylabel('$J_{IE}$','interpreter','latex','FontSize',18)


for p=1:length(alpha_arr)
    
    alpha=alpha_arr(p);
    Jei=Jei_arr(p);
    Jie=Jie_arr(p);    
    K_Ibar=1-alpha;
    K_Iptilpsi=cos(thetapI)*cos(thetapI+psi);
    K_Imtilpsi=cos(thetamI)*cos(thetamI-psi);
    K_Itilpsi=K_Iptilpsi-alpha*K_Imtilpsi;
    K_Eptilpsi=cos(thetapE)*cos(thetapE+psi);
    K_Emtilpsi=cos(thetamE)*cos(thetamE-psi);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t=1:10000;
    dtlearn=1;
    
    J_ie_final=Jiemax*(1-((1-K_Ibar)*(K_Ibar*K_Emtilpsi-K_Itilpsi)/(K_Ibar*K_Eptilpsi-K_Itilpsi))^(1/mu)); % the place the    
    plot(JbarD^2/J_ie_final,J_ie_final,'*','Color',[0 0 1],'MarkerSize',14,'LineWidth',1.6)

    
    
    % Normalized interpolation factor
    colgrad = (min(p,4) - 1) / (length(alpha_arr) - 1);  % Ranges from 0 to 1
    color_k = (1 - colgrad) * orange + colgrad * green;
    
    plot(Jei,Jie,'.','Color',color_k)
    
    Jeimeanarr(1)=Jei;
    Jiemeanarr(1)=Jie;
    
    
    
    for i=2:length(t)
        
        xlim([0 2])
        ylim([0 20])
        set(gca,'FontSize',14)
        set(gca,'TickLabelInterpreter','latex')
        grid on
        
        drawnow
        
              
        if Jei*Jie<JbarD^2 % in the FP regime we solve exactly
            me=(1+Jii-Jei)/((1-Jee)*(1+Jii)+Jei*Jie);
            mi=(1+Jie-Jee)/((1-Jee)*(1+Jii)+Jei*Jie);
            Corr_ei=me*mi;
            Corr_ie=Corr_ei';
            J_ei_dot=lambda_i*((1-alpha)*Corr_ei); % NOTICE that I put Corr_ei instead of Corr_ie because there is symmetry Cei,ij=Cie,ji in the FP regime
            J_ie_dot=lambda_e*(((Th_li_full(1-Jie/Jiemax)).^mu-alpha).*Corr_ie); % NOTICE that I put Corr_ei instead of Corr_ie because there is symmetry Cei,ij=Cie,ji in the FP regime
        else % The R regime is solved with simulations
            [Corr_ie,Corr_ei,~,~,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended]=proj.common.Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf/(3*lamf));
            for n=1:N_e
                for k=1:N_i
                    J_ei_dot(n,k)=lambda_i*lamf*((K_pI-alpha*K_mI)*squeeze(Corr_ie(k,n,:))*dt);
                    J_ie_dot(k,n)=lambda_e*lamf*(((Th_li_full(1-Jie(k,n)/Jiemax))^mu*K_pE-alpha*K_mE)*squeeze(Corr_ei(n,k,:))*dt);
                end
            end
            lamf=1;
        end
        
        Jei=Jei+dtlearn*J_ei_dot;
        Jie=Jie+dtlearn*J_ie_dot;
        
        Jeimeanarr(i)=Jei;
        Jiemeanarr(i)=Jie;
        
        
        plot(Jei,Jie,'.','Color',color_k ,'HandleVisibility', 'off')
        
        if ((Jeimeanarr(i)-Jeimeanarr(i-1))^2+(Jiemeanarr(i)-Jiemeanarr(i-1))^2)^0.5<0.0001
            break;
        end
    end
end
%%% Saving of figure %%%
fig=figure(1);
folderPath = 'figures\OrderParametersDynamics';  % change this to your target folder
fileName = 'STDPOrderParametersPhDiag.fig';  % name of your .fig file
fullPath = fullfile(folderPath, fileName);  % combine folder and file name
saveas(fig,fullPath)
%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(1)
lgd=legend({'','','','$0.95$','','$0.96$','','$0.97$','','$0.98$'},'Interpreter','latex','Location','NorthEast');
lgd.Title.String = '$\alpha$';
lgd.Title.Interpreter = 'latex';
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
grid on

% % create a new pair of axes inside current figure
% axes('position',[.65 .175 .25 .25])
% box on % put box around new pair of axes
% indexOfInterestJei = (Jei < 0.45 & Jei > 0.4); % range of t near perturbation
% indexOfInterestJie = (Jie < 7.7 & Jie > 7.35); % range of t near perturbation
% 
% plot(Jei(indexOfInterestJei),Jie(indexOfInterestJie)) % plot on new axes
% axis tight

%% The STDP dynamics of order parameters - continuous dynamics %%
T=1; % time constant tau in 5msec
dt=0.01;
Tunits=5*10^-3; % 5ms for T=1
tf=200;
D=0.4; % delay in msec
Jiimean=0.4;
Jeemean=0.6;
Jhat=(Jeemean+Jiimean)/2;
syms wD JbarD
if Jeemean>=Jiimean
range=[0.1 5 ;0.01 pi/(2*D)];
Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), T*wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
elseif Jeemean<Jiimean
range=[0.1 5 ;0.01 pi/D];    
Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), T*wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
end    
wD=double(Y.wD);
fD=wD/(2*pi); % in
JbarD=double(Y.JbarD);
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);
%%%%%%%%% STDP Parameters %%%%%%%%%
N_e=1; % It is the same for analysing the order parameters in small fluctuations in individual synapses
N_i=1;

alpha=0.98; % relative depression %
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
lambda_e=0.3;%0*200*15/1000; % excitatory synapses learning rate
lambda_i=0.03; % inhibitory synapses learning rate
lamf=1; % This is a scaling factor for the lambdas - in the FP regiion its 1, in the R region its lamf*some factor (0<some factor<1), it decreases tf as well
H_E=-1;
H_I=1;
%%% zero order synapses %%%
Jeimean=0.54;%0.32;%0.4;
Jiemean=5.8;%10;%5;%JbarD^2/Jeimean+0.1;%10^-5;
f=(1-Jiemean/Jiemax)^mu;
ftag=-mu/Jiemax*(1-Jiemean/Jiemax)^(mu-1);
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
%
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

gabs=(JbarD^2-Jeemean*Jiimean)^0.5;

Jei=Jeimean;
Jie=Jiemean;
Jee=Jeemean*ones(N_e,N_e);
Jii=Jiimean*ones(N_i,N_i);
J=[[Jee/N_e -Jei/N_i];[Jie/N_e -Jii/N_i]]; % connectivity matrix
diagonal=eye(N_e+N_i);
mfa=(diagonal-J)\ones(N_e+N_i,1); % The exact fixed point solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Corr_ie_final,Corr_ei_final,~,~,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended]=proj.common.Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf);
K_pE=1/tau_pE*exp(-Delta_extended*H_E/tau_pE).*heaviside(H_E*Delta_extended);
K_mE=1/tau_mE*exp(Delta_extended*H_E/tau_mE).*heaviside(-H_E*Delta_extended);
K_pI=1/tau_pI*exp(-Delta_extended*H_I/tau_pI).*heaviside(H_I*Delta_extended);
K_mI=1/tau_mI*exp(Delta_extended*H_I/tau_mI).*heaviside(-H_I*Delta_extended);


t=1:500;
%check_arr=[reshape(check_ei.',1,[]) reshape(check_ie.',1,[])];
J_ei_dot=zeros(N_e,N_i);
J_ie_dot=zeros(N_i,N_e);
dtlearn=1;

J_ie_final=Jiemax*(1-((1-K_Ibar)*(K_Ibar*K_Emtilmpsi-K_Itilpsi)/(K_Ibar*K_Eptilmpsi-K_Itilpsi))^(1/mu)); % the place the

syms f(x)
f(x)=Jiemax*(1-(x*((1-x)*K_Emtilmpsi-(K_Iptilpsi-x*K_Imtilpsi))/((1-x)*K_Eptilmpsi-(K_Iptilpsi-x*K_Imtilpsi)))^(1/mu))-JbarD*((JbarD-2*gabs^2*(1-Jeemean)*(1-x)/(K_Iptilpsi-x*K_Imtilpsi))/(1+Jiimean+2*JbarD*gabs^2*(1-x)/(K_Iptilpsi-x*K_Imtilpsi)));
alpha_cr=double(vpasolve(f,x,[0.1 0.9999]));
J_ie_null_ei_leave=JbarD*((JbarD-2*gabs^2*(1-Jeemean)*K_Ibar/K_Itilpsi)/(1+Jiimean+2*JbarD*gabs^2*K_Ibar/K_Itilpsi)); % departue place of the Jei nullcline off the bifurcation line

plot(JbarD^2/J_ie_final,J_ie_final,'*','Color',[0 0 1], 'HandleVisibility', 'off')
%hold on

Jeimeanarr=[];
Jiemeanarr=[];

m_e_history=0.2;
m_i_history=1.3;

m_e_full=[];
m_i_full=[];
time_full=[];

tmaxSTDP=10*max([tau_mE tau_mI tau_pI tau_pE]); % history time for m for calculating STDP
tSTDP=dt:dt:(tf + tmaxSTDP);
[TSTDP,~]=meshgrid(tSTDP,tSTDP);

Delta=TSTDP-TSTDP';
KpI=1/tau_pI*exp(-Delta*H_I/tau_pI).*heaviside(H_I*Delta);
KmI=1/tau_mI*exp(Delta*H_I/tau_mI).*heaviside(-H_I*Delta);
KpE=1/tau_pE*exp(-Delta*H_E/tau_pE).*heaviside(H_E*Delta);
KmE=1/tau_mE*exp(Delta*H_E/tau_mE).*heaviside(-H_E*Delta);

KI=KpI-alpha*KmI;

%%%% Bifurcation lines etc... %%%%
%plot(JbarD^2./(0:0.1:20),0:0.1:20,'Color','Black','LineWidth',3)
%hold on
%plot((1+Jiimean)*ones(1,length(0:20)),0:20,'Color','Black','LineWidth',3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Transtime=[];    

for i=2:length(t)
    
    %%% construcing the eigenvectors on each time step %%%
    Jeimean=mean(Jei(:));
    Jiemean=mean(Jie(:));
    
    f=(1-Jiemean/Jiemax)^mu;
    

    if i==4
        plot(Jeimean,Jiemean,'.','Color',[0.9 0.1 0.4])
        hold on
    elseif i>4
        h=plot(Jeimean,Jiemean,'.','Color',[0.9 0.1 0.4], 'HandleVisibility', 'off'); % No legend entry
    end    
    Jeimeanarr(i)=Jeimean;
    Jiemeanarr(i)=Jiemean;
    
    plot(JbarD^2/J_ie_final,J_ie_final,'*','Color',[0 0 1], 'HandleVisibility', 'off')
    
    xlim([0 2])
    ylim([0 20])
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    grid on
    
    drawnow
    
    [m_e_now,m_i_now,T_mean_m_e,~,tnetwork]=proj.common.Two_populations_full_rate_model_history(m_e_history,m_i_history,Jee,Jei,Jie,Jii,dt,tf);
    
    m_e_now=m_e_now(1:end-1);
    m_i_now=m_i_now(1:end-1);
    
    
    m_e_history=m_e_now(end-round(D/dt):end); % that is the history for the m_e and m_i for the calculation of the network dynamics
    m_i_history=m_i_now(end-round(D/dt):end); % that is the history for the m_e and m_i for the calculation of the network dynamics
    
    KE=f*KpE-alpha*KmE;
    
    if i==2
        time_full=tnetwork;
        m_e_bef=m_e_now(1)*ones(1,length(0:dt:tmaxSTDP-dt));
        m_i_bef=m_i_now(1)*ones(1,length(0:dt:tmaxSTDP-dt));
    else
        time_full=[time_full time_full(end)+tnetwork];
        m_e_bef=m_e_full((end-tmaxSTDP/dt+1):end);
        m_i_bef=m_i_full((end-tmaxSTDP/dt+1):end);
    end
    
    Transtime(i)=time_full(end);
    
    Jeidot=lambda_i*dt*( [m_i_bef m_i_now]*KI*[m_e_bef m_e_now]' - [m_i_bef 0*m_i_now]*KI*[m_e_bef 0*m_e_now]' )/tf;
    Jiedot=lambda_e*dt*( [m_e_bef m_e_now]*KE*[m_i_bef m_i_now]' - [m_e_bef 0*m_e_now]*KE*[m_i_bef 0*m_i_now]' )/tf;
    
    
    m_e_full=[m_e_full m_e_now];
    m_i_full=[m_i_full m_i_now];
    
    Jei=Jei+dtlearn*Jeidot;
    Jie=Jie+dtlearn*Jiedot;
    
end
%% Plot the network activity of last rounds in Critical Rhythmogenesis %%
figure;
plot(dt*(1:1:length(m_e_full))*Tunits,m_e_full,'Color',	[0.90, 0.40, 0.35]) % in units of s
hold on
plot(dt*(1:1:length(m_i_full))*Tunits,m_i_full,'Color',	[0.30, 0.60, 0.90]) % in units of s
xlabel('$\mathrm{Time} \ [s]$','interpreter','latex','FontSize',18)
ylabel('$m_X(t)$','interpreter','latex','FontSize',18)
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
grid on

for i=1:50
plot(ones(1,10)*Transtime(end-i)*Tunits,linspace(0,max(m_i_full(:)),10),'--','Color','Black')
hold on
end

xlim(Tunits*[Transtime(end-10)-1 Transtime(end-4)+1])
ylim([0 2.5])
lgd=legend({'$\mathrm{E}$','$\mathrm{I}$'},'Interpreter','latex','Location','NorthEast');
legendJiiTitle=sprintf('$X$','Interpreter','latex');
lgd.Title.String = legendJiiTitle;