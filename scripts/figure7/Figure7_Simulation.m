%% Figure 7 Simulation %%
% This script generates **Figure 7a-7d** from the paper:
% "Functional stability and recurrent STDP in rhythmogenesis"
%
% Description:

%       - Computes the STDP dynamics parallely (parfor) for the neural,
%       when at a certain time 3 of the neurons are killed (both for the
%       abscence and presence of intra-connections).
%       - Saves the dynamics in subfolder "data".
%       - Plots the synaptic dynamics with insets to the dyanmics close to
%       the "death time" (no intra - figure 7a, with intra - figure 7c).
%       - Plots the firing rates of the neurons right berfore the death of
%       the neruons, and after the death of the neruons (no intra - figure
%       7b, figure 7d).
%       - Computes the freuqnecy (with SD) of the network in critical rhythmogenesis
%       both before and after the death of the neurons


% Dependencies:

%       - Correlations_2D_full_diff.m - Computes the cross-correlations
%       - FullSynDynMultiTrailsPerturbAndFreq - Computes the STDP dynamics
%       of the network per trial

% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Definitions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Running on abscence (p=0) and presence (p=1) %%%%%
for p=0:1 % p=0 for J_{XX}=0 otherwise JEE=0.6 and JII=0.4
    if ~p % not intra
        filename='data\DynamicswithPerturbation3deadNoIntraBETTER.mat'; % filename for no-intra
        Jiimean=0; % order parmeter of Jii
        Jeemean=0; % order parmeter of Jee
        
    else % with intra
        filename='data\DynamicswithPerturbation3deadwithIntraBETTER.mat'; % filename for with intra
        Jiimean=0.4; % order parmeter of Jii
        Jeemean=0.6; % order parmeter of Jee
        
    end
    
    %%%%%%%%%%%%%%%%% Basic Parameters %%%%%%%%%%%%%%%%%%%%%
    oldparam = sympref('HeavisideAtOrigin',1/2); % heaviside is 0.5 on zero
    dt=0.01; % time bin
    tf=200; % final time of simulation for network dynamics
    %%%%%%%%% Phase diagram features (bif. etc.) %%%%%%%
    Tunits=5*10^-3; % 5 msec is T=1
    T=1; % time constant 5msec tau
    D=0.4; % delay 
    Jhat=(Jeemean+Jiimean)/2;
    syms wD JbarD
    if Jeemean>=Jiimean
        range=[0.1 5 ;0.01 pi/(2*D)];
        Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), T*wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
    elseif Jeemean<Jiimean
        range=[0.1 5 ;0.01 pi/D];
        Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), T*wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
    end
    wD=double(Y.wD); % angular frequency on bif. line
    fD=wD/(2*pi); % frequency on bif. line
    JbarD=double(Y.JbarD); % Jbar on bif. line
    phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
    psi=acos(Jhat/JbarD);
    %%%%%%%%% STDP Parameters %%%%%%%%%
    alpha=0.98; % relative depression
    mu=0.01; % measure of linearity
    Jiemax=20; % J_ie_max
    tau_pE=2; % typical potentiation time of excitatory synapses
    tau_pI=2; % typical potentiation time of inhibitory synapses
    tau_mE=7; % typical depression time of excitatory synapses
    tau_mI=3; % typical depression time of inhibitory synapses
    thetapI=acos((1+(wD*tau_pI)^2)^-0.5);
    thetamI=acos((1+(wD*tau_mI)^2)^-0.5);
    thetapE=acos((1+(wD*tau_pE)^2)^-0.5);
    thetamE=acos((1+(wD*tau_mE)^2)^-0.5);
    lambda_e=10; % excitatory synapses learning rate
    lambda_i=1; % inhibitory synapses learning rate
    H_E=-1;
    H_I=1;
    
    %%%%%%%%% Generating the synapses matrix %%%%%%%%%
    N_e=40; % Size of excitatory population
    N_i=10; % Size of inhibitory population
    
    %%% Starting point for inter-synapses %%%
    if ~p % no intra
        Jiemean=7.37; % Jie order parameter
        Jeimean=0.43; % Jei order parameter
    else % with intra
        Jiemean=8.14;  % Jie order parameter
        Jeimean=0.39; % Jei order parameter
    end
    f=(1-Jiemean/Jiemax)^mu; % synaptic dependent function
    ftag=-mu/Jiemax*(1-Jiemean/Jiemax)^(mu-1); % f derivative by Jie
    %%% Analytical expressions from paper (p stands for plus and m for minus)
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
    
    K_Eptilpsi=cos(thetapE)*cos(thetapE+psi);
    K_Emtilpsi=cos(thetamE)*cos(thetamE-psi);
    K_Etilpsi=f*K_Eptilpsi-alpha*K_Emtilpsi;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Building the vector around the zero order %%%%%%%%%
    meb=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean); % m_e bar
    mib=(1-Jeemean+Jiemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean); % m_i bar
    gabs=(JbarD^2-Jeemean*Jiimean)^0.5; % absolute value of g
    %%%% initialy create the coeffieicnts of the fluctuations (see paper) %%%%
    if Jeimean*Jiemean<JbarD^2 % FP Regime
        a=-lambda_i*K_Ibar*mib^2;
        b=lambda_i*K_Ibar*meb^2;
        c=-lambda_e*K_Ebar*mib^2;
        d=lambda_e*K_Ebar*meb^2;
        q=lambda_e*ftag*meb*mib;
    else % R regime
        mit=(Jiemean*(1+Jiimean)-JbarD^2)/(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2));
        met=JbarD/Jiemean*mit;
        a=-lambda_i*(K_Ibar*mib^2+K_Itilmphi*mit^2/(2*gabs));
        b=lambda_i*(K_Ibar*meb^2+K_Itilphi*met^2/(2*gabs));
        c=-lambda_e*(K_Ebar*mib^2+K_Etilphi*mit^2/(2*gabs));
        d=lambda_e*(K_Ebar*meb^2+K_Etilmphi*met^2/(2*gabs));
        q=lambda_e*ftag*(meb*mib+met*mit/2*K_Eptilpsi); % self deperssion coefficient
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Make the inital connectivity matrix %%%%%%%%%%%%%
    J_ei_std=Jeimean/10; % SD of Jei
    J_ei_noise=normrnd(0,J_ei_std,N_e,N_i); % noise part of Jei
    Jei=(Jeimean*ones(N_e,N_i))-mean(J_ei_noise(:))+J_ei_noise; % Jei connectivity matrix
    
    J_ie_std=Jiemean/10; % SD of Jie
    J_ie_noise=normrnd(0,J_ie_std,N_i,N_e); % noise part of Jie
    Jie=Jiemean*ones(N_i,N_e)-mean(J_ie_noise(:))+J_ie_noise; % Jie connectivity matrix
    
    Jee=Jeemean*ones(N_e,N_e); % make the intra a scalar
    Jii=Jiimean*ones(N_i,N_i);
    %%% Analytical solution %%%
    J=[[Jee/N_e -Jei/N_i];[Jie/N_e -Jii/N_i]]; % connectivity matrix
    diagonal=eye(N_e+N_i);
    mfa=(diagonal-J)\ones(N_e+N_i,1); % The exact fixed point solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Create temporal kernels %%%%%%%
    m_e_history=0.2;
    m_i_history=1.3;
    [Corr_ie,Corr_ei,~,~,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended]=Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf); % create cross correlations and Delta_extended which is the Delta for the temporal kernels
    K_pE=1/tau_pE*exp(-Delta_extended*H_E/tau_pE).*heaviside(H_E*Delta_extended);
    K_mE=1/tau_mE*exp(Delta_extended*H_E/tau_mE).*heaviside(-H_E*Delta_extended);
    K_pI=1/tau_pI*exp(-Delta_extended*H_I/tau_pI).*heaviside(H_I*Delta_extended);
    K_mI=1/tau_mI*exp(Delta_extended*H_I/tau_mI).*heaviside(-H_I*Delta_extended);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    t=1:8000; % time of STDP simulation
    numtrails=1; % number of trials
    num_dead=3; % number of dead inhibitory neurons on death_time
    death_time=3000;
    dtlearn=1; % timebin for STDP dynamics
    
    mulvn=N_e*(N_i-1); % size of family I
    mulve=N_i*(N_e-1); % size of family II
    mulva=N_e-1; % III
    mulvd=N_i-1; % IV
    
    %%%%%%%%%%%%%%%%%%% Create the eigenvectors %%%%%%%%%%%%%%%%%%%%
    % Family I
    vn=[ones(N_i-1,1) -diag(ones(1,N_i-1))];
    zvn=mat2cell(repmat(vn,1,N_e),N_i-1,N_i*ones(1,N_e));
    vn=[blkdiag(zvn{:}) zeros(size((blkdiag(zvn{:})),1),2*N_e*N_i-size((blkdiag(zvn{:})),2))]'/2^0.5;
    % Family II
    ve=[ones(N_e-1,1) -diag(ones(1,N_e-1))];
    zve=mat2cell(repmat(ve,1,N_i),N_e-1,N_e*ones(1,N_i));
    ve=[zeros(size((blkdiag(zve{:})),1),2*N_e*N_i-size((blkdiag(zve{:})),2)) blkdiag(zve{:})]'/2^0.5;
    % Family III
    vaie=repmat([ones(N_e-1,1) -diag(ones(1,N_e-1))],1,N_i);
    zvaei=mat2cell(repmat(ones(1,N_i),1,N_e-1),1,N_i*ones(1,N_e-1));
    vaei=[ones(N_e-1,N_i) -blkdiag(zvaei{:})];
    va=[(a-q)*vaei c*vaie]'/(((a-q)^2+c^2)*2*N_i)^0.5;
    % Family IV
    vdei=repmat([ones(N_i-1,1) -diag(ones(1,N_i-1))],1,N_e);
    zvdie=mat2cell(repmat(ones(1,N_e),1,N_i-1),1,N_e*ones(1,N_i-1));
    vdie=[ones(N_i-1,N_e) -blkdiag(zvdie{:})];
    vd=[b*vdei (d+q)*vdie]'/((b^2+(d+q)^2)*2*N_e)^0.5;
    % Uniform eigenvalues
    vu=[ones(1,2*N_e*N_i) ; ones(1,N_e*N_i) -ones(1,N_e*N_i)]'/(2*N_e*N_i)^0.5;
    % Computing the STDP dynamics
    [stdmebar,stdmibar,stdmetil,stdmitil,Jeimeandyn,Jiemeandyn,frequency,m_e_bef_per,m_i_bef_per,m_e_aft_per,m_i_aft_per,timem,JbarDyesint,Jeibef,Jiebef,Jeiaft,Jieaft,Jeiper,Jieper] = FullSynDynMultiTrailsPerturbAndFreq(m_e_history,m_i_history,D,Jei,Jie,Jee,Jii,t,Jiemax,mu,alpha,JbarD,Jeemean,Jiimean,vn,ve,vu,vaei,vaie,vdei,vdie,mulvn,mulve,mulva,mulvd,dtlearn,lambda_e,lambda_i,K_Ibar,K_Itilmphi,K_Itilphi,K_Eptilpsi,K_Etilphi,K_Etilmphi,K_pI,K_pE,K_mI,K_mE,N_e,N_i,dt,tf,num_dead,death_time);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the mean and SD of frequency of the network 
    f1=frequency(1:death_time)/Tunits; % frequencies before death
    f2=frequency(death_time:end)/Tunits; % frequencies after death
    
    f1mean = mean(f1(~isnan(f1))); % mean freuqnecy before death
    f1error = std(f1(~isnan(f1))); % SD of freuqnecy before death
    f2mean = mean(f2(~isnan(f2))); % mean freuqnecy after death
    f2error = std(f2(~isnan(f2))); % SD of freuqnecy after death
    
    % save data
    save(filename,"stdmebar","stdmibar","stdmetil","stdmitil","Jeimeandyn","Jiemeandyn","frequency","m_e_bef_per","m_i_bef_per","m_e_aft_per","m_i_aft_per","timem","JbarDyesint","Jeibef","Jiebef","Jeiaft","Jieaft","Jeiper","Jieper","f1mean","f2mean","f1error","f2error");
    
    %% Plot all subfigures of figure 7 %%
    figure;
    plot(JbarD^2./(0:0.01:20),0:0.01:20,'Color','Black','LineWidth',3,'LineStyle','-- ')
    hold on
    plot(JbarDyesint^2./(0:0.01:20),0:0.01:20,'Color','Black','LineWidth',3)
    plot((1+Jiimean)*ones(1,length(0:20)),0:20,'Color','Black','LineWidth',3)
    plot((1+Jiimean*(N_i/(N_i-num_dead)))*ones(1,length(0:20)),0:20,'--','Color','Black','LineWidth',3)
    %plot(JbarD^2/J_ie_final,J_ie_final,'*','Color',[0 0 1])
    plot(Jeimeandyn,Jiemeandyn,'.','Color',colo)
    plot(Jeimeandyn(1),Jiemeandyn(1),'+','Color',[1, 0, 0],'MarkerSize', 8, 'LineWidth', 2);
    plot(Jeimeandyn(end),Jiemeandyn(end),'x','Color',[0, 1, 0],'MarkerSize', 8, 'LineWidth', 2);
    plot(Jeimeandyn(death_time),Jiemeandyn(death_time),'.','Color','cyan','Markersize',8)
    plot(Jeibef,Jiebef,'.','Color',[1 0.5 0],'Markersize',8)
    plot(Jeiaft,Jieaft,'.','Color',[1, 1, 0.13],'Markersize',8)
    xlim([0 1+Jiimean*(N_i/(N_i-num_dead))+0.5])
    ylim([0 15])
    xlabel('$J_{EI}$','interpreter','latex','FontSize',18)
    ylabel('$J_{IE}$','interpreter','latex','FontSize',18)
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    grid on
    
    % === Inset axes ===
    inset_pos = [0.6 0.6 0.25 0.25];  % [x y width height] in normalized units
    inset_axes = axes('Position', inset_pos);
    
    % === Content of the inset (completely independent) ===
    t = 1:length(Jeimeandyn);
    plot(death_time*ones(1,length(0:5)),0:5,'--','Color','cyan','LineWidth',1)
    hold on
    y = (Jeimeandyn.*Jiemeandyn).^0.5;
    plot(inset_axes, t, y, '.','Color',colo)
    box(inset_axes, 'on')
    xlabel(inset_axes, '$t \ [\mathrm{a.u.}]$', 'Interpreter', 'latex', 'FontSize', 10)
    ylabel(inset_axes, '$\bar{J}$', 'Interpreter', 'latex', 'FontSize', 10)
    xlim([death_time-100 death_time+100])
    ylim([min(y)-0.1 max(y)+0.1])
    grid on
    
    
    limbot=0.95;
    limtop=1;
    
    figure;
    subplot(2,2,1)
    plot(timem*Tunits,m_e_bef_per,'Color',[0.90, 0.40, 0.35])
    m_e_bef_bar_mean=mean((max(m_e_bef_per,[],2)+min(m_e_bef_per,[],2))/2);
    m_e_bef_bar_std=std((max(m_e_bef_per,[],2)+min(m_e_bef_per,[],2))/2);
    m_e_bef_tilde_mean=mean((max(m_e_bef_per,[],2)-min(m_e_bef_per,[],2))/2);
    m_e_bef_tilde_std=std((max(m_e_bef_per,[],2)-min(m_e_bef_per,[],2))/2);
    % Place a text box in axes coordinates
    % Format text in LaTeX
    str = sprintf('\m_{E,i}: $%.2f$\\\\\\textbf{Std}: $%.2f$', m_e_bef_bar_mean, m_e_bef_bar_std);
    
    % Add the box with LaTeX rendering
    text(0.05, 0.95, str, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'Interpreter', 'latex', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k')
    xlim([limbot limtop])
    ylabel('$m_{E,i} \ \mathrm{before} \ t_{0}$','interpreter','latex','FontSize',18)
    set(gca, 'XTickLabel', []);
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    grid on
    
    subplot(2,2,2)
    plot(timem*Tunits,m_e_aft_per,'Color',[0.90, 0.40, 0.35])
    m_e_aft_bar_mean=mean((max(m_e_aft_per,[],2)+min(m_e_aft_per,[],2))/2);
    m_e_aft_bar_std=std((max(m_e_aft_per,[],2)+min(m_e_aft_per,[],2))/2);
    m_e_aft_tilde_mean=mean((max(m_e_aft_per,[],2)-min(m_e_aft_per,[],2))/2);
    m_e_aft_tilde_std=std((max(m_e_aft_per,[],2)-min(m_e_aft_per,[],2))/2);
    xlim([limbot limtop])
    ylabel('$m_{E,i} \ \mathrm{after} \ t_{0}$','interpreter','latex','FontSize',18)
    set(gca, 'XTickLabel', []);
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    
    grid on
    
    subplot(2,2,3)
    plot(timem*Tunits,m_i_bef_per,'Color',[0.30, 0.60, 0.90])
    m_i_bef_bar_mean=mean((max(m_i_bef_per,[],2)+min(m_i_bef_per,[],2))/2);
    m_i_bef_bar_std=std((max(m_i_bef_per,[],2)+min(m_i_bef_per,[],2))/2);
    m_i_bef_tilde_mean=mean((max(m_i_bef_per,[],2)-min(m_i_bef_per,[],2))/2);
    m_i_bef_tilde_std=std((max(m_i_bef_per,[],2)-min(m_i_bef_per,[],2))/2);
    xlim([limbot limtop])
    ylabel('$m_{I,i} \ \mathrm{before} \ t_{0}$','interpreter','latex','FontSize',18)
    xlabel('$t \ [\mathrm{sec}]$','interpreter','latex','FontSize',18)
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    
    
    grid on
    subplot(2,2,4)
    plot(timem*Tunits,m_i_aft_per((num_dead+1):end,:),'Color',[0.30, 0.60, 0.90])
    m_i_aft_bar_mean=mean((max(m_i_aft_per,[],2)+min(m_i_aft_per,[],2))/2);
    m_i_aft_bar_std=std((max(m_i_aft_per,[],2)+min(m_i_aft_per,[],2))/2);
    m_i_aft_tilde_mean=mean((max(m_i_aft_per,[],2)-min(m_i_aft_per,[],2))/2);
    m_i_aft_tilde_std=std((max(m_i_aft_per,[],2)-min(m_i_aft_per,[],2))/2);
    xlim([limbot limtop])
    ylabel('$m_{I,i} \ \mathrm{after} \ t_{0}$','interpreter','latex','FontSize',18)
    xlabel('$t \ [\mathrm{sec}]$','interpreter','latex','FontSize',18)
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    
    grid on
    
    
    % Reduce space manually
    ax1 = subplot(2,2,1);
    ax2 = subplot(2,2,3);
    
    
    pos1 = get(ax1, 'Position');
    pos2 = get(ax2, 'Position');
    
    % Move ax2 upward
    pos2(2) = pos2(2) + 0.1;
    set(ax2, 'Position', pos2);
    
    % Reduce space manually
    ax1 = subplot(2,2,2);
    ax2 = subplot(2,2,4);
    
    
    pos1 = get(ax1, 'Position');
    pos2 = get(ax2, 'Position');
    
    % Move ax2 upward
    pos2(2) = pos2(2) + 0.1;
    set(ax2, 'Position', pos2);
    
    
end


