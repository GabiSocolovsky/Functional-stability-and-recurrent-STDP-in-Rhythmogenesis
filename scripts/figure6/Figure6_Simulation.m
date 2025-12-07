%% Figure 6 Simulation %%
% This script generates **Figure 6a-6d** from the paper:
% "Functional stability and recurrent STDP in rhythmogenesis"
%
% Description:

%       - Computes the STDP dynamics parallely (parfor) of the fluctuations in each family (fig.
%       6a)
%       - Saves the fluctuation in subfolder "data"
%       - Computes the correlation coefficient in family IV (fig. 6b)
%       - Computes the standar deviation across a population in m_{X,i} bar
%       and m_{X,i} tilde (fig. 6c and fig. 6d)


% Dependencies:

%       - Correlations_2D_full_diff.m - Computes the cross-correlations
%       - FullSynDynMultiTrails.m - computes the STDP dynamics of the full
%       network
%       - parsave - Saves parallely the dynamics of different trials

% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Definitions %%%%%%%%%%%%%
oldparam = sympref("HeavisideAtOrigin",1/2); % 0.5 at the origin for heaviside
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runs the STDP dynamics of a network %%
%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%
oldparam = sympref('HeavisideAtOrigin',1/2); % heaviside is 0.5 on zero
dt=0.01; % time bin
tf=200; % final time of simulation for network dynamics
%%%%%%%%% Phase diagram features (bif. etc.) %%%%%%%
T=1; % time constant 5msec tau
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
wD=double(Y.wD); % angular frequency on the bifurcation line
fD=wD/(2*pi); % frequency on the bifurcation line
JbarD=double(Y.JbarD); % Jbar on bif. line
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% STDP Parameters %%%%%%%%%
alpha=0.98; % relative depression
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
lambda_e=10; % excitatory synapses learning rate
lambda_i=1; % inhibitory synapses learning rate
H_E=-1;
H_I=1;
%%%%%%%%% Generating the synapses matrix %%%%%%%%%
N_e=40; % Size of excitatory population
N_i=10; % Size of inhibitory population
Jiemean=6;
Jeimean=0.52;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    q=lambda_e*ftag*meb*mib;  % self deperssion coefficient
else % R regime
    mit=(Jiemean*(1+Jiimean)-JbarD^2)/(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2));
    met=JbarD/Jiemean*mit; % m_e_tilde
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
[Corr_ie,Corr_ei,~,~,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended]=proj.common.Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf); % create cross correlations and Delta_extended which is the Delta for the temporal kernels
K_pE=1/tau_pE*exp(-Delta_extended*H_E/tau_pE).*heaviside(H_E*Delta_extended);
K_mE=1/tau_mE*exp(Delta_extended*H_E/tau_mE).*heaviside(-H_E*Delta_extended);
K_pI=1/tau_pI*exp(-Delta_extended*H_I/tau_pI).*heaviside(H_I*Delta_extended);
K_mI=1/tau_mI*exp(Delta_extended*H_I/tau_mI).*heaviside(-H_I*Delta_extended);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
t=1:14000; % time of STDP simulation
numtrails=100; % number of trials
check_arr=[reshape((Jei-mean(Jei(:))).',1,[]) reshape((Jie-mean(Jie(:))).',1,[])]; % reorganize the connectivity matrix in a vector of first Jei fluctuations and below Jie fluctuations
J_ei_dot=zeros(N_e,N_i);
J_ie_dot=zeros(N_i,N_e);
dtlearn=1; % timebin for STDP dynamics

mulvn=N_e*(N_i-1); % size of family I
mulve=N_i*(N_e-1); % size of family II
mulva=N_e-1; % size of family III
mulvd=N_i-1; % size of family IV

%%%%%%%%%%%%%%%%%%% Create the eigenvectors %%%%%%%%%%%%%%%%%%%%
% Family I
vn=[ones(N_i-1,1) -diag(ones(1,N_i-1))];
zvn=mat2cell(repmat(vn,1,N_e),N_i-1,N_i*ones(1,N_e));
vn=[blkdiag(zvn{:}) zeros(size((blkdiag(zvn{:})),1),2*N_e*N_i-size((blkdiag(zvn{:})),2))]'/2^0.5; % final
% Family II
ve=[ones(N_e-1,1) -diag(ones(1,N_e-1))];
zve=mat2cell(repmat(ve,1,N_i),N_e-1,N_e*ones(1,N_i));
ve=[zeros(size((blkdiag(zve{:})),1),2*N_e*N_i-size((blkdiag(zve{:})),2)) blkdiag(zve{:})]'/2^0.5; % final
% Family III
vaie=repmat([ones(N_e-1,1) -diag(ones(1,N_e-1))],1,N_i); % final of Jie
zvaei=mat2cell(repmat(ones(1,N_i),1,N_e-1),1,N_i*ones(1,N_e-1));
vaei=[ones(N_e-1,N_i) -blkdiag(zvaei{:})]; % final of Jei
va=[(a-q)*vaei c*vaie]'/(((a-q)^2+c^2)*2*N_i)^0.5; 
% Family IV
vdei=repmat([ones(N_i-1,1) -diag(ones(1,N_i-1))],1,N_e); % final of Jei
zvdie=mat2cell(repmat(ones(1,N_e),1,N_i-1),1,N_e*ones(1,N_i-1));
vdie=[ones(N_i-1,N_e) -blkdiag(zvdie{:})]; % final of Jie
vd=[b*vdei (d+q)*vdie]'/((b^2+(d+q)^2)*2*N_e)^0.5;
% Uniform eigenvalues
vu=[ones(1,2*N_e*N_i) ; ones(1,N_e*N_i) -ones(1,N_e*N_i)]'/(2*N_e*N_i)^0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Run dynamics of different trials parallely %%%%%%%%%%%%%%%%%%%%%%%%%
parfor trail=1:numtrails
    [Vnvol,Vevol,Vavol,Vdvol,stdmebar,stdmibar,stdmetil,stdmitil,Jeimeandyn,Jiemeandyn,dJ_arr_t] = FullSynDynMultiTrails(m_e_history,m_i_history,Jei,Jie,Jee,Jii,t,Jiemax,mu,alpha,JbarD,Jeemean,Jiimean,vn,ve,vu,vaei,vaie,vdei,vdie,mulvn,mulve,mulva,mulvd,dtlearn,lambda_e,lambda_i,K_Ibar,K_Itilmphi,K_Itilphi,K_Eptilpsi,K_Etilphi,K_Etilmphi,K_pI,K_pE,K_mI,K_mE,N_e,N_i,dt,tf) % Run dynamics STDP of a trial
    parsave(sprintf('Figure6_simulationsdataNEW%d.mat',trail),Vnvol,Vevol,Vavol,Vdvol,stdmebar,stdmibar,stdmetil,stdmitil,Jeimeandyn,Jiemeandyn,dJ_arr_t); % save it
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract data from subfolder \data to numbered (by trials) variables %%
base_filename='Figure6_simulationsdataNEW'; % the common name to all trials
numtrails=100; % number of trials
%%%% making variables with all trials %%% 
% Not relevant for figure %
Vnvol_alltrails=[];
Vevol_alltrails=[];
Vavol_alltrails=[];
Vdvol_alltrails=[];
% Relevant for figures %
stdmebar_alltrails=[];
stdmibar_alltrails=[];
stdmetil_alltrails=[];
stdmitil_alltrails=[];
Jeimeandynalltrails=[];
Jiemeandynalltrails=[];
dJarralltrails=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numtrails % Make these variables which have all the information across trials
    changing_index=i;
    full_filename = strcat(base_filename, num2str(changing_index), '.mat');
    data = load(full_filename);
    Vnvol_alltrails(i,:)=data.DynVar.Vnvol;
    Vevol_alltrails(i,:)=data.DynVar.Vevol;
    Vavol_alltrails(i,:)=data.DynVar.Vavol;
    Vdvol_alltrails(i,:)=data.DynVar.Vdvol;
    stdmebar_alltrails(i,:)=data.DynVar.stdmebar;
    stdmibar_alltrails(i,:)=data.DynVar.stdmibar;
    stdmetil_alltrails(i,:)=data.DynVar.stdmetil;
    stdmitil_alltrails(i,:)=data.DynVar.stdmitil;
    Jeimeandynalltrails(i,:)=data.DynVar.Jeimeandyn;
    Jiemeandynalltrails(i,:)=data.DynVar.Jiemeandyn;
    dJarr=data.DynVar.dJ_arr_t;
    dJarralltrails(i,:,:)=dJarr';
end
%% Compute the variance of each family (not normalized yet - see section below for normalization) %% 
% Important! - the parfor in the second loop is valid because we are not running the dynamics here, we are just taking the data of the already computed dynamics and use this data to compute the variances
% save in subfolder \data the variances in variable "VolsqFULL"
oldparam = sympref('HeavisideAtOrigin',1/2); % heaviside is 0.5 on zero
dt=0.01; % time bin
tf=200; % final time of simulation for network dynamics
%%%%%%%%% Phase diagram features (bif. etc.) %%%%%%%
T=1; % time constant 5msec tau
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
wD=double(Y.wD); % angular frequency on the bifurcation line
fD=wD/(2*pi); % frequency on the bifurcation line
JbarD=double(Y.JbarD);
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% STDP Parameters %%%%%%%%%
alpha=0.98; % relative depression
mu=0.01; % measure of linearity
Jiemax=20; % J_ie_max
tau_pE=2; % typical potentiation time of excitatory synapses
tau_pI=2; % typical potentiation time of inhibitory synapses
tau_mE=5; % typical depression time of excitatory synapses
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

K_Eptilphi=cos(thetapE)*cos(thetapE-phi);
K_Emtilphi=cos(thetamE)*cos(thetamE+phi);
K_Eptilmphi=cos(thetapE)*cos(thetapE+phi);
K_Emtilmphi=cos(thetamE)*cos(thetamE-phi);

K_Eptilpsi=cos(thetapE)*cos(thetapE+psi);
K_Emtilpsi=cos(thetamE)*cos(thetamE-psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_ie_final=Jiemax*(1-((1-K_Ibar)*(K_Ibar*K_Emtilpsi-K_Itilpsi)/(K_Ibar*K_Eptilpsi-K_Itilpsi))^(1/mu)); % Jie asterisk

numtrails=100; % number of trials
tflearn=14000; % final time for STDP simulation

%%%%%%  Variances of families (not normalized yet) %%%%%%%
Volsq1=zeros(numtrails,tflearn); % I
Volsq2=zeros(numtrails,tflearn); % II
Volsq3=zeros(numtrails,tflearn); % III
Volsq4=zeros(numtrails,tflearn); % IV
Volsq5=zeros(numtrails,tflearn); % Uniform eigenvectors (they are two)

%%%%% computing the variances of each family in each trial at each time %%%%%
for i=1:numtrails % run on different trials
    parfor j=1:tflearn % use parfor to compute the variances at each time (we can use parfor because we are not computing the dynamics here, it is already computed)
        
        Jiemean=Jiemeandynalltrails(i,j);
        Jeimean=Jeimeandynalltrails(i,j);
        
        f=(1-Jiemean/Jiemax)^mu;
        ftag=-mu/Jiemax*(1-Jiemean/Jiemax)^(mu-1);
        
        K_Ebar=f-alpha;
        K_Etilpsi=f*K_Eptilpsi-alpha*K_Emtilpsi;
        K_Etilmphi=f*K_Eptilmphi-alpha*K_Emtilmphi;
        K_Etilphi=f*K_Eptilphi-alpha*K_Emtilphi;
        %%%%%%%% Building the vector around the zero order %%%%%%%%%
        meb=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean); % me bar
        mib=(1-Jeemean+Jiemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean); % mi bar
        
        if Jeimean*Jiemean<JbarD^2
            a=-lambda_i*K_Ibar*mib^2;
            b=lambda_i*K_Ibar*meb^2;
            c=-lambda_e*K_Ebar*mib^2;
            d=lambda_e*K_Ebar*meb^2;
            q=lambda_e*ftag*meb*mib;
        else
            gabs=(JbarD^2-Jeemean*Jiimean)^0.5; % absolute value of g 
            mit=(JbarD^2/Jeimean*(1+Jiimean)-JbarD^2)/(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2)); % mi tilde
            met=JbarD/Jiemean*mit; % me tilde
            a=-lambda_i*(K_Ibar*mib^2+K_Itilmphi*mit^2/(2*gabs));
            b=lambda_i*(K_Ibar*meb^2+K_Itilphi*met^2/(2*gabs));
            c=-lambda_e*(K_Ebar*mib^2+K_Etilphi*mit^2/(2*gabs));
            d=lambda_e*(K_Ebar*meb^2+K_Etilmphi*met^2/(2*gabs));
            q=lambda_e*ftag*(meb*mib+met*mit/2*K_Eptilpsi); % self deperssion coefficient 
        end
        % create the vectors
        vn=[ones(N_i-1,1) -diag(ones(1,N_i-1))];
        zvn=mat2cell(repmat(vn,1,N_e),N_i-1,N_i*ones(1,N_e));
        vn=[blkdiag(zvn{:}) zeros(size((blkdiag(zvn{:})),1),2*N_e*N_i-size((blkdiag(zvn{:})),2))]'/2^0.5; % eigenvector family I
        
        ve=[ones(N_e-1,1) -diag(ones(1,N_e-1))];
        zve=mat2cell(repmat(ve,1,N_i),N_e-1,N_e*ones(1,N_i));
        ve=[zeros(size((blkdiag(zve{:})),1),2*N_e*N_i-size((blkdiag(zve{:})),2)) blkdiag(zve{:})]'/2^0.5; % eigenvector family II
        
        vaie=repmat([ones(N_e-1,1) -diag(ones(1,N_e-1))],1,N_i);
        zvaei=mat2cell(repmat(ones(1,N_i),1,N_e-1),1,N_i*ones(1,N_e-1));
        vaei=[ones(N_e-1,N_i) -blkdiag(zvaei{:})];
        va=[(a-q)*vaei c*vaie]'/(((a-q)^2+c^2)*2*N_i)^0.5; % eigenvector family III
        
        vdei=repmat([ones(N_i-1,1) -diag(ones(1,N_i-1))],1,N_e);
        zvdie=mat2cell(repmat(ones(1,N_e),1,N_i-1),1,N_e*ones(1,N_i-1));
        vdie=[ones(N_i-1,N_e) -blkdiag(zvdie{:})];
        vd=[b*vdei (d+q)*vdie]'/((b^2+(d+q)^2)*2*N_e)^0.5; % eigenvector family IV
        
        vu=[ones(1,2*N_e*N_i) ; ones(1,N_e*N_i) -ones(1,N_e*N_i)]'/(2*N_e*N_i)^0.5; % eigenvector unifrom family
        
        %%% Graham-Schmidt %%%
        vn=gsog(vn);
        va=gsog(va);
        vd=gsog(vd);
        vu=gsog(vu);
        ve=gsog(ve);
        %%% Normalize the eigenvectors %%%
        vn=vn./repmat(sum(vn.^2,1).^0.5,N_e*N_i*2,1);
        va=va./repmat(sum(va.^2,1).^0.5,N_e*N_i*2,1);
        vd=vd./repmat(sum(vd.^2,1).^0.5,N_e*N_i*2,1);
        vu=vu./repmat(sum(vu.^2,1).^0.5,N_e*N_i*2,1);
        ve=ve./repmat(sum(ve.^2,1).^0.5,N_e*N_i*2,1);
        
        % make them a set %
        V=[vn  ve  va  vd  vu];
        % size of subspaces (families)
        vn_size=size(vn,2);
        ve_size=size(ve,2);
        va_size=size(va,2);
        vd_size=size(vd,2);
        vu_size=size(vu,2);
        
        
        dJ_arr=squeeze(dJarralltrails(i,:,j)); % take the fluctuations of this time
        
        % Compute coefficients of dJ_arr on families % 
        coef=V\dJ_arr';
        coef1=diag([ones(1,vn_size) zeros(1,ve_size) zeros(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef;
        coef2=diag([zeros(1,vn_size) ones(1,ve_size) zeros(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef;
        coef3=diag([zeros(1,vn_size) zeros(1,ve_size) ones(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef;
        coef4=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) ones(1,vd_size) zeros(1,vu_size)])*coef;
        coef5=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) zeros(1,vd_size) ones(1,vu_size)])*coef; % uniform family
        % Variance %
        Volsq1(i,j)=sum(coef1.^2);
        Volsq2(i,j)=sum(coef2.^2);
        Volsq3(i,j)=sum(coef3.^2);
        Volsq4(i,j)=sum(coef4.^2);
        Volsq5(i,j)=sum(coef5.^2);
        
    end
end

save('data/VolsqFULL','Volsq1','Volsq2','Volsq3','Volsq4','Volsq5'); % Save in subfolder \data
%% Plotting the Var of every family (figure 6a) and the SD in each population (figure 6c-6d) %%

figure(1)
stdshade(stdmebar_alltrails,0.3,[0.90, 0.40, 0.35],[],[])
hold on
stdshade(stdmibar_alltrails,0.3,[0.30, 0.60, 0.90],[],[])

ylabel('SD of $\bar{m}_{X,i}$ [a.u.]','interpreter','latex','FontSize',18)
xlabel('$t \ [\mathrm{a.u.}]$','interpreter','latex','FontSize',18)
lgd=legend({'','$\mathrm{E}$','','$\mathrm{I}$'},'Interpreter','latex','Location','Northeast');
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
xlim([0 14000])
grid on

% Plot SD of me and mi tilde %
figure(2);
stdshade_nanfriendly(stdmetil_alltrails, 0.2, [0.90, 0.40, 0.35], 2);
hold on
stdshade_nanfriendly(stdmitil_alltrails, 0.2, [0.30, 0.60, 0.90], 2);%ylim([0 0.0008])
ylabel('SD of $\tilde{m}_{X,i}$ [a.u.]','interpreter','latex','FontSize',18)
xlabel('$t \ [\mathrm{a.u.}]$','interpreter','latex','FontSize',18)
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
lgd=legend({'','$\mathrm{E}$','','$\mathrm{I}$'},'Interpreter','latex','Location','Northeast');
xlim([0 14000])
grid on

% Plot Normalized Variances %
figure(3);
stdshade((Volsq1)/(N_e*(N_i-1)),0.3,[1 0 0],[],[])
hold on
stdshade(Volsq2/(N_i*(N_e-1)),0.1,[0 1 0],[],[])
stdshade(Volsq3/(N_e-1),0.3,'magenta',[],[])
stdshade(Volsq4/(N_i-1),0.3,[0 0 1],[],[])
grid on
xlabel('$t \ [\mathrm{a.u.}]$','interpreter','latex','FontSize',18)
ylabel('$\mathrm{Fluctuations \ variance} \ [\mathrm{a.u.}]$','interpreter','latex','FontSize',18)
%lgd=legend({'$V_{0}$','$V_{E}$','$V_{A}$','$V_{D+E}$'},'Interpreter','latex','Location','Northeast');
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
grid on
xlim([0 14000])
%% Add asymptotic variance of family IV in figure of variances (dashed line in Fig 6a) %%%
oldparam = sympref('HeavisideAtOrigin',1/2); % heaviside is 0.5 on zero
%%%%%%%%% Phase diagram features (bif. etc.) %%%%%%%
T=1; % time constant 5msec tau
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
wD=double(Y.wD); % frequency on bif. line
fD=wD/(2*pi); 
JbarD=double(Y.JbarD); % Jbar on bif. line
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% STDP Parameters %%%%%%%%%%%%%%%%%%%%%
alpha=0.98; % relative depression
mu=0.01; % measure of linearity
Jiemax=20; % J_ie_max
tau_pE=2; % typical potentiation time of excitatory synapses
tau_pI=2; % typical potentiation time of inhibitory synapses
tau_mE=5; % typical depression time of excitatory synapses
tau_mI=3; % typical depression time of inhibitory synapses
thetapI=acos((1+(wD*tau_pI)^2)^-0.5);
thetamI=acos((1+(wD*tau_mI)^2)^-0.5);
thetapE=acos((1+(wD*tau_pE)^2)^-0.5);
thetamE=acos((1+(wD*tau_mE)^2)^-0.5);
lambda_e=10; % excitatory synapses learning rate
lambda_i=1; % inhibitory synapses learning rate
%%%%%%%%% Generating the synapses matrix %%%%%%%%%
N_e=40; % Size of excitatory population
N_i=10; % Size of inhibitory population
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

K_Eptilphi=cos(thetapE)*cos(thetapE-phi);
K_Emtilphi=cos(thetamE)*cos(thetamE+phi);
K_Eptilmphi=cos(thetapE)*cos(thetapE+phi);
K_Emtilmphi=cos(thetamE)*cos(thetamE-phi);

K_Eptilpsi=cos(thetapE)*cos(thetapE+psi);
K_Emtilpsi=cos(thetamE)*cos(thetamE-psi);

J_ie_final=Jiemax*(1-((1-K_Ibar)*(K_Ibar*K_Emtilpsi-K_Itilpsi)/(K_Ibar*K_Eptilpsi-K_Itilpsi))^(1/mu)); % Jie asterisk
Jiemean=J_ie_final; % start from J asterisk (start from asymptotic value)
Jeimean=JbarD^2/J_ie_final;

f=(1-Jiemean/Jiemax)^mu; % synaptic depndent function
ftag=-mu/Jiemax*(1-Jiemean/Jiemax)^(mu-1);

K_Ebar=f-alpha;
K_Etilmphi=f*K_Eptilmphi-alpha*K_Emtilmphi;
K_Etilphi=f*K_Eptilphi-alpha*K_Emtilphi;
%%%%%%%% Building the vector around the zero order %%%%%%%%%
meb=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean); % me bar
mib=(1-Jeemean+Jiemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean); % mi bar

% coeffiecnients in FP region %
aFP=-lambda_i*K_Ibar*mib^2;
bFP=lambda_i*K_Ibar*meb^2;
cFP=-lambda_e*K_Ebar*mib^2;
dFP=lambda_e*K_Ebar*meb^2;
qFP=lambda_e*ftag*meb*mib;

gabs=(JbarD^2-Jeemean*Jiimean)^0.5; % absolute value of g
mit=(JbarD^2/Jeimean*(1+Jiimean)-JbarD^2)/(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2)); % mi tilde
met=JbarD/Jiemean*mit; % me tilde
% coeffiecnients in R region %
aR=-lambda_i*(K_Ibar*mib^2+K_Itilmphi*mit^2/(2*gabs));
bR=lambda_i*(K_Ibar*meb^2+K_Itilphi*met^2/(2*gabs));
cR=-lambda_e*(K_Ebar*mib^2+K_Etilphi*mit^2/(2*gabs));
dR=lambda_e*(K_Ebar*meb^2+K_Etilmphi*met^2/(2*gabs));
qR=lambda_e*ftag*(meb*mib+met*mit/2*K_Eptilpsi); % self deperssion coefficient 
% create the vectors
for p=1:2
    if p==1 % first vectors in FP region
        a=aFP;
        b=bFP;
        c=cFP;
        d=dFP;
        q=qFP;
    else % second vectors in R region
        a=aR;
        b=bR;
        c=cR;
        d=dR;
        q=qR;
    end
    % Eigenvectors of fmailies
    vn=[ones(N_i-1,1) -diag(ones(1,N_i-1))];
    zvn=mat2cell(repmat(vn,1,N_e),N_i-1,N_i*ones(1,N_e));
    vn=[blkdiag(zvn{:}) zeros(size((blkdiag(zvn{:})),1),2*N_e*N_i-size((blkdiag(zvn{:})),2))]'/2^0.5; % family I
    
    ve=[ones(N_e-1,1) -diag(ones(1,N_e-1))];
    zve=mat2cell(repmat(ve,1,N_i),N_e-1,N_e*ones(1,N_i));
    ve=[zeros(size((blkdiag(zve{:})),1),2*N_e*N_i-size((blkdiag(zve{:})),2)) blkdiag(zve{:})]'/2^0.5; % II
    
    vaie=repmat([ones(N_e-1,1) -diag(ones(1,N_e-1))],1,N_i);
    zvaei=mat2cell(repmat(ones(1,N_i),1,N_e-1),1,N_i*ones(1,N_e-1));
    vaei=[ones(N_e-1,N_i) -blkdiag(zvaei{:})];
    va=[(a-q)*vaei c*vaie]'/(((a-q)^2+c^2)*2*N_i)^0.5; % III
    
    vdei=repmat([ones(N_i-1,1) -diag(ones(1,N_i-1))],1,N_e);
    zvdie=mat2cell(repmat(ones(1,N_e),1,N_i-1),1,N_e*ones(1,N_i-1));
    vdie=[ones(N_i-1,N_e) -blkdiag(zvdie{:})];
    vd=[b*vdei (d+q)*vdie]'/((b^2+(d+q)^2)*2*N_e)^0.5; % IV
    
    vu=[ones(1,2*N_e*N_i) ; ones(1,N_e*N_i) -ones(1,N_e*N_i)]'/(2*N_e*N_i)^0.5; % Uniform eigenvectors
    
    % Graham Schmidt %
    vn=gsog(vn);
    va=gsog(va);
    vd=gsog(vd);
    vu=gsog(vu);
    ve=gsog(ve);
    
    % Normalize eigenvectors %
    vn=vn./repmat(sum(vn.^2,1).^0.5,N_e*N_i*2,1);
    va=va./repmat(sum(va.^2,1).^0.5,N_e*N_i*2,1);
    vd=vd./repmat(sum(vd.^2,1).^0.5,N_e*N_i*2,1);
    vu=vu./repmat(sum(vu.^2,1).^0.5,N_e*N_i*2,1);
    ve=ve./repmat(sum(ve.^2,1).^0.5,N_e*N_i*2,1);
    
    % put all of them in a set
    if p==1
        VFP=[vn  ve  va  vd  vu];
    else
        VR=[vn  ve  va  vd  vu];
    end
end
% size of families
vn_size=size(vn,2);
ve_size=size(ve,2);
va_size=size(va,2);
vd_size=size(vd,2);
vu_size=size(vu,2);
% coefficients (anatical)
MEEbar=lambda_e*K_Ebar*meb^2;
MIEbar=lambda_i*K_Ibar*meb^2;
MEEtil=MEEbar+lambda_e*K_Etilmphi*met^2/(2*gabs);
MIEtil=MIEbar+lambda_i*K_Itilmphi*met^2/(2*gabs);
FEbar=lambda_e*ftag*meb*mib;
FEtil=FEbar+lambda_e*ftag*K_Eptilpsi*met*mit/2;

sigFP4=MEEbar+FEbar; % Sigma IV in FP
sigR4=MEEtil+FEtil; % Sigma IV in R 
etaFP4=(MEEbar+FEbar)/MIEbar; % Eta IV in FP
etaR4=(MEEtil+FEtil)/MIEtil; % Eta IV in R
Ta=etaFP4/etaR4*((1+etaR4^2)/(1+etaFP4^2))^0.5; % Coefficient named "a" matrix T in Appendix 6

noiseEIVar=10^-4; % noise in EI synapses
noiseIEVar=10^-2; % noise in IE synapses

VFPproj=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) ones(1,vd_size) zeros(1,vu_size)])/VFP;
Cxibarsq=diag([noiseEIVar*ones(1,N_e*N_i) noiseIEVar*ones(1,N_e*N_i)]);
xibarsqmean=trace(transpose(VFPproj)*VFPproj*Cxibarsq)/(N_i-1); % xi bar mean squared

VRproj=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) ones(1,vd_size) zeros(1,vu_size)])/VR;
Cxitilsq=diag([noiseEIVar*ones(1,N_e*N_i) noiseIEVar*ones(1,N_e*N_i)]);
xitilsqmean=trace(transpose(VRproj)*VRproj*Cxitilsq)/(N_i-1); % xi tilde mean squared


JeidotFP=lambda_i*K_Ibar.*meb.*mib; % Jei dot in the FP region
JeidotR=JeidotFP+lambda_i*K_Itilpsi.*met.*mit/2; % Jei dot in the R region

c4varFP=(xibarsqmean/abs(JeidotFP)+xitilsqmean/(abs(JeidotR)*Ta^2))/(-2*(sigFP4/JeidotFP+sigR4/abs(JeidotR))); % variance of family IV in FP region
c4varR=(xibarsqmean*Ta^2/abs(JeidotFP)+xitilsqmean/(abs(JeidotR)))/(-2*(sigFP4/JeidotFP+sigR4/abs(JeidotR))); % variance of family IV in R region

c4var=(c4varFP/abs(JeidotFP)+c4varR/abs(JeidotR))/(1/abs(JeidotFP)+1/abs(JeidotR)); % asymptotic effective variance 

% plot
figure(3)
plot(1:14000,c4var*ones(1,length(1:14000)),'--','Color','black','LineWidth',2)
hold on

lgd=legend({'','$\mathrm{I}$','','$\mathrm{II}$','','$\mathrm{III}$','','$\mathrm{IV}$','',''},'Interpreter','latex','Location','Northeast');
%% Computing the autocorrelations coefficients (fig 6b)  %%
% From the data we compute the autocorrelations coefficients
% We assume here that we arrived Cri. Rhythm. and compute all the coefficients a,b,c,d,q in FP and R only once
oldparam = sympref('HeavisideAtOrigin',1/2); % heaviside is 0.5 on zero
%%%%%%%%% Phase diagram features (bif. etc.) %%%%%%%
T=1; % time constant 5msec tau
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
wD=double(Y.wD); % angular frequency on the bif. line
fD=wD/(2*pi); % frequency on the bif. line
JbarD=double(Y.JbarD); % Jbar on the bif. line
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% STDP parameters %%%%%%%%%%%%%%
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
lambda_e=10; % excitatory synapses learning rate
lambda_i=1; % inhibitory synapses learning rate
H_E=-1;
H_I=1;
%%%%%%%%% Generating the synapses matrix %%%%%%%%%
N_e=40; % Size of excitatory population
N_i=10; % Size of inhibitory population
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

K_Eptilphi=cos(thetapE)*cos(thetapE-phi);
K_Emtilphi=cos(thetamE)*cos(thetamE+phi);
K_Eptilmphi=cos(thetapE)*cos(thetapE+phi);
K_Emtilmphi=cos(thetamE)*cos(thetamE-phi);

K_Eptilpsi=cos(thetapE)*cos(thetapE+psi);
K_Emtilpsi=cos(thetamE)*cos(thetamE-psi);

J_ie_final=Jiemax*(1-((1-K_Ibar)*(K_Ibar*K_Emtilpsi-K_Itilpsi)/(K_Ibar*K_Eptilpsi-K_Itilpsi))^(1/mu)); % Jie asterisk
Jiemean=J_ie_final; % start from J asterisk (start from asymptotic value)
Jeimean=JbarD^2/Jiemean;

f=(1-Jiemean/Jiemax).^mu; % synaptic depndent function
ftag=-mu/Jiemax*(1-Jiemean/Jiemax)^(mu-1);

K_Ebar=f-alpha;
K_Etilpsi=f*K_Eptilpsi-alpha*K_Emtilpsi;
K_Etilmphi=f*K_Eptilmphi-alpha*K_Emtilmphi;
K_Etilphi=f*K_Eptilphi-alpha*K_Emtilphi;
%%%%%%%% Building the vector around the zero order %%%%%%%%%
meb=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
mib=(1-Jeemean+Jiemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);

% coeffiecnients in FP region %
aFP=-lambda_i*K_Ibar*mib^2;
bFP=lambda_i*K_Ibar*meb^2;
cFP=-lambda_e*K_Ebar*mib^2;
dFP=lambda_e*K_Ebar*meb^2;
qFP=lambda_e*ftag*meb*mib;

gabs=(JbarD^2-Jeemean*Jiimean)^0.5; % absolute value of g
mit=(JbarD^2/Jeimean*(1+Jiimean)-JbarD^2)/(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2)); % mi tilde
met=JbarD/Jiemean*mit; % me tilde
% coeffiecnients in R region %
aR=-lambda_i*(K_Ibar*mib^2+K_Itilmphi*mit^2/(2*gabs));
bR=lambda_i*(K_Ibar*meb^2+K_Itilphi*met^2/(2*gabs));
cR=-lambda_e*(K_Ebar*mib^2+K_Etilphi*mit^2/(2*gabs));
dR=lambda_e*(K_Ebar*meb^2+K_Etilmphi*met^2/(2*gabs));
qR=lambda_e*ftag*(meb*mib+met*mit/2*K_Eptilpsi); % self deperssion coefficient
% create the vectors
for p=1:2
    
    if p==1 % first vectors in FP region
        a=aFP;
        b=bFP;
        c=cFP;
        d=dFP;
        q=qFP;
    else % second vectors in R region
        a=aR;
        b=bR;
        c=cR;
        d=dR;
        q=qR;
    end
    
    % Eigenvectors of families
    vn=[ones(N_i-1,1) -diag(ones(1,N_i-1))];
    zvn=mat2cell(repmat(vn,1,N_e),N_i-1,N_i*ones(1,N_e));
    vn=[blkdiag(zvn{:}) zeros(size((blkdiag(zvn{:})),1),2*N_e*N_i-size((blkdiag(zvn{:})),2))]'/2^0.5; % family I
    
    ve=[ones(N_e-1,1) -diag(ones(1,N_e-1))];
    zve=mat2cell(repmat(ve,1,N_i),N_e-1,N_e*ones(1,N_i));
    ve=[zeros(size((blkdiag(zve{:})),1),2*N_e*N_i-size((blkdiag(zve{:})),2)) blkdiag(zve{:})]'/2^0.5; % family II
    
    vaie=repmat([ones(N_e-1,1) -diag(ones(1,N_e-1))],1,N_i);
    zvaei=mat2cell(repmat(ones(1,N_i),1,N_e-1),1,N_i*ones(1,N_e-1));
    vaei=[ones(N_e-1,N_i) -blkdiag(zvaei{:})];
    va=[(a-q)*vaei c*vaie]'/(((a-q)^2+c^2)*2*N_i)^0.5; % III
    
    vdei=repmat([ones(N_i-1,1) -diag(ones(1,N_i-1))],1,N_e);
    zvdie=mat2cell(repmat(ones(1,N_e),1,N_i-1),1,N_e*ones(1,N_i-1));
    vdie=[ones(N_i-1,N_e) -blkdiag(zvdie{:})];
    vd=[b*vdei (d+q)*vdie]'/((b^2+(d+q)^2)*2*N_e)^0.5; % IV
    
    vu=[ones(1,2*N_e*N_i) ; ones(1,N_e*N_i) -ones(1,N_e*N_i)]'/(2*N_e*N_i)^0.5; % uniform family
    
    % Graham Schmidt %
    vn=gsog(vn);
    va=gsog(va);
    vd=gsog(vd);
    vu=gsog(vu);
    ve=gsog(ve);
    
    % Normalize eigenvectors %    
    vn=vn./repmat(sum(vn.^2,1).^0.5,N_e*N_i*2,1);
    va=va./repmat(sum(va.^2,1).^0.5,N_e*N_i*2,1);
    vd=vd./repmat(sum(vd.^2,1).^0.5,N_e*N_i*2,1);
    vu=vu./repmat(sum(vu.^2,1).^0.5,N_e*N_i*2,1);
    ve=ve./repmat(sum(ve.^2,1).^0.5,N_e*N_i*2,1);
    
    % put all of them in a set    
    if p==1
        VFP=[vn  ve  va  vd  vu];
    else
        VR=[vn  ve  va  vd  vu];
    end
    
end

VR=sign(VR.*VFP).*VR; % We want the eigenvectors in each region to have the same signs

% Sizes of families
vn_size=size(vn,2);
ve_size=size(ve,2);
va_size=size(va,2);
vd_size=size(vd,2);
vu_size=size(vu,2);

numtrails=100; % number of trials
tflearn=14000; % Time of STDP simulation

t0acorr=2000; % inital time for computing autocorrelation
tfacorr=14000; % final time for computing autocorrelation
maxlag=900; % maximum difference in time for autocorr
acorr4=zeros(1,length(0:maxlag)); % autocorrelation (empty)


for i=1:numtrails % computing autocorrelations in each trial
    
    dJ_arr=squeeze(dJarralltrails(i,:,t0acorr:tfacorr))'; % taking synaptic fluctuaitons when system arrives CR state    
    % Fluctuations in FP region %
    yFP=heaviside(JbarD^2-squeeze(Jeimeandynalltrails(i,t0acorr:tfacorr)).*squeeze(Jiemeandynalltrails(i,t0acorr:tfacorr)))'; % making a vector saying when in FP region (1) and when not (0)
    repyFP=repmat(yFP,1,N_e*N_i*2); % replicating this vector to the size of a network
    dJ_arrFP=dJ_arr.*repyFP; % multiplying the fluctuations with repyFP, so only fluctuations in the FP remain
    
    % Fluctuatuions in R region
    yR=heaviside(-JbarD^2+squeeze(Jeimeandynalltrails(i,t0acorr:tfacorr)).*squeeze(Jiemeandynalltrails(i,t0acorr:tfacorr)))'; % making a vector saying when in FP region (1) and when not (0)
    repyR=repmat(yR,1,N_e*N_i*2); % replicating this vector to the size of a network
    dJ_arrR=dJ_arr.*repyR; % multiplying the fluctuations with repyR, so only fluctuations in the R remain
    
    % Find coefficients in each family
    coef=VFP\dJ_arrFP'+VR\dJ_arrR'; % all coefficients
    coef1=diag([ones(1,vn_size) zeros(1,ve_size) zeros(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef; % only coeff. in family I appear
    coef2=diag([zeros(1,vn_size) ones(1,ve_size) zeros(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef; % II
    coef3=diag([zeros(1,vn_size) zeros(1,ve_size) ones(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef; % III
    coef4=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) ones(1,vd_size) zeros(1,vu_size)])*coef; % IV
    coef5=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) zeros(1,vd_size) ones(1,vu_size)])*coef; % uniform family
    
    % Computing autocorrrelation
    coef4st=coef4(:,1:end-1000); % take the coefficients to up to a certain point in time
    Tint=size(coef4st,2); % the overall time
    
    parfor j=0:maxlag % parallely compute the autocorrelations
        acorr4parfor(1+j)=trace(coef4st*coef4(:,(1+j):(j+Tint))')/(Tint*vd_size);
    end
    
    acorr4(i,:)=acorr4parfor; % put the autocorrelation for a given tau (lag) in a certain row
end

save('data/AutoCorrV4','acorr4'); % save
%% Plot the autocorrelation coefficient (numerical and analytical - figure 6b) %%
figure(4)

stdshade(acorr4./acorr4(:,1),0.3,[0.90, 0.40, 0.35],[],[]) % autocorrelation coefficient numerically
hold on
acorr0=c4var;
tau=0:900; % lag time
sig4=(sigFP4/JeidotFP+sigR4/abs(JeidotR))/(1/JeidotFP+1/abs(JeidotR));
typtime=1/abs(sig4);
plot(tau,exp(tau*(sigFP4/JeidotFP+sigR4/abs(JeidotR))/(1/JeidotFP+1/abs(JeidotR))),'--','Color','black') % plot autocorrelation coefficient analytically
xlabel('$\tau \ [\mathrm{a.u.}]$','interpreter','latex','FontSize',18)
ylabel('$\mathrm{\mathcal{\rho}_{IV}(\tau)}$','interpreter','latex','FontSize',18)
%lgd=legend({'$V_{0}$','$V_{E}$','$V_{A}$','$V_{D+E}$'},'Interpreter','latex','Location','Northeast');
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
grid on
xlim([0 750])

grid on