%% The full synaptic dynamics parfor %%
oldparam = sympref('HeavisideAtOrigin',1/2); % heaviside is 0.5 on zero
dt=0.01; % the accuracy through all the script and auxiliary functions (except for the nullclines part)
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
wD=double(Y.wD);
fD=wD/(2*pi); % in
JbarD=double(Y.JbarD);
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);
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
lambda_e=10;%0*200*15/1000; % excitatory synapses learning rate
lambda_i=1; % inhibitory synapses learning rate
H_E=-1; 
H_I=1;
%%%%%%%%% Generating the synapses matrix %%%%%%%%%
N_e=40; % Size of excitatory population
N_i=10; % Size of inhibitory population
%%% zero order synapses %%% 
Jiemean=6;%J_ie_final;
Jeimean=0.52;%JbarD^2/J_ie_final;
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
%%%%%%%% Building the vector around the zero order %%%%%%%%%
meb=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
mib=(1-Jeemean+Jiemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
gabs=(JbarD^2-Jeemean*Jiimean)^0.5;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ordmag=10^-2;%10^-11;

J_ei_std=Jeimean/10;
va_ei=[(q-a) (q-a) ;(-q+a) (-q+a); 0 0];
vq_ei=[0 0 ; 0 0 ; 0 0];
vd_ei=[-b b ; -b b ; -b b];
%check_ei=ordmag*vd_ei;%0.0001*[-1 1 ; -1 1 ; -1 1];%zeros(N_e-1,N_i)];%[ones(1,N_i) ; -ones(N_e-1,N_i)/(N_e-1)];%0.01*[1 -1 ; zeros(N_e-1,N_i)];%+normrnd(0,0.001,N_e,N_i);%% %;% ;%0.01*[[1 -1 zeros(1,N_i-2)] ; zeros(N_e-1,N_i)];
%Jei=Jeimean+check_ei;%

J_ei_noise=normrnd(0,J_ei_std,N_e,N_i);
Jei=(Jeimean*ones(N_e,N_i))-mean(J_ei_noise(:))+J_ei_noise;


J_ie_std=Jiemean/10;
va_ie=[-c c 0; -c c 0];
vq_ie=[-1 1 0 ; 0 0 0];
vd_ie=[-(d+q) -(d+q) -(d+q) ; d+q d+q d+q];
%check_ie=ordmag*vd_ie;%0.001*[-1 -1 -1 ; 1 1 1];%[ones(1,N_e) ; -ones(N_i-1,N_e)/(N_i-1)];% 0.0000001*[[1 -0.5 -0.5] ; [1 -0.5 -0.5]];%
%Jie=Jiemean+check_ie;
J_ie_noise=normrnd(0,J_ie_std,N_i,N_e);
Jie=Jiemean*ones(N_i,N_e)-mean(J_ie_noise(:))+J_ie_noise;

Jee=Jeemean*ones(N_e,N_e);
%Jee=Jee-diag(diag(Jee));

Jii=Jiimean*ones(N_i,N_i);
%Jii=Jii-diag(diag(Jii));
%%% Analytical solution %%%
J=[[Jee/N_e -Jei/N_i];[Jie/N_e -Jii/N_i]]; % connectivity matrix
diagonal=eye(N_e+N_i);
mfa=(diagonal-J)\ones(N_e+N_i,1); % The exact fixed point solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_e_history=0.2;
m_i_history=1.3;
[Corr_ie,Corr_ei,~,~,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended]=proj.common.Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf);
K_pE=1/tau_pE*exp(-Delta_extended*H_E/tau_pE).*heaviside(H_E*Delta_extended);
K_mE=1/tau_mE*exp(Delta_extended*H_E/tau_mE).*heaviside(-H_E*Delta_extended);
K_pI=1/tau_pI*exp(-Delta_extended*H_I/tau_pI).*heaviside(H_I*Delta_extended);
K_mI=1/tau_mI*exp(Delta_extended*H_I/tau_mI).*heaviside(-H_I*Delta_extended);
tic
t=1:14000;
%t_chunk=100;
numtrails=100;
%check_arr=[reshape(check_ei.',1,[]) reshape(check_ie.',1,[])];
check_arr=[reshape((Jei-mean(Jei(:))).',1,[]) reshape((Jie-mean(Jie(:))).',1,[])];
J_ei_dot=zeros(N_e,N_i);
J_ie_dot=zeros(N_i,N_e);
dtlearn=1;

mulvn=N_e*(N_i-1);
mulve=N_i*(N_e-1);
mulva=N_e-1;
mulvd=N_i-1;

% create the vectors
vn=[ones(N_i-1,1) -diag(ones(1,N_i-1))];
zvn=mat2cell(repmat(vn,1,N_e),N_i-1,N_i*ones(1,N_e));
vn=[blkdiag(zvn{:}) zeros(size((blkdiag(zvn{:})),1),2*N_e*N_i-size((blkdiag(zvn{:})),2))]'/2^0.5;

ve=[ones(N_e-1,1) -diag(ones(1,N_e-1))];
zve=mat2cell(repmat(ve,1,N_i),N_e-1,N_e*ones(1,N_i));
ve=[zeros(size((blkdiag(zve{:})),1),2*N_e*N_i-size((blkdiag(zve{:})),2)) blkdiag(zve{:})]'/2^0.5;

vaie=repmat([ones(N_e-1,1) -diag(ones(1,N_e-1))],1,N_i);
zvaei=mat2cell(repmat(ones(1,N_i),1,N_e-1),1,N_i*ones(1,N_e-1));
vaei=[ones(N_e-1,N_i) -blkdiag(zvaei{:})];
va=[(a-q)*vaei c*vaie]'/(((a-q)^2+c^2)*2*N_i)^0.5;

vdei=repmat([ones(N_i-1,1) -diag(ones(1,N_i-1))],1,N_e);
zvdie=mat2cell(repmat(ones(1,N_e),1,N_i-1),1,N_e*ones(1,N_i-1));
vdie=[ones(N_i-1,N_e) -blkdiag(zvdie{:})];
vd=[b*vdei (d+q)*vdie]'/((b^2+(d+q)^2)*2*N_e)^0.5;

vu=[ones(1,2*N_e*N_i) ; ones(1,N_e*N_i) -ones(1,N_e*N_i)]'/(2*N_e*N_i)^0.5;

parfor trail=1:numtrails
    [Vnvol,Vevol,Vavol,Vdvol,stdmebar,stdmibar,stdmetil,stdmitil,Jeimeandyn,Jiemeandyn,dJ_arr_t] = FullSynDynMultiTrails(m_e_history,m_i_history,Jei,Jie,Jee,Jii,t,Jiemax,mu,alpha,JbarD,Jeemean,Jiimean,vn,ve,vu,vaei,vaie,vdei,vdie,mulvn,mulve,mulva,mulvd,dtlearn,lambda_e,lambda_i,K_Ibar,K_Itilmphi,K_Itilphi,K_Eptilpsi,K_Etilphi,K_Etilmphi,K_pI,K_pE,K_mI,K_mE,N_e,N_i,dt,tf)
    parsave(sprintf('Figure6_simulationsdataNEW%d.mat',trail),Vnvol,Vevol,Vavol,Vdvol,stdmebar,stdmibar,stdmetil,stdmitil,Jeimeandyn,Jiemeandyn,dJ_arr_t);
end
toc
%%%% BELOW is the name of the previous working simulation!!!!!! %%%%
%%%%!!!! Figure5_simulationsdataNEW !!!!%%%%

%% The graphs of the previous sections for 100 trails
base_filename='Figure6_simulationsdataNEW';
numtrails=100;
Vnvol_alltrails=[];
Vevol_alltrails=[];
Vavol_alltrails=[];
Vdvol_alltrails=[];
stdmebar_alltrails=[];
stdmibar_alltrails=[];
stdmetil_alltrails=[];
stdmitil_alltrails=[];
Jeimeandynalltrails=[];
Jiemeandynalltrails=[];
dJarralltrails=[];
for i=1:numtrails
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
%% Computing the mean size squared (var ?) of each population
oldparam = sympref('HeavisideAtOrigin',1/2); % heaviside is 0.5 on zero
dt=0.01; % the accuracy through all the script and auxiliary functions (except for the nullclines part)
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
wD=double(Y.wD);
fD=wD/(2*pi); % in
JbarD=double(Y.JbarD);
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);

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
lambda_e=10;%0*200*15/1000; % excitatory synapses learning rate
lambda_i=1; % inhibitory synapses learning rate
H_E=-1; 
H_I=1;
%%%%%%%%% Generating the synapses matrix %%%%%%%%%
N_e=40; % Size of excitatory population
N_i=10; % Size of inhibitory population
%%% zero order synapses %%% 
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

J_ie_final=Jiemax*(1-((1-K_Ibar)*(K_Ibar*K_Emtilpsi-K_Itilpsi)/(K_Ibar*K_Eptilpsi-K_Itilpsi))^(1/mu)); % the place the


numtrails=100;
tflearn=14000;

Volsq1=zeros(numtrails,tflearn);
Volsq2=zeros(numtrails,tflearn);
Volsq3=zeros(numtrails,tflearn);
Volsq4=zeros(numtrails,tflearn);
Volsq5=zeros(numtrails,tflearn);


for i=1:numtrails
    parfor j=1:tflearn
        
        Jiemean=Jiemeandynalltrails(i,j);
        Jeimean=Jeimeandynalltrails(i,j);
        
        f=(1-Jiemean/Jiemax)^mu;
        ftag=-mu/Jiemax*(1-Jiemean/Jiemax)^(mu-1);
        
        K_Ebar=f-alpha;
        K_Etilpsi=f*K_Eptilpsi-alpha*K_Emtilpsi;
        K_Etilmphi=f*K_Eptilmphi-alpha*K_Emtilmphi;
        K_Etilphi=f*K_Eptilphi-alpha*K_Emtilphi;
        %%%%%%%% Building the vector around the zero order %%%%%%%%%
        meb=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
        mib=(1-Jeemean+Jiemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
        
        if Jeimean*Jiemean<JbarD^2
            a=-lambda_i*K_Ibar*mib^2;
            b=lambda_i*K_Ibar*meb^2;
            c=-lambda_e*K_Ebar*mib^2;
            d=lambda_e*K_Ebar*meb^2;
            q=lambda_e*ftag*meb*mib;
        else
            gabs=(JbarD^2-Jeemean*Jiimean)^0.5;
            mit=(JbarD^2/Jeimean*(1+Jiimean)-JbarD^2)/(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2));
            met=JbarD/Jiemean*mit;
            a=-lambda_i*(K_Ibar*mib^2+K_Itilmphi*mit^2/(2*gabs));
            b=lambda_i*(K_Ibar*meb^2+K_Itilphi*met^2/(2*gabs));
            c=-lambda_e*(K_Ebar*mib^2+K_Etilphi*mit^2/(2*gabs));
            d=lambda_e*(K_Ebar*meb^2+K_Etilmphi*met^2/(2*gabs));
            q=lambda_e*ftag*(meb*mib+met*mit/2*K_Eptilpsi); % self deperssion coefficient    end
        end
        % create the vectors
        vn=[ones(N_i-1,1) -diag(ones(1,N_i-1))];
        zvn=mat2cell(repmat(vn,1,N_e),N_i-1,N_i*ones(1,N_e));
        vn=[blkdiag(zvn{:}) zeros(size((blkdiag(zvn{:})),1),2*N_e*N_i-size((blkdiag(zvn{:})),2))]'/2^0.5;
        
        ve=[ones(N_e-1,1) -diag(ones(1,N_e-1))];
        zve=mat2cell(repmat(ve,1,N_i),N_e-1,N_e*ones(1,N_i));
        ve=[zeros(size((blkdiag(zve{:})),1),2*N_e*N_i-size((blkdiag(zve{:})),2)) blkdiag(zve{:})]'/2^0.5;
        
        vaie=repmat([ones(N_e-1,1) -diag(ones(1,N_e-1))],1,N_i);
        zvaei=mat2cell(repmat(ones(1,N_i),1,N_e-1),1,N_i*ones(1,N_e-1));
        vaei=[ones(N_e-1,N_i) -blkdiag(zvaei{:})];
        va=[(a-q)*vaei c*vaie]'/(((a-q)^2+c^2)*2*N_i)^0.5;
        
        vdei=repmat([ones(N_i-1,1) -diag(ones(1,N_i-1))],1,N_e);
        zvdie=mat2cell(repmat(ones(1,N_e),1,N_i-1),1,N_e*ones(1,N_i-1));
        vdie=[ones(N_i-1,N_e) -blkdiag(zvdie{:})];
        vd=[b*vdei (d+q)*vdie]'/((b^2+(d+q)^2)*2*N_e)^0.5;
        
        vu=[ones(1,2*N_e*N_i) ; ones(1,N_e*N_i) -ones(1,N_e*N_i)]'/(2*N_e*N_i)^0.5;
        
        
        vn=gsog(vn);
        va=gsog(va);
        vd=gsog(vd);
        vu=gsog(vu);
        ve=gsog(ve);
        
        vn=vn./repmat(sum(vn.^2,1).^0.5,N_e*N_i*2,1);
        va=va./repmat(sum(va.^2,1).^0.5,N_e*N_i*2,1);
        vd=vd./repmat(sum(vd.^2,1).^0.5,N_e*N_i*2,1);
        vu=vu./repmat(sum(vu.^2,1).^0.5,N_e*N_i*2,1);
        ve=ve./repmat(sum(ve.^2,1).^0.5,N_e*N_i*2,1);
        
        
        V=[vn  ve  va  vd  vu];
        
        vn_size=size(vn,2);
        ve_size=size(ve,2);
        va_size=size(va,2);
        vd_size=size(vd,2);
        vu_size=size(vu,2);
        
        
        dJ_arr=squeeze(dJarralltrails(i,:,j));
        
        
        coef=V\dJ_arr';
        coef1=diag([ones(1,vn_size) zeros(1,ve_size) zeros(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef;
        coef2=diag([zeros(1,vn_size) ones(1,ve_size) zeros(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef;
        coef3=diag([zeros(1,vn_size) zeros(1,ve_size) ones(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef;
        coef4=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) ones(1,vd_size) zeros(1,vu_size)])*coef;
        coef5=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) zeros(1,vd_size) ones(1,vu_size)])*coef;
        Volsq1(i,j)=sum(coef1.^2);
        Volsq2(i,j)=sum(coef2.^2);
        Volsq3(i,j)=sum(coef3.^2);
        Volsq4(i,j)=sum(coef4.^2);
        Volsq5(i,j)=sum(coef5.^2);
        
    end
end

save('data/VolsqFULL','Volsq1','Volsq2','Volsq3','Volsq4','Volsq5');
%% Plotting the Var of every family and the SD in each population
%load('VolsqFULL.mat', 'Volsq1', 'Volsq2', 'Volsq3', 'Volsq4', 'Volsq5')

Var1norm=mean(Volsq1/(N_e*(N_i-1)),1);
Var2norm=mean(Volsq2/(N_i*(N_e-1)),1);
Var3norm=mean(Volsq3/(N_e-1)^2,1);
Var4norm=mean(Volsq4/(N_i-1)^2,1);

% Plots of the previous section
figure(1)
stdshade(stdmebar_alltrails,0.3,[0.90, 0.40, 0.35],[],[])
hold on
stdshade(stdmibar_alltrails,0.3,[0.30, 0.60, 0.90],[],[])

%ylim([0 0.005])
ylabel('SD of $\bar{m}_{X,i}$ [a.u.]','interpreter','latex','FontSize',18)
xlabel('$t \ [\mathrm{a.u.}]$','interpreter','latex','FontSize',18)
lgd=legend({'','$\mathrm{E}$','','$\mathrm{I}$'},'Interpreter','latex','Location','Northeast');
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
xlim([0 14000])
grid on

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



%   figure;
%   plot(JbarD^2./(0:0.1:20),0:0.1:20,'Color','Black','LineWidth',3)
%   hold on
%   plot((1+Jiimean)*ones(1,length(0:20)),0:20,'Color','Black','LineWidth',3)
%   plot(JbarD^2/J_ie_final,J_ie_final,'*','Color',[0 0 1])
%   plot(Jeimeandynalltrails(1,:),Jiemeandynalltrails(1,:),'r.')
%   xlim([0 1.5])
%   grid on

 %figure;
% % plot(JbarD^2./(0:0.1:20),0:0.1:20,'Color','Black','LineWidth',3)
% % hold on
% % plot((1+Jiimean)*ones(1,length(0:20)),0:20,'Color','Black','LineWidth',3)
% % plot(JbarD^2/J_ie_final,J_ie_final,'*','Color',[0 0 1])
 %plot(Jeimeandynalltrails(1,:))
% grid on

%%% Analytical Results %%
oldparam = sympref('HeavisideAtOrigin',1/2); % heaviside is 0.5 on zero
dt=0.01; % the accuracy through all the script and auxiliary functions (except for the nullclines part)
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
wD=double(Y.wD);
fD=wD/(2*pi); % in
JbarD=double(Y.JbarD);
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);

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
lambda_e=10;%0*200*15/1000; % excitatory synapses learning rate
lambda_i=1; % inhibitory synapses learning rate
H_E=-1;
H_I=1;
%%%%%%%%% Generating the synapses matrix %%%%%%%%%
N_e=40; % Size of excitatory population
N_i=10; % Size of inhibitory population
%%% zero order synapses %%%
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

J_ie_final=Jiemax*(1-((1-K_Ibar)*(K_Ibar*K_Emtilpsi-K_Itilpsi)/(K_Ibar*K_Eptilpsi-K_Itilpsi))^(1/mu)); % the place the
Jiemean=J_ie_final;
Jeimean=JbarD^2/J_ie_final;

f=(1-Jiemean/Jiemax)^mu;
ftag=-mu/Jiemax*(1-Jiemean/Jiemax)^(mu-1);

K_Ebar=f-alpha;
K_Etilpsi=f*K_Eptilpsi-alpha*K_Emtilpsi;
K_Etilmphi=f*K_Eptilmphi-alpha*K_Emtilmphi;
K_Etilphi=f*K_Eptilphi-alpha*K_Emtilphi;
%%%%%%%% Building the vector around the zero order %%%%%%%%%
meb=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
mib=(1-Jeemean+Jiemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);


aFP=-lambda_i*K_Ibar*mib^2;
bFP=lambda_i*K_Ibar*meb^2;
cFP=-lambda_e*K_Ebar*mib^2;
dFP=lambda_e*K_Ebar*meb^2;
qFP=lambda_e*ftag*meb*mib;

gabs=(JbarD^2-Jeemean*Jiimean)^0.5;
mit=(JbarD^2/Jeimean*(1+Jiimean)-JbarD^2)/(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2));
met=JbarD/Jiemean*mit;
aR=-lambda_i*(K_Ibar*mib^2+K_Itilmphi*mit^2/(2*gabs));
bR=lambda_i*(K_Ibar*meb^2+K_Itilphi*met^2/(2*gabs));
cR=-lambda_e*(K_Ebar*mib^2+K_Etilphi*mit^2/(2*gabs));
dR=lambda_e*(K_Ebar*meb^2+K_Etilmphi*met^2/(2*gabs));
qR=lambda_e*ftag*(meb*mib+met*mit/2*K_Eptilpsi); % self deperssion coefficient    end
% create the vectors
for p=1:2
    if p==1
        a=aFP;
        b=bFP;
        c=cFP;
        d=dFP;
        q=qFP;
    else
        a=aR;
        b=bR;
        c=cR;
        d=dR;
        q=qR;
    end
    vn=[ones(N_i-1,1) -diag(ones(1,N_i-1))];
    zvn=mat2cell(repmat(vn,1,N_e),N_i-1,N_i*ones(1,N_e));
    vn=[blkdiag(zvn{:}) zeros(size((blkdiag(zvn{:})),1),2*N_e*N_i-size((blkdiag(zvn{:})),2))]'/2^0.5;
    
    ve=[ones(N_e-1,1) -diag(ones(1,N_e-1))];
    zve=mat2cell(repmat(ve,1,N_i),N_e-1,N_e*ones(1,N_i));
    ve=[zeros(size((blkdiag(zve{:})),1),2*N_e*N_i-size((blkdiag(zve{:})),2)) blkdiag(zve{:})]'/2^0.5;
    
    vaie=repmat([ones(N_e-1,1) -diag(ones(1,N_e-1))],1,N_i);
    zvaei=mat2cell(repmat(ones(1,N_i),1,N_e-1),1,N_i*ones(1,N_e-1));
    vaei=[ones(N_e-1,N_i) -blkdiag(zvaei{:})];
    va=[(a-q)*vaei c*vaie]'/(((a-q)^2+c^2)*2*N_i)^0.5;
    
    vdei=repmat([ones(N_i-1,1) -diag(ones(1,N_i-1))],1,N_e);
    zvdie=mat2cell(repmat(ones(1,N_e),1,N_i-1),1,N_e*ones(1,N_i-1));
    vdie=[ones(N_i-1,N_e) -blkdiag(zvdie{:})];
    vd=[b*vdei (d+q)*vdie]'/((b^2+(d+q)^2)*2*N_e)^0.5;
    
    vu=[ones(1,2*N_e*N_i) ; ones(1,N_e*N_i) -ones(1,N_e*N_i)]'/(2*N_e*N_i)^0.5;
    
    
    vn=gsog(vn);
    va=gsog(va);
    vd=gsog(vd);
    vu=gsog(vu);
    ve=gsog(ve);
    
    vn=vn./repmat(sum(vn.^2,1).^0.5,N_e*N_i*2,1);
    va=va./repmat(sum(va.^2,1).^0.5,N_e*N_i*2,1);
    vd=vd./repmat(sum(vd.^2,1).^0.5,N_e*N_i*2,1);
    vu=vu./repmat(sum(vu.^2,1).^0.5,N_e*N_i*2,1);
    ve=ve./repmat(sum(ve.^2,1).^0.5,N_e*N_i*2,1);
    
    if p==1
        VFP=[vn  ve  va  vd  vu];
    else
        VR=[vn  ve  va  vd  vu];
    end
end

vn_size=size(vn,2);
ve_size=size(ve,2);
va_size=size(va,2);
vd_size=size(vd,2);
vu_size=size(vu,2);

MEEbar=lambda_e*K_Ebar*meb^2;
MIEbar=lambda_i*K_Ibar*meb^2;
MEEtil=MEEbar+lambda_e*K_Etilmphi*met^2/(2*gabs);
MIEtil=MIEbar+lambda_i*K_Itilmphi*met^2/(2*gabs);
FEbar=lambda_e*ftag*meb*mib;
FEtil=FEbar+lambda_e*ftag*K_Eptilpsi*met*mit/2;

sigFP4=MEEbar+FEbar;
sigR4=MEEtil+FEtil;
etaFP4=(MEEbar+FEbar)/MIEbar;
etaR4=(MEEtil+FEtil)/MIEtil;
Ta=etaFP4/etaR4*((1+etaR4^2)/(1+etaFP4^2))^0.5;

noiseEIVar=10^-4;
noiseIEVar=10^-2;

xitrails=1000;
%coef4FPoneEig=0; % the component of one eigenv
%coef4RoneEig=0;
coef4FPall=0;
coef4Rall=0;

% for k=1:xitrails
% xiEI=normrnd(0,noiseEIVar^0.5,N_e*N_i,1);
% xiIE=normrnd(0,noiseIEVar^0.5,N_e*N_i,1);
% xi=[xiEI ; xiIE];
% 
% coefFP=VFP\xi;
% coef4FP=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) ones(1,vd_size) zeros(1,vu_size)])*coefFP;
% %coef4FPoneEig=coef4FPoneEig+coef4FP(end-2)^2;
% coef4FPall=coef4FPall+sum(coef4FP.^2,"all");
% 
% xiEI=normrnd(0,noiseEIVar^0.5,N_e*N_i,1);
% xiIE=normrnd(0,noiseIEVar^0.5,N_e*N_i,1);
% xi=[xiEI ; xiIE];
% 
% coefR=VR\xi;
% coef4R=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) ones(1,vd_size) zeros(1,vu_size)])*coefR;
% %coef4RoneEig=coef4RoneEig+coef4R(end-2)^2;
% coef4Rall=coef4Rall+sum(coef4R.^2,"all");
% end

% xibarsqmean=coef4FPall/(xitrails*(N_i-1));
% xitilsqmean=coef4Rall/(xitrails*(N_i-1));

VFPproj=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) ones(1,vd_size) zeros(1,vu_size)])/VFP;
Cxibarsq=diag([noiseEIVar*ones(1,N_e*N_i) noiseIEVar*ones(1,N_e*N_i)]);
xibarsqmean=trace(transpose(VFPproj)*VFPproj*Cxibarsq)/(N_i-1);

VRproj=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) ones(1,vd_size) zeros(1,vu_size)])/VR;
Cxitilsq=diag([noiseEIVar*ones(1,N_e*N_i) noiseIEVar*ones(1,N_e*N_i)]);
xitilsqmean=trace(transpose(VRproj)*VRproj*Cxitilsq)/(N_i-1);


JeidotFP=lambda_i*K_Ibar.*meb.*mib;
JeidotR=JeidotFP+lambda_i*K_Itilpsi.*met.*mit/2;

c4varFP=(xibarsqmean/abs(JeidotFP)+xitilsqmean/(abs(JeidotR)*Ta^2))/(-2*(sigFP4/JeidotFP+sigR4/abs(JeidotR)));
c4varR=(xibarsqmean*Ta^2/abs(JeidotFP)+xitilsqmean/(abs(JeidotR)))/(-2*(sigFP4/JeidotFP+sigR4/abs(JeidotR)));

c4var=(c4varFP/abs(JeidotFP)+c4varR/abs(JeidotR))/(1/abs(JeidotFP)+1/abs(JeidotR));


figure(3)
plot(1:14000,c4var*ones(1,length(1:14000)),'--','Color','black','LineWidth',2)
hold on

lgd=legend({'','$\mathrm{I}$','','$\mathrm{II}$','','$\mathrm{III}$','','$\mathrm{IV}$','',''},'Interpreter','latex','Location','Northeast');


%% Computing the autocorrelations 
% We assume here that we arrived Cri. Rhythm. and compute all the coefficients a,b,c,d,q in FP and R only once 
oldparam = sympref('HeavisideAtOrigin',1/2); % heaviside is 0.5 on zero
dt=0.01; % the accuracy through all the script and auxiliary functions (except for the nullclines part)
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
wD=double(Y.wD);
fD=wD/(2*pi); % in
JbarD=double(Y.JbarD);
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);

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
lambda_e=10;%0*200*15/1000; % excitatory synapses learning rate
lambda_i=1; % inhibitory synapses learning rate
H_E=-1; 
H_I=1;
%%%%%%%%% Generating the synapses matrix %%%%%%%%%
N_e=40; % Size of excitatory population
N_i=10; % Size of inhibitory population
%%% zero order synapses %%% 
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

J_ie_final=Jiemax*(1-((1-K_Ibar)*(K_Ibar*K_Emtilpsi-K_Itilpsi)/(K_Ibar*K_Eptilpsi-K_Itilpsi))^(1/mu)); % the place the
Jiemean=J_ie_final;
Jeimean=JbarD^2/Jiemean;
        
f=(1-Jiemean/Jiemax).^mu;
ftag=-mu/Jiemax*(1-Jiemean/Jiemax)^(mu-1);
        
K_Ebar=f-alpha;
K_Etilpsi=f*K_Eptilpsi-alpha*K_Emtilpsi;
K_Etilmphi=f*K_Eptilmphi-alpha*K_Emtilmphi;
K_Etilphi=f*K_Eptilphi-alpha*K_Emtilphi;
%%%%%%%% Building the vector around the zero order %%%%%%%%%
meb=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
mib=(1-Jeemean+Jiemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);

aFP=-lambda_i*K_Ibar*mib^2;
bFP=lambda_i*K_Ibar*meb^2;
cFP=-lambda_e*K_Ebar*mib^2;
dFP=lambda_e*K_Ebar*meb^2;
qFP=lambda_e*ftag*meb*mib;

gabs=(JbarD^2-Jeemean*Jiimean)^0.5;
mit=(JbarD^2/Jeimean*(1+Jiimean)-JbarD^2)/(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2));
met=JbarD/Jiemean*mit;
aR=-lambda_i*(K_Ibar*mib^2+K_Itilmphi*mit^2/(2*gabs));
bR=lambda_i*(K_Ibar*meb^2+K_Itilphi*met^2/(2*gabs));
cR=-lambda_e*(K_Ebar*mib^2+K_Etilphi*mit^2/(2*gabs));
dR=lambda_e*(K_Ebar*meb^2+K_Etilmphi*met^2/(2*gabs));
qR=lambda_e*ftag*(meb*mib+met*mit/2*K_Eptilpsi); % self deperssion coefficient    end


for p=1:2
    
    if p==1
       a=aFP;
       b=bFP;
       c=cFP;
       d=dFP;
       q=qFP;
    else
       a=aR;
       b=bR;
       c=cR;
       d=dR;
       q=qR;        
    end
    
    % create the vectors
    vn=[ones(N_i-1,1) -diag(ones(1,N_i-1))];
    zvn=mat2cell(repmat(vn,1,N_e),N_i-1,N_i*ones(1,N_e));
    vn=[blkdiag(zvn{:}) zeros(size((blkdiag(zvn{:})),1),2*N_e*N_i-size((blkdiag(zvn{:})),2))]'/2^0.5;
    
    ve=[ones(N_e-1,1) -diag(ones(1,N_e-1))];
    zve=mat2cell(repmat(ve,1,N_i),N_e-1,N_e*ones(1,N_i));
    ve=[zeros(size((blkdiag(zve{:})),1),2*N_e*N_i-size((blkdiag(zve{:})),2)) blkdiag(zve{:})]'/2^0.5;
    
    vaie=repmat([ones(N_e-1,1) -diag(ones(1,N_e-1))],1,N_i);
    zvaei=mat2cell(repmat(ones(1,N_i),1,N_e-1),1,N_i*ones(1,N_e-1));
    vaei=[ones(N_e-1,N_i) -blkdiag(zvaei{:})];
    va=[(a-q)*vaei c*vaie]'/(((a-q)^2+c^2)*2*N_i)^0.5;
    
    vdei=repmat([ones(N_i-1,1) -diag(ones(1,N_i-1))],1,N_e);
    zvdie=mat2cell(repmat(ones(1,N_e),1,N_i-1),1,N_e*ones(1,N_i-1));
    vdie=[ones(N_i-1,N_e) -blkdiag(zvdie{:})];
    vd=[b*vdei (d+q)*vdie]'/((b^2+(d+q)^2)*2*N_e)^0.5;
    
    vu=[ones(1,2*N_e*N_i) ; ones(1,N_e*N_i) -ones(1,N_e*N_i)]'/(2*N_e*N_i)^0.5;
    
    
    vn=gsog(vn);
    va=gsog(va);
    vd=gsog(vd);
    vu=gsog(vu);
    ve=gsog(ve);
    
    vn=vn./repmat(sum(vn.^2,1).^0.5,N_e*N_i*2,1);
    va=va./repmat(sum(va.^2,1).^0.5,N_e*N_i*2,1);
    vd=vd./repmat(sum(vd.^2,1).^0.5,N_e*N_i*2,1);
    vu=vu./repmat(sum(vu.^2,1).^0.5,N_e*N_i*2,1);
    ve=ve./repmat(sum(ve.^2,1).^0.5,N_e*N_i*2,1);
    
    if p==1
        VFP=[vn  ve  va  vd  vu];
    else
        VR=[vn  ve  va  vd  vu];
    end
    
end

VR=sign(VR.*VFP).*VR; % We want the eigenvectors in each region to have the same signs

vn_size=size(vn,2);
ve_size=size(ve,2);
va_size=size(va,2);
vd_size=size(vd,2);
vu_size=size(vu,2);

numtrails=100;
tflearn=14000;

t0acorr=2000;
tfacorr=14000;
maxlag=900;
acorr4=zeros(1,length(0:maxlag));


for i=1:numtrails
        

        dJ_arr=squeeze(dJarralltrails(i,:,t0acorr:tfacorr))';
        yFP=heaviside(JbarD^2-squeeze(Jeimeandynalltrails(i,t0acorr:tfacorr)).*squeeze(Jiemeandynalltrails(i,t0acorr:tfacorr)))';
        repyFP=repmat(yFP,1,N_e*N_i*2);
        dJ_arrFP=dJ_arr.*repyFP;
        
        yR=heaviside(-JbarD^2+squeeze(Jeimeandynalltrails(i,t0acorr:tfacorr)).*squeeze(Jiemeandynalltrails(i,t0acorr:tfacorr)))';
        repyR=repmat(yR,1,N_e*N_i*2);
        dJ_arrR=dJ_arr.*repyR;
        
        coef=VFP\dJ_arrFP'+VR\dJ_arrR';
        coef1=diag([ones(1,vn_size) zeros(1,ve_size) zeros(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef;
        coef2=diag([zeros(1,vn_size) ones(1,ve_size) zeros(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef;
        coef3=diag([zeros(1,vn_size) zeros(1,ve_size) ones(1,va_size) zeros(1,vd_size) zeros(1,vu_size)])*coef;
        coef4=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) ones(1,vd_size) zeros(1,vu_size)])*coef;
        coef5=diag([zeros(1,vn_size) zeros(1,ve_size) zeros(1,va_size) zeros(1,vd_size) ones(1,vu_size)])*coef;
        
        coef4st=coef4(:,1:end-1000);
        Tint=size(coef4st,2);
        
        parfor j=0:maxlag
        acorr4parfor(1+j)=trace(coef4st*coef4(:,(1+j):(j+Tint))')/(Tint*vd_size);
        end
        
        acorr4(i,:)=acorr4parfor;
end

save('data/AutoCorrV4','acorr4');

%%
figure(4)

stdshade(acorr4./acorr4(:,1),0.3,[0.90, 0.40, 0.35],[],[])
hold on
acorr0=c4var;
tau=0:900;
sig4=(sigFP4/JeidotFP+sigR4/abs(JeidotR))/(1/JeidotFP+1/abs(JeidotR));
typtime=1/abs(sig4);
plot(tau,exp(tau*(sigFP4/JeidotFP+sigR4/abs(JeidotR))/(1/JeidotFP+1/abs(JeidotR))),'--','Color','black')
xlabel('$\tau \ [\mathrm{a.u.}]$','interpreter','latex','FontSize',18)
ylabel('$\mathrm{\mathcal{\rho}_{IV}(\tau)}$','interpreter','latex','FontSize',18)
%lgd=legend({'$V_{0}$','$V_{E}$','$V_{A}$','$V_{D+E}$'},'Interpreter','latex','Location','Northeast');
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
grid on
xlim([0 750])

grid on