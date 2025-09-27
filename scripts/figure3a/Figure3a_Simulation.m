%% Nullclines in 2D  for order parameters - parallel computing %%
D=0.4;
dt=0.01;
tf=200;
Jeemean=0.6;
Jiimean=0.4;
Jiemean=6;%J_ie_final;
Jeimean=0.52;%JbarD^2/J_ie_final;
Jhat=(Jeemean+Jiimean)/2;
syms wD JbarD
if Jeemean>Jiimean
range=[0.1 5 ;0.01 pi/(2*D)];
Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
elseif Jeemean<Jiimean
range=[0.1 5 ;0.01 pi/D];    
Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
%else 
%range=[0.1 5 ;0.01 pi/D];  
%Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), T*wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
end    
wD=double(Y.wD);
fD=wD/(2*pi); % in
JbarD=double(Y.JbarD);
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);
%%%%%%%%% STDP Parameters %%%%%%%%%
alpha=0.96; % relative depression
mu=0.01;%0.015;%0.15;%0.015; % measure of linearity
Jiemax=20; % J_ie_max
tau_pE=2; % typical potentiation time of excitatory synapses
tau_pI=2; % typical potentiation time of inhibitory synapses
tau_mE=7; % typical depression time of excitatory synapses
tau_mI=3; % typical depression time of inhibitory synapses
thetapI=acos((1+(wD*tau_pI)^2)^-0.5);
thetamI=acos((1+(wD*tau_mI)^2)^-0.5);
thetapE=acos((1+(wD*tau_pE)^2)^-0.5);
thetamE=acos((1+(wD*tau_mE)^2)^-0.5);
lambda_e=1;%0*200*15/1000; % excitatory synapses learning rate
lambda_i=0.1; % inhibitory synapses learning rate
H_E=-1; 
H_I=1;

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

[Corr_ie_final,Corr_ei_final,~,~,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended]=proj.common.Correlations_2D_full_diff(Jeemean,Jeimean,Jiemean,Jiimean,dt,tf);
K_pE=1/tau_pE*exp(-Delta_extended*H_E/tau_pE).*heaviside(H_E*Delta_extended);
K_mE=1/tau_mE*exp(Delta_extended*H_E/tau_mE).*heaviside(-H_E*Delta_extended);
K_pI=1/tau_pI*exp(-Delta_extended*H_I/tau_pI).*heaviside(H_I*Delta_extended);
K_mI=1/tau_mI*exp(Delta_extended*H_I/tau_mI).*heaviside(-H_I*Delta_extended);
[Nullei,Nullie,dJei,dJie,~,~]=NullclinesParallelPt1(JbarD,lambda_i,lambda_e,alpha,tau_pE,tau_mE,tau_pI,tau_mI,H_E,H_I,K_pE,K_mE,K_pI,K_mI,mu,Jiemax,Jeemean,Jiimean,dt,tf);

%%% Vector flow %%%
dt=0.001;
tic
[Jeidot,Jiedot,Jei_arr,Jie_arr]=JdotPhDiagPt1(JbarD,lambda_i,lambda_e,alpha,tau_mE,tau_mI,tau_pE,tau_pI,H_E,H_I,mu,Jiemax,Jeemean,Jiimean,dt,tf);
toc

[Jei_grid,Jie_grid]=meshgrid(Jei_arr,Jie_arr);

Xei=Nullei(:,1);
Xei=Xei(~isnan(Xei));
[Xeisorted,Ind]=sort(Xei);
Yei=Nullei(:,2);
Yei=Yei(~isnan(Yei));
Yeisorted=Yei(Ind);
JeinullX=linspace(Xeisorted(1),Xeisorted(end),length(Xeisorted)*5)';
JeinullY=interp1(Xeisorted,Yeisorted,JeinullX,'spline');

Xie=Nullie(:,1);
Xie=Xie(~isnan(Xie));
[Xiesorted,Ind]=sort(Xie);
Yie=Nullie(:,2);
Yie=Yie(~isnan(Yie));
Yiesorted=Yie(Ind);
JienullX=linspace(Xiesorted(1),Xiesorted(end),length(Xiesorted)*5)';
JienullY=interp1(Xiesorted,Yiesorted,JienullX,'spline');
%%
figure;
plot(JeinullX,JeinullY,'red','LineWidth',2.5)
hold on
plot(JienullX,JienullY,'green','LineWidth',2.5)
%plot(Nullie(:,1),Nullie(:,2))
plot(0:0.01:1+Jiimean,JbarD^2./(0:0.01:1+Jiimean),'Color','Black','LineWidth',2);
%plot(JbarD^2/J_ie_final,J_ie_final,'*','Color',[0 0 1])
%plot(JbarD^2/J_ie_null_ei_leave,J_ie_null_ei_leave,'*','Color',[1 0 0])
%plot(JbarD^2/J_ie_null_ie_leave,J_ie_null_ie_leave,'*','Color',[0 1 0])
q=quiver(Jei_grid',Jie_grid',Jeidot,Jiedot,0.3);
%q=quiver(Jei_grid',Jie_grid',Jeidot,Jiedot,1);
q.Color = '[0 0.4470 0.7410]';
xlim([0 1+Jiimean]);%1+mean(J_ii(:))])%xlim([0 1+mean(J_ii(:))])
ylim([0 20])

%%% Saving for figure %%%
Jeidotall=1000*(Jeidot(:));
Jiedotall=1000*(Jiedot(:));
Jieall=round(Jie_grid',7);
Jeiall=round(Jei_grid',7);
Jeiall=Jeiall(:);
Jieall=Jieall(:);

save('vectorflow/Jeidot','Jeidotall','-double');
save('vectorflow/Jiedot','Jiedotall','-double');
save('vectorflow/Jeiarr','Jeiall','-double');
save('vectorflow/Jiearr','Jieall','-double');
save('vectorflow/JienullX','JienullX','-double');
save('vectorflow/JienullY','JienullY','-double');
save('vectorflow/JeinullX','JeinullX','-double');
save('vectorflow/JeinullY','JeinullY','-double');
