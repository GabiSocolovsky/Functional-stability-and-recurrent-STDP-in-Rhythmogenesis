%% What are the eigenvalues? run of tau_mE %% Figure 3c
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
gabs=(JbarD^2-Jeemean*Jiimean)^0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_e=1;
lambda_i=1;
alpha_i=0.93;
dalpha=0.0001;%0.001;
alpha_f=0.999;%0.999;
alpha_vec=alpha_i:dalpha:alpha_f; % relative depression %
mu=0.01;%0.015;%0.15;%0.015; % measure of linearity
Jiemax=20; % J_ie_max
tau_pE=2; % typical potentiation time of excitatory synapses
tau_pI=2; % typical potentiation time of inhibitory synapses

tau_mI_i=2;
dtau_mI=0.5;
tau_mI_f=5;
tau_mI_vec=tau_mI_i:dtau_mI:tau_mI_f;
q=1; % index for different tau_mI

alpha_cr=nan(length(tau_mI_vec),length(alpha_vec));
tau_mE_cr=nan(length(tau_mI_vec),length(alpha_vec));
for tau_mI=tau_mI_i:dtau_mI:tau_mI_f
tau_mE_i=tau_mI+0.001;
dtau_mE=0.005;
tau_mE_f=12;
tau_mE_vec=tau_mE_i:dtau_mE:tau_mE_f; % typical depression time of excitatory synapses

    
thetapI=acos((1+(wD*tau_pI)^2)^-0.5);
thetamI=acos((1+(wD*tau_mI)^2)^-0.5);
thetapE=acos((1+(wD*tau_pE)^2)^-0.5);
K_Iptilpsi=cos(thetapI)*cos(thetapI+psi);
K_Imtilpsi=cos(thetamI)*cos(thetamI-psi);
K_Eptilpsi=cos(thetapE)*cos(thetapE+psi);


[alpha_arr,tau_mE_arr]=meshgrid(alpha_vec,tau_mE_vec);



eq=Jiemax*(1-(alpha_arr.*((1-alpha_arr).*cos(acos((1+(wD*tau_mE_arr).^2).^-0.5)).*cos(acos((1+(wD*tau_mE_arr).^2).^-0.5)-psi)-(K_Iptilpsi-alpha_arr.*K_Imtilpsi))./((1-alpha_arr).*K_Eptilpsi-(K_Iptilpsi-alpha_arr.*K_Imtilpsi))).^(1/mu))-JbarD*((JbarD-2*gabs^2*(1-Jeemean)*(1-alpha_arr)./(K_Iptilpsi-alpha_arr*K_Imtilpsi))./(1+Jiimean+2*JbarD*gabs^2*(1-alpha_arr)./(K_Iptilpsi-alpha_arr.*K_Imtilpsi)));
[val,loc]=min(abs(eq),[],1);

for i=1:length(val)
    if val(i)>0.1
        continue
    end
    tau_mE_cr(q,i)=tau_mE_arr(loc(i));
    alpha_cr(q,i)=alpha_vec(i);
end


thetamE=acos((1+(wD*tau_mE_arr).^2).^-0.5);
K_Ibar=1-alpha_arr;
K_Iptilphi=cos(thetapI)*cos(thetapI+phi);
K_Imtilphi=cos(thetamI)*cos(thetamI-phi);
K_Itilphi=K_Iptilphi-alpha_arr.*K_Imtilphi;
K_Iptilmphi=cos(thetapI)*cos(thetapI-phi);
K_Imtilmphi=cos(thetamI)*cos(thetamI+phi);
K_Itilmphi=K_Iptilmphi-alpha_arr.*K_Imtilmphi;

K_Iptilpsi=cos(thetapI)*cos(thetapI+psi);
K_Imtilpsi=cos(thetamI)*cos(thetamI-psi);
K_Itilpsi=K_Iptilpsi-alpha_arr.*K_Imtilpsi;

K_Eptilphi=cos(thetapE)*cos(thetapE-phi);
K_Emtilphi=cos(thetamE).*cos(thetamE+phi);
K_Eptilmphi=cos(thetapE)*cos(thetapE+phi);
K_Emtilmphi=cos(thetamE).*cos(thetamE-phi);

K_Eptilpsi=cos(thetapE)*cos(thetapE+psi);
K_Emtilpsi=cos(thetamE).*cos(thetamE-psi);

J_ie_final=Jiemax*(1-((1-K_Ibar).*(K_Ibar.*K_Emtilpsi-K_Itilpsi)./(K_Ibar.*K_Eptilpsi-K_Itilpsi)).^(1/mu));
J_ie_final(J_ie_final<JbarD^2/(1+Jiimean))=nan;
J_ei_final=JbarD^2./J_ie_final;
f=(1-J_ie_final/Jiemax).^mu;
ftag=-mu/Jiemax*(1-J_ie_final/Jiemax).^(mu-1);

K_Ebar=f-alpha_arr;
K_Etilphi=f.*K_Eptilphi-alpha_arr.*K_Emtilphi;
K_Etilmphi=f.*K_Eptilmphi-alpha_arr.*K_Emtilmphi;
K_Etilpsi=f.*K_Eptilpsi-alpha_arr.*K_Emtilpsi;

gabs=(JbarD^2-Jeemean*Jiimean)^0.5;
meb=(1+Jiimean-J_ei_final)./((1-Jeemean)*(1+Jiimean)+J_ei_final.*J_ie_final);
mib=(1-Jeemean+J_ie_final)./((1-Jeemean)*(1+Jiimean)+J_ei_final.*J_ie_final);
mit=(J_ie_final*(1+Jiimean)-JbarD^2)./(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2));    
met=JbarD./J_ie_final.*mit;
Atil=-lambda_i*(K_Ibar.*mib.^2+K_Itilmphi.*mit.^2/(2*gabs));
FbarE=lambda_e*ftag.*meb.*mib;
FtilE=FbarE+lambda_e*ftag.*met.*mit/2.*K_Eptilpsi;


MbarEE=lambda_e*K_Ebar.*meb.^2;
Sigbar4=MbarEE+FbarE;

MtilEE=MbarEE+lambda_e*K_Etilmphi.*met.^2/(2*gabs);
Sigtil4=MtilEE+FtilE;



JeidotFP=lambda_i*K_Ibar.*meb.*mib;
JeidotR=JeidotFP+lambda_i*K_Itilpsi.*met.*mit/2;


q=q+1;

end
for q=1:length(tau_mI_vec)
plot(alpha_cr(q,:),tau_mE_cr(q,:)*5,'color',[0.8500 0.3250 0.0980]*(q/length(tau_mI_vec)),'LineWidth',3)
hold on
end
lgd=legend({'$10$','$12.5$','$15$','$17.5$','$20$','$22.5$','$25$'},'Interpreter','latex','Location','NorthWest');
lgd.Title.String = '$\tau_{I,-} \ [\mathrm{ms}]$';
xlabel('$\alpha$','interpreter','latex','FontSize',18)
ylabel('$\tau_{E,-} \ [\mathrm{ms]}$','interpreter','latex','FontSize',18)
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')
grid on