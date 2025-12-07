%% Figures 5a-5b Simulation %%
% This script computes the data for **Figures 5a-5b** from the paper:
% "Functional stability and recurrent STDP in rhythmogenesis"
%
% Description:

%       - Computes the eigenvalues of family IV when in CR as a function of alpha for
%       different values of tau_Em, for both the FP and R regions (figure 5a)
%       - Computes the effective eigenvalue of family IV as a function
%       alpha for different values of tauEm

% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Definitions %%%%%%%%%%%%%
oldparam = sympref("HeavisideAtOrigin",1/2); % 0.5 at the origin for heaviside
%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%
%%%%% Phase diagram features (bif. etc.) %%%%%
T=1; % time constant 5msec tau 
D=0.4; % delay in msec
Jiimean=0.4; % Jee order parameter
Jeemean=0.6; % Jii order parameter
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
fD=wD/(2*pi); % % frequency on bif. line
JbarD=double(Y.JbarD); % J bar on bif. line
phi=acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5));
psi=acos(Jhat/JbarD);
gabs=(JbarD^2-Jeemean*Jiimean)^0.5; % absolute value of g (see paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% STDP parameter %%%%%%%%%%%%%%%%%%
lambda_e=1;
lambda_i=1;
alpha_i=0.93; % initial value of alpha
dalpha=0.0001;
alpha_f=0.999; % final value of alpha
alpha_vec=alpha_i:dalpha:alpha_f; % relative depression %
mu=0.01; % measure of linearity
Jiemax=20; % J_ie_max
tau_pE=2; % typical potentiation time of excitatory synapses
tau_pI=2; % typical potentiation time of inhibitory synapses
tau_mI_i=3; % inital value tau_mI
dtau_mI=0.5; % final value of tau_mI
tau_mI_f=3;
tau_mI_vec=tau_mI_i:dtau_mI:tau_mI_f;
q=1; % index for different tau_mI

alpha_cr=nan(length(tau_mI_vec),length(alpha_vec)); % critical value of alpha
tau_mE_cr=nan(length(tau_mI_vec),length(alpha_vec)); % critical value of tau_mE
for tau_mI=tau_mI_i:dtau_mI:tau_mI_f % run on different tau_mI
tau_mE_i=tau_mI;
dtau_mE=0.005;
tau_mE_f=12;
tau_mE_vec=tau_mE_i:dtau_mE:tau_mE_f; % typical depression time of excitatory synapses

%%% analytical expression from paper %%%%        
thetapI=acos((1+(wD*tau_pI)^2)^-0.5);
thetamI=acos((1+(wD*tau_mI)^2)^-0.5);
thetapE=acos((1+(wD*tau_pE)^2)^-0.5);
K_Iptilpsi=cos(thetapI)*cos(thetapI+psi);
K_Imtilpsi=cos(thetamI)*cos(thetamI-psi);
K_Eptilpsi=cos(thetapE)*cos(thetapE+psi);


[alpha_arr,tau_mE_arr]=meshgrid(alpha_vec,tau_mE_vec);

%%% Solve equation for these values %%%
eq=Jiemax*(1-(alpha_arr.*((1-alpha_arr).*cos(acos((1+(wD*tau_mE_arr).^2).^-0.5)).*cos(acos((1+(wD*tau_mE_arr).^2).^-0.5)-psi)-(K_Iptilpsi-alpha_arr.*K_Imtilpsi))./((1-alpha_arr).*K_Eptilpsi-(K_Iptilpsi-alpha_arr.*K_Imtilpsi))).^(1/mu))-JbarD*((JbarD-2*gabs^2*(1-Jeemean)*(1-alpha_arr)./(K_Iptilpsi-alpha_arr*K_Imtilpsi))./(1+Jiimean+2*JbarD*gabs^2*(1-alpha_arr)./(K_Iptilpsi-alpha_arr.*K_Imtilpsi)));
[val,loc]=min(abs(eq),[],1);

for i=1:length(val)
    if val(i)>0.1
        continue
    end
    tau_mE_cr(q,i)=tau_mE_arr(loc(i)); % found tau_mE critical
    alpha_cr(q,i)=alpha_vec(i); % found alpha critical
end


% More analytical expressions from paper (m (p) is minus (plus)) %
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

J_ie_final=Jiemax*(1-((1-K_Ibar).*(K_Ibar.*K_Emtilpsi-K_Itilpsi)./(K_Ibar.*K_Eptilpsi-K_Itilpsi)).^(1/mu)); % Jie asteerisk
J_ie_final(J_ie_final<JbarD^2/(1+Jiimean))=nan;
J_ei_final=JbarD^2./J_ie_final; % Jei asterisk
f=(1-J_ie_final/Jiemax).^mu; % synaptic depndent function
ftag=-mu/Jiemax*(1-J_ie_final/Jiemax).^(mu-1); % f derivative in Jie

K_Ebar=f-alpha_arr;
K_Etilphi=f.*K_Eptilphi-alpha_arr.*K_Emtilphi;
K_Etilmphi=f.*K_Eptilmphi-alpha_arr.*K_Emtilmphi;
K_Etilpsi=f.*K_Eptilpsi-alpha_arr.*K_Emtilpsi;

gabs=(JbarD^2-Jeemean*Jiimean)^0.5; % absolute value of g
meb=(1+Jiimean-J_ei_final)./((1-Jeemean)*(1+Jiimean)+J_ei_final.*J_ie_final); % m_e bar
mib=(1-Jeemean+J_ie_final)./((1-Jeemean)*(1+Jiimean)+J_ei_final.*J_ie_final); % m_i bar
mit=(J_ie_final*(1+Jiimean)-JbarD^2)./(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2)); % m_i_tilde    
met=JbarD./J_ie_final.*mit; % m_e_tilde
FbarE=lambda_e*ftag.*meb.*mib; % F bar E
FtilE=FbarE+lambda_e*ftag.*met.*mit/2.*K_Eptilpsi; % F tilde E


MbarEE=lambda_e*K_Ebar.*meb.^2; % M bar EE
Sigbar4=(MbarEE+FbarE); % eigenvalue from family IV in FP region


MtilEE=MbarEE+lambda_e*K_Etilmphi.*met.^2/(2*gabs); % M tilde EE
Sigtil4=MtilEE+FtilE;  % eigenvalue from family IV in FP region

JeidotFP=lambda_i*K_Ibar.*meb.*mib; % Jei time derivative in FP region
JeidotR=JeidotFP+lambda_i*K_Itilpsi.*met.*mit/2; % Jei time derivative in R region

q=q+1;

end

%% Plot figure 5a %%
figure;

selected_taus = [4 5 6 7 8 9 10];
nCurves = length(selected_taus);
cmap = parula(nCurves);
legend_handles = gobjects(1, nCurves);  % to store handles for legend
curve_idx = 1;

for i = 1:length(tau_mE_vec)
    if ismember(tau_mE_vec(i), selected_taus)
        % Clean invalid values
        for j = 1:length(alpha_vec)
            if tau_mE_arr(i,j) < tau_mE_cr(j) || isnan(tau_mE_cr(j))
                Sigbar4(i,j) = nan;
                Sigtil4(i,j) = nan;
            end
        end

        % Plot Sigbar4 (solid), and store handle for legend
        h1 = plot(alpha_arr(i,:), Sigbar4(i,:), ...
            'LineWidth', 2, ...
            'Color', cmap(curve_idx,:));
        hold on

        % Plot Sigtil4 (dashed), without adding to legend
        plot(alpha_arr(i,:), Sigtil4(i,:), ...
            'LineWidth', 2, ...
            'LineStyle', '--', ...
            'Color', cmap(curve_idx,:));

        legend_handles(curve_idx) = h1;  % store only solid line handle
        curve_idx = curve_idx + 1;
    end
end

xlabel('$\alpha$', 'interpreter','latex','FontSize',18)
ylabel('$\bar{\sigma}_{\mathrm{IV}}$ (solid) \& $\tilde{\sigma}_{\mathrm{IV}}$ (dashed)', ...
       'interpreter','latex','FontSize',18)

legend_labels = arrayfun(@(x) sprintf('$%d$', x*5), selected_taus, 'UniformOutput', false);
lgd = legend(legend_handles, legend_labels, 'Interpreter','latex','Location','NorthEast');
lgd.Title.String = '$\tau_{E,-} \ [\mathrm{ms}]$';

set(gca,'FontSize',13)
set(gca,'TickLabelInterpreter','latex')
grid on

%% Plot Figure 5b: Normalized sum of Sigbar4 and Sigtil4 %%
figure;

selected_taus = [4 5 6 7 8 9 10];
nCurves = length(selected_taus);
cmap = parula(nCurves);
legend_handles = gobjects(1, nCurves);
curve_idx = 1;

% Clean invalid values
Sigbar4((Sigbar4<0) & (Sigtil4<0)) = nan;
Sigtil4((Sigbar4<0) & (Sigtil4<0)) = nan;

for i = 1:length(tau_mE_vec)
    if ismember(tau_mE_vec(i), selected_taus)
        % Compute normalized sum
        norm_sum = (Sigbar4(i,:)./abs(JeidotFP(i,:)) + Sigtil4(i,:)./abs(JeidotR(i,:)))./(1./abs(JeidotFP(i,:))+1./abs(JeidotR(i,:)));

        % Plot and store solid line for legend
        h = plot(alpha_arr(i,:), norm_sum, ...
            'LineWidth', 2, ...
            'Color', cmap(curve_idx,:));
        hold on

        legend_handles(curve_idx) = h;
        curve_idx = curve_idx + 1;
    end
end

ylim([-3.5*10^-3 0])

xlabel('$\alpha$', 'interpreter','latex','FontSize',18)
ylabel('$\sigma_{\mathrm{IV}}$', ...
       'interpreter','latex','FontSize',22)

legend_labels = arrayfun(@(x) sprintf('$%d$', x*5), selected_taus, 'UniformOutput', false);
lgd = legend(legend_handles, legend_labels, 'Interpreter','latex','Location','SouthEast');
lgd.Title.String = '$\tau_{E,-} \ [\mathrm{ms}]$';

set(gca,'FontSize',13)
set(gca,'TickLabelInterpreter','latex')
grid on