%% The MEAN synaptic dynamics of inter and intra %%
%%% plot of the phase diagram in
T=1;
dt=0.01;
Tunits=5*10^-3; % 5ms for T=1
tf=200;
D=0.4; % delay in msec
Jeimean=0.2;%0.5304;%0.32;%0.4;
Jiemean=1;%4.1907;%10;%5;%JbarD^2/Jeimean+0.1;%10^-5;
Jiimean=0;
Jeemean=0.2;%0.8635;
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

alpha=0.98; % relative depression %%
muie=0.01; % measure of linearity
muee=0.02;
Jiemax=20; % J_ie_max
Jeemax=1;
Jiimax=0;
tau_pE=2; % typical potentiation time of excitatory synapses
tau_pI=2; % typical potentiation time of inhibitory synapses
tau_mE=5; % typical depression time of excitatory synapses
tau_mI=3; % typical depression time of inhibitory synapses
thetapI=acos((1+(wD*tau_pI)^2)^-0.5);
thetamI=acos((1+(wD*tau_mI)^2)^-0.5);
thetapE=acos((1+(wD*tau_pE)^2)^-0.5);
thetamE=acos((1+(wD*tau_mE)^2)^-0.5);
lambdascale=3;%120;
lambda_ie=lambdascale*1;%0*200*15/1000; % excitatory synapses learning rate
lambda_ei=lambdascale*0.1; % inhibitory synapses learning rate
lambda_ee=lambdascale*0.1;
H_E=-1;
H_I=1;
%%% zero order synapses %%%
K_Ibar=1-alpha;

K_Iptilpsi=cos(thetapI)*cos(thetapI+psi);
K_Imtilpsi=cos(thetamI)*cos(thetamI-psi);
K_Itilpsi=K_Iptilpsi-alpha*K_Imtilpsi;

gabs=(JbarD^2-Jeemean*Jiimean)^0.5;
J=[[Jeemean/N_e -Jeimean/N_i];[Jiemean/N_e -Jiimean/N_i]]; % connectivity matrix
diagonal=eye(N_e+N_i);
mfa=(diagonal-J)\ones(N_e+N_i,1); % The exact fixed point solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Corr_ie_final,Corr_ei_final,Corr_ee_final,Corr_ii_final,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended]=proj.common.Correlations_2D_full_diff(Jeemean,Jeimean,Jiemean,Jiimean,dt,tf);

t=1:11500;
%check_arr=[reshape(check_ei.',1,[]) reshape(check_ie.',1,[])];
dtlearn=1;

Jeimeanarr=nan(1,length(t));
Jiemeanarr=nan(1,length(t));
Jeemeanarr=nan(1,length(t));
Jiimeanarr=nan(1,length(t));
JbarDarr=nan(1,length(t));

m_e_history=0.2;
m_i_history=1.3;

KpI=1/tau_pI*exp(-Delta_extended*H_I/tau_pI).*heaviside(H_I*Delta_extended);
KmI=1/tau_mI*exp(Delta_extended*H_I/tau_mI).*heaviside(-H_I*Delta_extended);
KpE=1/tau_pE*exp(-Delta_extended*H_E/tau_pE).*heaviside(H_E*Delta_extended);
KmE=1/tau_mE*exp(Delta_extended*H_E/tau_mE).*heaviside(-H_E*Delta_extended);

frequencyarr=nan(1,length(t));

KEI=KpI-alpha*KmI;
figure;
for i=2:length(t)

    %%% construcing the eigenvectors on each time step %%%

    fie=(1-Jiemean/Jiemax)^muie;
    fee=(1-Jeemean/Jeemax)^muee;
    %fii=(1-Jiimean/Jiimax)^mu;


    if i==2
        PhDiagJeeJbar=figure(1);
        plot((Jeimean*Jiemean)^0.5,Jeemean,'.','Color',[0.2 0.8 0.2])
        hold on
        T=1;
        D=0.4;
        Jii_init=0;
        dJii=0.2;
        Jii_final=0;
        Jii_arr=Jii_init:dJii:Jii_final;
        Jee_init=0;
        dJee=0.05;
        Jee_final=5;
        Jee_arr=(Jee_init:dJee:Jee_final)+0.001;%(2+Jii_init-0.05);
        JbarD_arr=nan(length(Jii_arr),length(Jee_arr));
        wD_arr=nan(length(Jii_arr),length(Jee_arr));
        for k=1:length(Jii_arr)
            Jii=Jii_arr(k);
            for j=1:length(Jee_arr)
                Jee=Jee_arr(j);
                syms wD JbarD
                if Jee>=2+Jii
                    break
                elseif Jee>Jii
                    range=[0.1 5 ;0.2 pi/(2*D)];
                    Y=vpasolve([(JbarD^2-Jee*Jii)^0.5==1/cos(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5))), T*wD==-tan(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
                    JbarD_arr(k,j)=double(Y.JbarD);
                    wD_arr(k,j)=double(Y.wD);
                elseif Jee<Jii
                    range=[0.1 5 ;0.2 pi/D];
                    Y=vpasolve([(JbarD^2-Jee*Jii)^0.5==1/cos(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5))), T*wD==-tan(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
                    JbarD_arr(k,j)=double(Y.JbarD);
                    wD_arr(k,j)=double(Y.wD);
                elseif Jee==Jii
                    range=[0.1 pi/(2*D)-0.01];
                    Y=vpasolve(T*wD==cot(wD*D),wD,range); % the frequency on the bifurcation line
                    JbarD_arr(k,j)=(1+(double(Y)*T)^2)^0.5;
                    wD_arr(k,j)=double(Y);
                end
            end
            JbarD_arr(k,j)=1+Jii;
            JbarD_arrlast=JbarD_arr(end,:);
            Jbar_arrlast=0:0.004:4;
            f=nan(length(Jbar_arrlast),length(Jee_arr));
            parfor w=1:length(Jee_arr)
                Jee=Jee_arr(w);
                JbarD=JbarD_arrlast(w);
                if Jee<=2+Jii
                    f(:,w)=proj.common.FrequencyContourJbarJee(JbarD,Jbar_arrlast,Jee,Jii,dt,tf);
                else
                    f(:,w)=proj.common.FrequencyContourJbarJee(nan,Jbar_arrlast,Jee,Jii,dt,tf);
                end
            end
            f=f/(Tunits); % units in Hz
            [JBAR_arrlast,JEE_arr]=meshgrid(Jbar_arrlast,Jee_arr);
            [C,h]=contourf(JBAR_arrlast,JEE_arr,f',100, 'HandleVisibility', 'off');
            hold on
            set(h,'LineColor','none')
            plot(JbarD_arr(k,:),Jee_arr,'Color','black','LineWidth',3,'Color',[(k-1)/(length(Jii_arr)+1) (k-1)/(length(Jii_arr)+1) (k-1)/(length(Jii_arr)+1)])
            plot(0:0.01:(1+Jii),1+(0:0.01:(1+Jii)).^2/(1+Jii),'LineWidth',3,'Color',[(k-1)/(length(Jii_arr)+1) (k-1)/(length(Jii_arr)+1) (k-1)/(length(Jii_arr)+1)], 'HandleVisibility', 'off')
            plot(((1+Jii):0.1:4),2*((1+Jii):0.1:4)-Jii,'LineWidth',3,'Color',[(k-1)/(length(Jii_arr)+1) (k-1)/(length(Jii_arr)+1) (k-1)/(length(Jii_arr)+1)], 'HandleVisibility', 'off')
            text(0.8,0.4,'$\mathrm{FP}$','interpreter','latex','FontSize',24)
            text(3,3,'$\mathrm{R}$','interpreter','latex','FontSize',24)
            text(0.6,2.5,'$\mathrm{D}$','interpreter','latex','FontSize',24)
            grid on
            cb = colorbar;                   % create colorbar
            % Get the position of the colorbar in normalized units
            cbPos = cb.Position;
            
            % Define title text and position it above the colorbar
            titleStr = '$\mathrm{f} $ [Hz]';  % LaTeX string
            
            % Create a text object in the same axes (parent of colorbar)
            ax = gca;
            xPos = cbPos(1) + cbPos(3)/2;     % center of colorbar
            yPos = cbPos(2) + cbPos(4) + 0.02; % slightly above the colorbar
            
            % Add the title as a text object
            annotation('textbox', [xPos-0.05, yPos, 0.1, 0.05], ...
                'String', titleStr, ...
                'Interpreter', 'latex', ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 12);
        end
        
        xlabel('$\bar{J}$','interpreter','latex','FontSize',18)
        ylabel('$J_{EE}$','interpreter','latex','FontSize',18)
        xlim([0 4])
        ylim([0 5])
        set(gca,'FontSize',14)
        set(gca,'TickLabelInterpreter','latex')
        grid on

        %JbarDarrNullofJee=linspace(nanmin(JbarD_arr(:)),nanmax(JbarD_arr(:)),30);
        %Hm=1+1./(2*JbarDarrNullofJee.^2.*(1+(JbarDarrNullofJee.^2-1)*tau_mE^2));
        %Hp=1+1./(2*JbarDarrNullofJee.^2.*(1+(JbarDarrNullofJee.^2-1)*tau_pE^2));

        JeenullFP=Jeemax*(1-alpha^(1/muee));
        %JeenullR=Jeemax*(1-(alpha*Hm./Hp).^(1/muee));
        %plot(0:0.2:1.5,JeenullFP*ones(1,length(0:0.2:1.5)),'LineWidth',3,'Color','Red')
        %plot(JbarDarrNullofJee,JeenullR,'LineWidth',3,'Color','Red')
        saveas(PhDiagJeeJbar,'PhDiagJeeJbar.fig')
    end


    figure(1)
    if i==2
        h1=plot((Jeimean*Jiemean)^0.5,Jeemean,'+','Color',[0.8500, 0.3250, 0.0980],'MarkerSize', 12,'LineWidth',2,'HandleVisibility', 'off'); % No legend entry
    elseif i==length(t)
        h2=plot((Jeimean*Jiemean)^0.5,Jeemean,'x','Color',[1, 1, 1],'MarkerSize', 12,'LineWidth',2,'HandleVisibility', 'off'); % No legend entry
    else
        h3=plot((Jeimean*Jiemean)^0.5,Jeemean,'.','Color',[0, 0.4470, 0.7410], 'HandleVisibility', 'off'); % No legend entry
    end
        syms wD JbarD
    if Jeemean>=Jiimean
        range=[0.1 5 ;0.01 pi/(2*D)];
        Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), T*wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
    elseif Jeemean<Jiimean
        range=[0.1 5 ;0.01 pi/D];
        Y=vpasolve([(JbarD^2-Jeemean*Jiimean)^0.5==1/cos(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5))), T*wD==-tan(wD*D-acos((Jeemean-Jiimean)/(2*(JbarD^2-Jeemean*Jiimean)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
    end
    wD=double(Y.wD);
    JbarD=double(Y.JbarD);
    thetapI=acos((1+(wD*tau_pI)^2)^-0.5);
    thetamI=acos((1+(wD*tau_mI)^2)^-0.5);
    K_Ibar=1-alpha;

    K_Iptilpsi=cos(thetapI)*cos(thetapI+psi);
    K_Imtilpsi=cos(thetamI)*cos(thetamI-psi);
    K_Itilpsi=K_Iptilpsi-alpha*K_Imtilpsi;


    K_Eptilpsi=cos(thetapE)*cos(thetapE+psi);
    K_Emtilpsi=cos(thetamE)*cos(thetamE-psi);

    Jhat=(Jeemean+Jiimean)/2;
    psi=acos(Jhat/JbarD);


    Jeimeanarr(i)=Jeimean;
    Jiemeanarr(i)=Jiemean;
    Jeemeanarr(i)=Jeemean;
    Jiimeanarr(i)=Jiimean;
    JbarDarr(i)=JbarD;





    if Jeimean*Jiemean<double(JbarD^2)
        me=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
        mi=(1+Jiemean-Jeemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
        Jeidot=lambda_ei*(1-alpha)*me*mi;
        Jiedot=lambda_ie*(fie-alpha)*me*mi;
        Jeedot=lambda_ee*(fee-alpha)*me^2;
        %Jiidot=lambda_i*(fii-alpha)*mi*mi;
    else
        [Corr_ie_final,Corr_ei_final,Corr_ee_final,Corr_ii_final,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended]=proj.common.Correlations_2D_full_diff(Jeemean,Jeimean,Jiemean,Jiimean,dt,tf);
        frequencyarr(i)=1/(T_mean_m_e*Tunits);
        KIE=fie*KpE-alpha*KmE;
        KEE=fee*KpE-alpha*KmE;
        %KII=fii*KpI-alpha*KmI;
        Jeidot=lambda_ei*dt*KEI*squeeze(Corr_ie_final)*dt;
        Jiedot=lambda_ie*dt*KIE*squeeze(Corr_ei_final)*dt;
        Jeedot=lambda_ee*dt*KEE*Corr_ee_final'*dt;
        %Jiidot=lambda_i*dt*KII'*squeeze(Corr_ii_final)*dt
    end


    if (Jeidot^2+Jiedot^2+Jeedot^2)^0.5<0.01
        lambdascale=lambdascale*2;%120;
        lambda_ie=lambdascale*1;%0*200*15/1000; % excitatory synapses learning rate
        lambda_ei=lambdascale*0.1; % inhibitory synapses learning rate
        lambda_ee=lambdascale*0.1;
        if Jeimean*Jiemean<double(JbarD^2)
            me=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
            mi=(1+Jiemean-Jeemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
            Jeidot=lambda_ei*(1-alpha)*me*mi;
            Jiedot=lambda_ie*(fie-alpha)*me*mi;
            Jeedot=lambda_ee*(fee-alpha)*me^2;
            %Jiidot=lambda_i*(fii-alpha)*mi*mi;
        else
            [Corr_ie_final,Corr_ei_final,Corr_ee_final,Corr_ii_final,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended]=proj.common.Correlations_2D_full_diff(Jeemean,Jeimean,Jiemean,Jiimean,dt,tf);
            frequencyarr(i)=1/(T_mean_m_e*Tunits);
%             if frequencyarr(i)>0 && frequencyarr(i)<5 
%                frequencyarr(i)=nan;  
%             end
            KIE=fie*KpE-alpha*KmE;
            KEE=fee*KpE-alpha*KmE;
            %KII=fii*KpI-alpha*KmI;
            Jeidot=lambda_ei*dt*KEI*squeeze(Corr_ie_final)*dt;
            Jiedot=lambda_ie*dt*KIE*squeeze(Corr_ei_final)*dt;
            Jeedot=lambda_ee*dt*KEE*Corr_ee_final'*dt;
            %Jiidot=lambda_i*dt*KII'*squeeze(Corr_ii_final)*dt
        end
    elseif (Jeidot^2+Jiedot^2+Jeedot^2)^0.5>0.1
        lambdascale=3;%120;
        lambda_ie=lambdascale*1;%0*200*15/1000; % excitatory synapses learning rate
        lambda_ei=lambdascale*0.1; % inhibitory synapses learning rate
        lambda_ee=lambdascale*0.1;
        if Jeimean*Jiemean<double(JbarD^2)
            me=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
            mi=(1+Jiemean-Jeemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean);
            Jeidot=lambda_ei*(1-alpha)*me*mi;
            Jiedot=lambda_ie*(fie-alpha)*me*mi;
            Jeedot=lambda_ee*(fee-alpha)*me^2;
            %Jiidot=lambda_i*(fii-alpha)*mi*mi;
        else
            [Corr_ie_final,Corr_ei_final,Corr_ee_final,Corr_ii_final,m_e_T,m_i_T,T_mean_m_e,T_mean_m_i,Delta_extended]=proj.common.Correlations_2D_full_diff(Jeemean,Jeimean,Jiemean,Jiimean,dt,tf);
            frequencyarr(i)=1/(T_mean_m_e*Tunits);
%             if frequencyarr(i)>0 && frequencyarr(i)<5
%                 frequencyarr(i)=nan;
%             end
            KIE=fie*KpE-alpha*KmE;
            KEE=fee*KpE-alpha*KmE;
            %KII=fii*KpI-alpha*KmI;
            Jeidot=lambda_ei*dt*KEI*squeeze(Corr_ie_final)*dt;
            Jiedot=lambda_ie*dt*KIE*squeeze(Corr_ei_final)*dt;
            Jeedot=lambda_ee*dt*KEE*Corr_ee_final'*dt;
            %Jiidot=lambda_i*dt*KII'*squeeze(Corr_ii_final)*dt
        end
    end

    Jeimean=Jeimean+dtlearn*Jeidot;
    Jiemean=Jiemean+dtlearn*Jiedot;
    Jeemean=Jeemean+dtlearn*Jeedot;
    %Jiimean=Jiimean+dtlearn*Jiidot;
    
    Jei_init=0;
    dJei=0.01;
    Jei_final=1+Jiimean;
    Jei_arr=Jei_init:dJei:Jei_final;


    drawnow

    if mod(i,1000)==0
        save("data\JeeJeiJiedynamics","Jeemeanarr","Jeimeanarr","Jiemeanarr","Jii","frequencyarr")
        %close all
        %loadfig('PhDiagJeeJbar.fig')
    end

end
%%
fig5=figure(5)
t=axes;
plot(t,frequencyarr,'.','Color',[0, 0.4470, 0.7410]	, 'HandleVisibility', 'off'); % No legend entry


% Create inset axes in fig1
figure(1);  % Make sure you're working in the main figure
insetAx = axes('Position', [0.18 0.65 0.18 0.18]);  % [left, bottom, width, height]

% Copy all children (i.e., lines, titles, etc.) from ax2 into insetAx
copyobj(allchild(t), insetAx);

% adjust inset axes limits, etc.
xlim(insetAx, t.XLim);
ylim(insetAx, t.YLim);
xlabel(insetAx,'$t \ [\mathrm{a.u.}]$','interpreter','latex','FontSize',12);
ylabel('$\mathrm{f} \ [\mathrm{Hz}]$','interpreter','latex','FontSize',12)
grid on
% Close the second figure if no longer needed
close(fig5);