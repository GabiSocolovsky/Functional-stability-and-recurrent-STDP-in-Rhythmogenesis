%% Jee Jbar phase diagram %%
T=1;
Tunits  = 5e-3;      % Time unit = Tau_m = 5 ms
dt=0.01;
tf=200;
D=0.4;
Jii_init=0;
dJii=0.2;
Jii_final=0.8;
Jii_arr=Jii_init:dJii:Jii_final;
Jee_init=0;
dJee=0.05;
Jee_final=5;
Jee_arr=(Jee_init:dJee:Jee_final)+0.001;%(2+Jii_init-0.05);
JbarD_arr=nan(length(Jii_arr),length(Jee_arr));
wD_arr=nan(length(Jii_arr),length(Jee_arr));
for i=1:length(Jii_arr)
    Jii=Jii_arr(i);
    for j=1:length(Jee_arr)
        Jee=Jee_arr(j);
        syms wD JbarD
        if Jee>=2+Jii
            break        
        elseif Jee>Jii
            range=[0.1 5 ;0.2 pi/(2*D)];
            Y=vpasolve([(JbarD^2-Jee*Jii)^0.5==1/cos(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5))), T*wD==-tan(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
            JbarD_arr(i,j)=double(Y.JbarD);
            wD_arr(i,j)=double(Y.wD);
        elseif Jee<Jii
            range=[0.1 5 ;0.2 pi/D];    
            Y=vpasolve([(JbarD^2-Jee*Jii)^0.5==1/cos(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5))), T*wD==-tan(wD*D-acos((Jee-Jii)/(2*(JbarD^2-Jee*Jii)^0.5)))], [JbarD,wD],range); % the frequency on the bifurcation line
            JbarD_arr(i,j)=double(Y.JbarD);
            wD_arr(i,j)=double(Y.wD);
        elseif Jee==Jii
            range=[0.1 pi/(2*D)-0.01];
            Y=vpasolve(T*wD==cot(wD*D),wD,range); % the frequency on the bifurcation line            
            JbarD_arr(i,j)=(1+(double(Y)*T)^2)^0.5;
            wD_arr(i,j)=double(Y);
        end
    end
    JbarD_arr(i,j)=1+Jii;

    %%%% Frequency Contour %%%%
    if Jii==0.8
        JbarD_arrlast=JbarD_arr(end,:);
        Jbar_arrlast=0:0.002:4;
        f=nan(length(Jbar_arrlast),length(Jee_arr));
        parfor k=1:length(Jee_arr)
            Jee=Jee_arr(k);
            JbarD=JbarD_arrlast(k);
            if Jee<=2+Jii
                f(:,k)=proj.common.FrequencyContourJbarJee(JbarD,Jbar_arrlast,Jee,Jii,dt,tf);
            else
                f(:,k)=proj.common.FrequencyContourJbarJee(nan,Jbar_arrlast,Jee,Jii,dt,tf);
            end
        end
            f=f/(Tunits); % units in Hz
            [JBAR_arrlast,JEE_arr]=meshgrid(Jbar_arrlast,Jee_arr);
        [C,h]=contourf(JBAR_arrlast,JEE_arr,f',100, 'HandleVisibility', 'off');
        set(h,'LineColor','none')
        grid on
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(JbarD_arr(i,:),Jee_arr,'Color','black','LineWidth',3,'Color',[(i-1)/(length(Jii_arr)+1) (i-1)/(length(Jii_arr)+1) (i-1)/(length(Jii_arr)+1)])
    hold on
    plot(0:0.01:(1+Jii),1+(0:0.01:(1+Jii)).^2/(1+Jii),'LineWidth',3,'Color',[(i-1)/(length(Jii_arr)+1) (i-1)/(length(Jii_arr)+1) (i-1)/(length(Jii_arr)+1)], 'HandleVisibility', 'off')
    plot(((1+Jii):0.1:4),2*((1+Jii):0.1:4)-Jii,'LineWidth',3,'Color',[(i-1)/(length(Jii_arr)+1) (i-1)/(length(Jii_arr)+1) (i-1)/(length(Jii_arr)+1)], 'HandleVisibility', 'off')
end
plot(0:0.5:4,ones(1,length(0:0.5:4)),'--','Color','red', 'HandleVisibility', 'off')
lgd=legend({'$0$','$0.2$','$ 0.4$','$0.6$','$0.8$'},'Interpreter','latex');
legendJiiTitle=sprintf('$J_{II}$','Interpreter','latex');
lgd.Title.String = legendJiiTitle;
xlabel('$\bar{J}$','interpreter','latex','FontSize',18)
ylabel('$J_{EE}$','interpreter','latex','FontSize',18)
text(0.8,0.4,'$\mathrm{FP}$','interpreter','latex','FontSize',24)
text(3,3,'$\mathrm{R}$','interpreter','latex','FontSize',24)
text(0.6,3.5,'$\mathrm{D}$','interpreter','latex','FontSize',24)
xlim([0 4])
ylim([0 5])
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')