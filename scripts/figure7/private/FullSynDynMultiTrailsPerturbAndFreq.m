function [stdmebar,stdmibar,stdmetil,stdmitil,Jeimeandyn,Jiemeandyn,frequency,m_e_bef_per,m_i_bef_per,m_e_aft_per,m_i_aft_per,timem,JbarDyesint,Jeibef,Jiebef,Jeiaft,Jieaft,Jeiper,Jieper] = FullSynDynMultiTrailsPerturbAndFreq(m_e_history,m_i_history,D,Jei,Jie,Jee,Jii,t,Jiemax,mu,alpha,JbarD,Jeemean,Jiimean,vn,ve,vu,vaei,vaie,vdei,vdie,mulvn,mulve,mulva,mulvd,dtlearn,lambda_e,lambda_i,K_Ibar,K_Itilmphi,K_Itilphi,K_Eptilmpsi,K_Etilphi,K_Etilmphi,K_pI,K_pE,K_mI,K_mE,N_e,N_i,dt,tf,num_dead,death_time)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

stdmebar=nan(1,length(t));
stdmibar=nan(1,length(t));
stdmetil=nan(1,length(t));
stdmitil=nan(1,length(t));

Jeimeandyn=nan(1,length(t));
Jiemeandyn=nan(1,length(t));

frequency=nan(1,length(t));

flagafter=0; % this flag becomes 1 ones it is 10 time steps closer to the end of the simulation and inside the R region!
Jiimeaninit=Jiimean;
Jeemeaninit=Jeemean;
comp=(N_i-num_dead)/N_i;

syms wDyesint JbarDyesint
if Jeemean>=Jiimean
    range=[0.1 5 ;0.01 pi/(2*D)];
    Y=vpasolve([(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5==1/cos(wDyesint*D-acos((Jeemean-Jiimean*comp)/(2*(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5))), wDyesint==-tan(wDyesint*D-acos((Jeemean-Jiimean*comp)/(2*(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5)))], [JbarDyesint,wDyesint],range); % the frequency on the bifurcation line
elseif Jeemean<Jiimean
    range=[0.1 5 ;0.01 pi/D];
    Y=vpasolve([(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5==1/cos(wDyesint*D-acos((Jeemean-Jiimean*comp)/(2*(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5))), wDyesint==-tan(wDyesint*D-acos((Jeemean-Jiimean*comp)/(2*(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5)))], [JbarDyesint,wDyesint],range); % the frequency on the bifurcation line
end
JbarDyesint=double(Y.JbarDyesint);

Jeibef=0;
Jiebef=0;

for i=1:length(t)
    
    if i>=death_time
        Jii(1:num_dead,:)=0;
        Jii(:,1:num_dead)=0;
        Jei(:,1:num_dead)=0;
        Jie(1:num_dead,:)=0;
        JbarD=JbarDyesint;
        if i==death_time
           Jeiper=mean(nonzeros(Jei(:)));
           Jieper=mean(nonzeros(Jie(:)));
        end
    end
    
    
    %%% construcing the eigenvectors on each time step %%%
    Jeimean=mean(nonzeros(Jei(:))); % I do not want the dead synapses to participate in the dynamics
    Jiemean=mean(nonzeros(Jie(:)));
    
    Jeimeandyn(i)=Jeimean;
    Jiemeandyn(i)=Jiemean;
    
    if Jeimean*Jiemean<JbarD^2 && i<death_time % in the FP regime we solve exactly
        J=[[Jee/N_e -Jei/N_i];[Jie/N_e -Jii/N_i]]; % connectivity matrix
        diagonal=eye(N_e+N_i);
        mfa=(diagonal-J)\ones(N_e+N_i,1); % The exact fixed point solution
        stdmebar(i)=std(mfa(1:N_e));
        if i<death_time
            stdmibar(i)=std(mfa(N_e+1:N_e+N_i));
%             if i==death_time-1
%                 [m_e_bef_per,m_i_bef_per,~,~,timem] = Two_populations_full_rate_model_history(m_e_history,m_i_history,Jee,Jei,Jie,Jii,dt,tf);
%             end
        else
            stdmibar(i)=std(mfa(N_e+1+num_dead:N_e+N_i));
        end
        stdmetil(i)=nan;
        stdmitil(i)=nan;
        Corr_ei=mfa(1:N_e)*mfa((N_e+1):(N_e+N_i))';
        Corr_ie=Corr_ei';
        J_ei_dot=lambda_i*((1-alpha)*Corr_ei); % NOTICE that I put Corr_ei instead of Corr_ie because there is symmetry Cei,ij=Cie,ji in the FP regime
        J_ie_dot=lambda_e*(((Th_li_full(1-Jie/Jiemax)).^mu-alpha).*Corr_ie); % NOTICE that I put Corr_ei instead of Corr_ie because there is symmetry Cei,ij=Cie,ji in the FP regime
    elseif Jeimean*Jiemean<JbarD^2 && i>=death_time 
        J=[[Jee/N_e -Jei/(N_i-num_dead)];[Jie/N_e -Jii/(N_i-num_dead)]]; % connectivity matrix
        diagonal=eye(N_e+N_i);
        mfa=(diagonal-J)\ones(N_e+N_i,1); % The exact fixed point solution
        stdmebar(i)=std(mfa(1:N_e));
        stdmibar(i)=std(mfa(N_e+1+num_dead:N_e+N_i));
        stdmetil(i)=nan;
        stdmitil(i)=nan;
        Corr_ei=mfa(1:N_e)*mfa((N_e+1):(N_e+N_i))';
        Corr_ie=Corr_ei';
        J_ei_dot=lambda_i*((1-alpha)*Corr_ei); % NOTICE that I put Corr_ei instead of Corr_ie because there is symmetry Cei,ij=Cie,ji in the FP regime
        J_ie_dot=lambda_e*(((Th_li_full(1-Jie/Jiemax)).^mu-alpha).*Corr_ie); % NOTICE that I put Corr_ei instead of Corr_ie because there is symmetry Cei,ij=Cie,ji in the FP regime
    else % The R regime is solved with simulations
        [Corr_ie,Corr_ei,~,~,m_e_T,m_i_T,~,~,~]=proj.common.Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf);
        if i<death_time && i>death_time-10 && Jeimean*Jiemean>Jeibef*Jiebef 
            [m_e_bef_per,m_i_bef_per,~,~,timem] = proj.common.Two_populations_full_rate_model_history(m_e_history,m_i_history,Jee,Jei,Jie,Jii,dt,tf);
            Jeibef=mean(Jei(:));
            Jiebef=mean(Jie(:));
        elseif i>=death_time
            m_i_T=m_i_T(num_dead+1:end,:);
            if i>length(t)-5 && ~flagafter
                [m_e_aft_per,m_i_aft_per,~,~,timem] = proj.common.Two_populations_full_rate_model_history(m_e_history,m_i_history,Jee,Jei,Jie,Jii,dt,tf);
                Jeiaft=mean(nonzeros(Jei(:)));
                Jieaft=mean(nonzeros(Jie(:)));
                flagafter=1;
            end
        end
        [~,~,T,~,~] = proj.common.Two_populations_full_rate_model_history(m_e_history,m_i_history,Jee,Jei,Jie,Jii,dt,tf);
        frequency(i)=1/mean(T);
        stdmebar(i)=std(mean(m_e_T,2));
        stdmibar(i)=std(mean(m_i_T,2));
        stdmetil(i)=std(max(m_e_T,[],2)-mean(m_e_T,2));
        stdmitil(i)=std(max(m_i_T,[],2)-mean(m_i_T,2));
        for n=1:N_e
            for k=1:N_i
                J_ei_dot(n,k)=lambda_i*((K_pI-alpha*K_mI)*squeeze(Corr_ie(k,n,:))*dt);
                J_ie_dot(k,n)=lambda_e*(((Th_li_full(1-Jie(k,n)/Jiemax))^mu*K_pE-alpha*K_mE)*squeeze(Corr_ei(n,k,:))*dt);
            end
        end
    end
    
    noiseI=0.01;%0.0003;
    noiseE=10*noiseI;
    Jei=(Jei+dtlearn*(J_ei_dot+noiseI*randn(size(Jei)))).*heaviside(Jei+dtlearn*(J_ei_dot+noiseI*randn(size(Jei))));
    Jie=(Jie+dtlearn*(J_ie_dot+noiseE*randn(size(Jie)))).*heaviside(Jie+dtlearn*(J_ie_dot+noiseE*randn(size(Jie))));
    
end
end