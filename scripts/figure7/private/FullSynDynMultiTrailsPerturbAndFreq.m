%% Full STDP dynamics with death of neurons multi trials %%
function [stdmebar,stdmibar,stdmetil,stdmitil,Jeimeandyn,Jiemeandyn,frequency,m_e_bef_per,m_i_bef_per,m_e_aft_per,m_i_aft_per,timem,JbarDyesint,Jeibef,Jiebef,Jeiaft,Jieaft,Jeiper,Jieper] = FullSynDynMultiTrailsPerturbAndFreq(m_e_history,m_i_history,D,Jei,Jie,Jee,Jii,t,Jiemax,mu,alpha,JbarD,Jeemean,Jiimean,vn,ve,vu,vaei,vaie,vdei,vdie,mulvn,mulve,mulva,mulvd,dtlearn,lambda_e,lambda_i,K_Ibar,K_Itilmphi,K_Itilphi,K_Eptilmpsi,K_Etilphi,K_Etilmphi,K_pI,K_pE,K_mI,K_mE,N_e,N_i,dt,tf,num_dead,death_time)
% This function computes the full STDP with death of neurons in the paper:
% "Functional stability and recurrent STDP in Rhythmogenesis"

%   Description:

%      - Computes the STDP dynamics of a network that at a certain time
%      loses 3 neurons
%      - Computes neural firing rates before and after the death of the neurons
%      - Computes the frequency of the system at each STDP time step (along
%      with the standard deviation)

%   Inputs:

%       m_e_history, m_i_history - history of each population
%       D - time delay
%       t - simulation time
%       Jiemax, alpha, mu - defined in paper
%       JbarD - Jbar on bif. line
%       Jeemean - mean value of Jee
%       Jiimean - mean value of Jii
%       vn - set of eigenvectors of family I
%       ve - set of eigenvectors of family II
%       vu - set of eigenvectots of uniform family
%       vaei - Jei eigenvectors set of family III
%       vaie - Jie eigenvectors set of family III
%       vdei - Jei eigenvectors set of family IV
%       vdie - Jie eigenvectors set of family IV
%       mulvn - size of family I
%       mulve - size of family II
%       mulva - size of family III
%       mulvd - size of family IV
%       dtlearn - time bin of STDP dynamics
%       lambda_e, lambda_i - learning rates
%       mu - non-linearity
%       alpha - relative strength of depression
%       K_Ibar,K_Itilmphi,K_Itilphi,K_Eptilmpsi,K_Etilphi,K_Etilmphi,K_pI,K_pE,K_mI,K_mE
%       - temporal kernels and analytical expression (see paper)
%       N_e,N_i - Size of populations
%       Jee  -   E-to-E connectivity matrix
%       Jei  -   I-to-E connectivity matrix
%       Jie  -   E-to-I connectivity matrix
%       Jii  -   I-to-I connectivity matrix
%       dt   -   time bin
%       tf   -   final time of simulation
%       num_dead   -   number of dead neurons
%       death_time  -  time when killing neurons

%   Outputs:

%       stdmebar, stdmibar - Standard deviation across population
%       of m_x bar (x=e,i)
%       stdmetil,stdmitil - Standard deviation across population
%       of m_x tilde (x=e,i)
%       Jeimeandyn, Jiemeandyn - order parameters of Jei and Jie
%       dJ_arr_t - The fluctuations of the synaptic weights Jei,ij and Jie,ij
%       frequency - The frequency of the network at each STDP time step
%       m_e_bef_per - m_E right before killing the neruons
%       m_i_bef_per - m_I right before killing the neruons
%       m_e_aft_per - m_E after killing the neruons (when in CR)
%       m_i_aft_per - m_I after killing the neruons (when in CR)
%       timem - time of network dynamics
%       JbarDyesint - Jbar on bif. after killing neurons
%       Jeibef - Jei before death of neurons
%       Jiebef - Jie before death of neurons
%       Jeiaft - Jei after death of neurons
%       Jieaft - Jei after death of neurons
%       Jeiper - Jei at the time of death of neurons
%       Jieper - Jie at the time of death of neurons

%   Dependencies:

%   Correlations_2D_full_diff.m - Computes cross-correlations
%   Two_populations_full_rate_model_history - Computes the firing rate
%   dynamics of the network
%
% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Defining variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definitions found in introduction %
stdmebar=nan(1,length(t));
stdmibar=nan(1,length(t));
stdmetil=nan(1,length(t));
stdmitil=nan(1,length(t));

Jeimeandyn=nan(1,length(t));
Jiemeandyn=nan(1,length(t));

frequency=nan(1,length(t));

flagafter=0; % this flag becomes 1 ones it is 10 time steps closer to the end of the simulation and inside the R region!
comp=(N_i-num_dead)/N_i; % compementary factor after neurons die

syms wDyesint JbarDyesint
if Jeemean>=Jiimean
    range=[0.1 5 ;0.01 pi/(2*D)];
    Y=vpasolve([(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5==1/cos(wDyesint*D-acos((Jeemean-Jiimean*comp)/(2*(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5))), wDyesint==-tan(wDyesint*D-acos((Jeemean-Jiimean*comp)/(2*(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5)))], [JbarDyesint,wDyesint],range); % the frequency on the bifurcation line
elseif Jeemean<Jiimean
    range=[0.1 5 ;0.01 pi/D];
    Y=vpasolve([(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5==1/cos(wDyesint*D-acos((Jeemean-Jiimean*comp)/(2*(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5))), wDyesint==-tan(wDyesint*D-acos((Jeemean-Jiimean*comp)/(2*(JbarDyesint^2*comp-Jeemean*Jiimean*comp)^0.5)))], [JbarDyesint,wDyesint],range); % the frequency on the bifurcation line
end
JbarDyesint=double(Y.JbarDyesint); % Jbar on bifurcation line after death of neurons

Jeibef=0;
Jiebef=0;

for i=1:length(t)
    
    if i>=death_time % after death_time the connections to 3 neurons are dead
        Jii(1:num_dead,:)=0;
        Jii(:,1:num_dead)=0;
        Jei(:,1:num_dead)=0;
        Jie(1:num_dead,:)=0;
        JbarD=JbarDyesint; % this is the new bif. line
        if i==death_time % at time of death what are the synaptic connections?
            Jeiper=mean(nonzeros(Jei(:)));
            Jieper=mean(nonzeros(Jie(:)));
        end
    end
    
    
    %%% construcing the eigenvectors on each time step %%%
    Jeimean=mean(nonzeros(Jei(:))); % I do not want the dead synapses to participate in the dynamics
    Jiemean=mean(nonzeros(Jie(:)));
    
    Jeimeandyn(i)=Jeimean; % dynamics of order parameters of Jei
    Jiemeandyn(i)=Jiemean; % and Jie
    
    if Jeimean*Jiemean<JbarD^2 && i<death_time % in the FP regime we solve exactly
        J=[[Jee/N_e -Jei/N_i];[Jie/N_e -Jii/N_i]]; % connectivity matrix
        diagonal=eye(N_e+N_i);
        mfa=(diagonal-J)\ones(N_e+N_i,1); % The exact fixed point solution
        stdmebar(i)=std(mfa(1:N_e));
        if i<death_time % compute the SD of mi bar differently before and after death
            stdmibar(i)=std(mfa(N_e+1:N_e+N_i));
        else
            stdmibar(i)=std(mfa(N_e+1+num_dead:N_e+N_i));
        end
        stdmetil(i)=nan;
        stdmitil(i)=nan;
        Corr_ei=mfa(1:N_e)*mfa((N_e+1):(N_e+N_i))'; % cross-correlations
        Corr_ie=Corr_ei'; % cross-correlations
        J_ei_dot=lambda_i*((1-alpha)*Corr_ei); % NOTICE that I put Corr_ei instead of Corr_ie because there is symmetry Cei,ij=Cie,ji in the FP regime
        J_ie_dot=lambda_e*(((Th_li_full(1-Jie/Jiemax)).^mu-alpha).*Corr_ie); % NOTICE that I put Corr_ei instead of Corr_ie because there is symmetry Cei,ij=Cie,ji in the FP regime
    elseif Jeimean*Jiemean<JbarD^2 && i>=death_time  % After death of neurons the bif. line has changed and there is sometime in the new FP region
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
        [Corr_ie,Corr_ei,~,~,m_e_T,m_i_T,~,~,~]=proj.common.Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf); % find cross-correlations
        if i<death_time && i>death_time-10 && Jeimean*Jiemean>Jeibef*Jiebef % find the neural dynamics ans synapses before death of neurons
            [m_e_bef_per,m_i_bef_per,~,~,timem] = proj.common.Two_populations_full_rate_model_history(m_e_history,m_i_history,Jee,Jei,Jie,Jii,dt,tf);
            Jeibef=mean(Jei(:));
            Jiebef=mean(Jie(:));
        elseif i>=death_time % find the neural dynamics and synapses after death of neurons
            m_i_T=m_i_T(num_dead+1:end,:);
            if i>length(t)-5 && ~flagafter
                [m_e_aft_per,m_i_aft_per,~,~,timem] = proj.common.Two_populations_full_rate_model_history(m_e_history,m_i_history,Jee,Jei,Jie,Jii,dt,tf);
                Jeiaft=mean(nonzeros(Jei(:)));
                Jieaft=mean(nonzeros(Jie(:)));
                flagafter=1;
            end
        end
        [~,~,T,~,~] = proj.common.Two_populations_full_rate_model_history(m_e_history,m_i_history,Jee,Jei,Jie,Jii,dt,tf); % find the time period of the system
        frequency(i)=1/mean(T); % frequency of the system
        stdmebar(i)=std(mean(m_e_T,2)); % SD of me bar
        stdmibar(i)=std(mean(m_i_T,2)); % SD of mi bar
        stdmetil(i)=std(max(m_e_T,[],2)-mean(m_e_T,2)); % SD of me tilde
        stdmitil(i)=std(max(m_i_T,[],2)-mean(m_i_T,2)); % SD of mi tilde
        for n=1:N_e % Compute STDP dynamics
            for k=1:N_i
                J_ei_dot(n,k)=lambda_i*((K_pI-alpha*K_mI)*squeeze(Corr_ie(k,n,:))*dt);
                J_ie_dot(k,n)=lambda_e*(((Th_li_full(1-Jie(k,n)/Jiemax))^mu*K_pE-alpha*K_mE)*squeeze(Corr_ei(n,k,:))*dt);
            end
        end
    end
    noiseI=0.01; % SD of noise of Jei synspses
    noiseE=10*noiseI; % SD of noise of Jie synspses
    Jei=(Jei+dtlearn*(J_ei_dot+noiseI*randn(size(Jei)))).*heaviside(Jei+dtlearn*(J_ei_dot+noiseI*randn(size(Jei)))); % STDP update
    Jie=(Jie+dtlearn*(J_ie_dot+noiseE*randn(size(Jie)))).*heaviside(Jie+dtlearn*(J_ie_dot+noiseE*randn(size(Jie)))); % STDP update
    
end
end