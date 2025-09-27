function [Vnvol,Vevol,Vavol,Vdvol,stdmebar,stdmibar,stdmetil,stdmitil,Jeimeandyn,Jiemeandyn,dJ_arr_t] = FullSynDynMultiTrails(m_e_history,m_i_history,Jei,Jie,Jee,Jii,t,Jiemax,mu,alpha,JbarD,Jeemean,Jiimean,vn,ve,vu,vaei,vaie,vdei,vdie,mulvn,mulve,mulva,mulvd,dtlearn,lambda_e,lambda_i,K_Ibar,K_Itilmphi,K_Itilphi,K_Eptilmpsi,K_Etilphi,K_Etilmphi,K_pI,K_pE,K_mI,K_mE,N_e,N_i,dt,tf)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Vnvol=nan(1,length(t));
Vevol=nan(1,length(t));
Vavol=nan(1,length(t));
Vdvol=nan(1,length(t));

stdmebar=nan(1,length(t));
stdmibar=nan(1,length(t));
stdmetil=nan(1,length(t));
stdmitil=nan(1,length(t));

Jeimeandyn=nan(1,length(t));
Jiemeandyn=nan(1,length(t));
dJ_arr_t=nan(length(t),N_e*N_i*2);

for i=1:length(t)
    dJ_arr=[reshape(Jei.',1,[])-mean(Jei(:)) reshape(Jie.',1,[])-mean(Jie(:))];
    dJ_arr_t(i,:)=dJ_arr;
    %%% construcing the eigenvectors on each time step %%%
    Jeimean=mean(Jei(:));
    Jiemean=mean(Jie(:));
    
    Jeimeandyn(i)=Jeimean;
    Jiemeandyn(i)=Jiemean;    
    
    meb=(1+Jiimean-Jeimean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean); % DC part of me (needed in both cases of the condition)
    mib=(1-Jeemean+Jiemean)/((1-Jeemean)*(1+Jiimean)+Jeimean*Jiemean); % DC part of mi (needed in both cases of the condition)
    
    f=(1-Jiemean/Jiemax)^mu;
    ftag=-mu/Jiemax*(1-Jiemean/Jiemax)^(mu-1);
    K_Ebar=f-alpha;
    
    if Jeimean*Jiemean<JbarD^2 % FP Regime
        a=-lambda_i*K_Ibar*mib^2;
        b=lambda_i*K_Ibar*meb^2;
        c=-lambda_e*K_Ebar*mib^2;
        d=lambda_e*K_Ebar*meb^2;
        q=lambda_e*ftag*meb*mib;
    else % R regime
        gabs=(JbarD^2-Jeemean*Jiimean)^0.5;
        mit=(JbarD^2/Jeimean*(1+Jiimean)-JbarD^2)/(gabs*JbarD*((1-Jeemean)*(1+Jiimean)+JbarD^2));
        met=JbarD/Jiemean*mit;
        a=-lambda_i*(K_Ibar*mib^2+K_Itilmphi*mit^2/(2*gabs));
        b=lambda_i*(K_Ibar*meb^2+K_Itilphi*met^2/(2*gabs));
        c=-lambda_e*(K_Ebar*mib^2+K_Etilphi*mit^2/(2*gabs));
        d=lambda_e*(K_Ebar*meb^2+K_Etilmphi*met^2/(2*gabs));
        q=lambda_e*ftag*(meb*mib+met*mit/2*K_Eptilmpsi); % self deperssion coefficient    end
    end
    
    %vn=[ 0 0 0 0 -1 1 0 0 0 0 0 0 ; 0 0 -1 1 0 0 0 0 0 0 0 0 ; -1 1 0 0 0 0 0 0 0 0 0 0]/2^0.5;
    %ve=[0 0 0 0 0 0 0 0 0 -1 0 1 ; 0 0 0 0 0 0 0 0 0 -1 1 0 ; 0 0 0 0 0 0 -1 0 1 0 0 0 ; 0 0 0 0 0 0 -1 1 0 0 0 0]/2^0.5;
    %va1=[(q-a) (q-a) (-q+a) (-q+a) 0 0 -c c 0 -c c 0]/norm([(q-a) (q-a) (-q+a) (-q+a) 0 0 -c c 0 -c c 0]);
    %va2=[0 0 (q-a) (q-a) (-q+a) (-q+a) 0 -c c 0 -c c]/norm([0 0 (q-a) (q-a) (-q+a) (-q+a) 0 -c c 0 -c c]);
    %va=[va1 ; va2];
    %vd=[-b b -b b -b b -(d+q) -(d+q) -(d+q) d+q d+q d+q]/norm([-b b -b b -b b -(d+q) -(d+q) -(d+q) d+q d+q d+q]);
    %vu=[1 1 1 1 1 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1]/12^0.5;
    
    va=[(a-q)*vaei c*vaie]'/(((a-q)^2+c^2)*2*N_i)^0.5;
    vd=[b*vdei (d+q)*vdie]'/((b^2+(d+q)^2)*2*N_e)^0.5;
    
    V=[vn  ve  va  vd  vu];
    coef=V\dJ_arr';
    abcoef=abs(coef);
    loc=find(heaviside(abcoef-10^-6));
    
    
    %%%% The volume of each subspace %%%%
    Vnvol(i)=(coef(1:mulvn)'*vn'*vn*coef(1:mulvn))^0.5;
    Vevol(i)=(coef((mulvn+1):(mulvn+mulve))'*ve'*ve*coef((mulvn+1):(mulvn+mulve)))^0.5;
    Vavol(i)=(coef((mulvn+mulve+1):(mulvn+mulve+mulva))'*va'*va*coef((mulvn+mulve+1):(mulvn+mulve+mulva)))^0.5;
    Vdvol(i)=(coef((mulvn+mulve+mulva+1):(mulvn+mulve+mulva+mulvd))'*vd'*vd*coef((mulvn+mulve+mulva+1):(mulvn+mulve+mulva+mulvd)))^0.5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if Jeimean*Jiemean<JbarD^2 % in the FP regime we solve exactly
        J=[[Jee/N_e -Jei/N_i];[Jie/N_e -Jii/N_i]]; % connectivity matrix
        diagonal=eye(N_e+N_i);
        mfa=(diagonal-J)\ones(N_e+N_i,1); % The exact fixed point solution
        stdmebar(i)=std(mfa(1:N_e));
        stdmibar(i)=std(mfa(N_e+1:N_e+N_i));
        stdmetil(i)=nan;
        stdmitil(i)=nan;
        Corr_ei=mfa(1:N_e)*mfa((N_e+1):(N_e+N_i))';
        Corr_ie=Corr_ei';
        J_ei_dot=lambda_i*((1-alpha)*Corr_ei); % NOTICE that I put Corr_ei instead of Corr_ie because there is symmetry Cei,ij=Cie,ji in the FP regime
        J_ie_dot=lambda_e*(((Th_li_full(1-Jie/Jiemax)).^mu-alpha).*Corr_ie); % NOTICE that I put Corr_ei instead of Corr_ie because there is symmetry Cei,ij=Cie,ji in the FP regime
    else % The R regime is solved with simulations
        [Corr_ie,Corr_ei,~,~,m_e_T,m_i_T,~,~,~]=proj.common.Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf);
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