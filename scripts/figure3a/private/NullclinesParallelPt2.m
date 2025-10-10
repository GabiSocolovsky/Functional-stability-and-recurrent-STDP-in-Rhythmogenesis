%% Nullclines Pt. 2 %%
function [NullclineJei,NullclineJie] = NullclinesParallelPt2(JbarD,lambda_i,lambda_e,alpha,tau_pE,tau_mE,tau_pI,tau_mI,H_E,H_I,KpE,KmE,KpI,KmI,mu,J_ie_max,dt,tf,Jee,Jii,Jei,Jie_arr,Jie_init,dJie,Jie_fin)
% This function computes the nullclines of Jei and Jie, belonges to:
% "Functional stability and recurrent STDP in Rhythmogenesis"

%   Description:

%       The function computes for different Jie
%       the nullcline (checks if there is a nullcline and finds it). The
%       way of knowing if there is a nullcline is if Jeidot or Jiedot
%       changes its sign. If we find these two values where the change of sign happens
%       than I increase the resoultion in dJie and find between these two
%       values a more accurate value

%   Inputs: JbarD,lambda_i,lambda_e,alpha,tau_pE,tau_mE,tau_pI,tau_mI,H_E,H_I,KpE,KmE,KpI,KmI,mu,J_ie_max,Jee,Jii,dt,tf

%       JbarD   -   Jbar on bifurcation
%       lambda_i, lambda_e   -  learning rate for synapse Jei and Jie
%       respectively
%       tau_pE, tau_mE, tau_pI, tau_mI - typical times synapse Jei (I) or
%       H_E,H_I - Hebbianity of Jei (I) and Jie (E)
%       KpE, KmE, KpI, KmI - K of potentiation (p) or depression (m) for
%       each synapse
%       mu   -   non-linearity
%       J_ie_max  -  JIEmax
%       Jie (E) for potentiation (p) or depression (m)
%       Jee          -   Jee order parameter
%       Jii          -   Jii order parameter
%       dt            -   time bin
%       tf            -   final time of simulation

%   Outputs:

%       m_e                      -   Dynamics of Excitatory neurons firing rates
%       m_i                      -   Dynamics of Inhibitory neurons firing rates
%       T_mean_m_e,T_mean_m_i    -   Time period of m_e and m_i (sometimes
%       different when in only the inhibitory neurons are rhythmic)
%       time                     -   Vector of time

%   Dependencies:

%   Correlations_2D_full_diff.m  - computes correlations
%
% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qei=0; % the cue that dt is smaller (q=0 - no , q=1 - yes) for null of Jei
qie=0; % the cue that dt is smaller (q=0 - no , q=1 - yes) for null of Jie
doneei=0; % the cue for a point in nullcline of Jei that is already made
doneie=0; % the cue for a point in nullcline of Jie that is already made

dtfirst=dt;
dt=dt/10;
[~,~,~,~,~,~,~,~,Delta_extended]=proj.common.Correlations_2D_full_diff(Jee,0.5,10,Jii,dt,tf);
K_pE_smallerdt=1/tau_pE*exp(-Delta_extended*H_E/tau_pE).*heaviside(H_E*Delta_extended);
K_mE_smallerdt=1/tau_mE*exp(Delta_extended*H_E/tau_mE).*heaviside(-H_E*Delta_extended);
K_pI_smallerdt=1/tau_pI*exp(-Delta_extended*H_I/tau_pI).*heaviside(H_I*Delta_extended);
K_mI_smallerdt=1/tau_mI*exp(Delta_extended*H_I/tau_mI).*heaviside(-H_I*Delta_extended);
dt=dtfirst;

Jiedot=zeros(1,length(Jie_init:dJie:Jie_fin));
Jeidot=zeros(1,length(Jie_init:dJie:Jie_fin));
Jieback=1;
NullclineJei=[nan nan];
NullclineJie=[nan nan];

j=1;

while j<=length(Jie_arr)
    Jie=Jie_arr(j);
    f=(1-Jie/J_ie_max)^mu;
    if Jei*Jie<JbarD^2 % FP Regime
        me=(1+mean(Jii)-Jei)/((1-mean(Jee))*(1+mean(Jii))+Jei*Jie);
        mi=(1-mean(Jee)+Jie)/((1-mean(Jee))*(1+mean(Jii))+Jei*Jie);
        Jeidot(j)=lambda_i*(1-alpha)*me*mi;
        Jiedot(j)=lambda_e*(f-alpha)*me*mi;
    else % R Regime
        [Corr_ie,Corr_ei,~,~,~,~,~,~,~]=proj.common.Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf);
        if dt==0.01
            Jeidot(j)=lambda_i*(KpI-alpha*KmI)*squeeze(Corr_ie)*dt;
            Jiedot(j)=lambda_e*(f*KpE-alpha*KmE)*squeeze(Corr_ei)*dt;
        else
            Jeidot(j)=lambda_i*(K_pI_smallerdt-alpha*K_mI_smallerdt)*squeeze(Corr_ie)*dt;
            Jiedot(j)=lambda_e*(f*K_pE_smallerdt-alpha*K_mE_smallerdt)*squeeze(Corr_ei)*dt;
        end
    end
    if j==1
        j=2;
        continue
    end
    if sign(Jeidot(j))~=sign(Jeidot(j-1)) && sign(Jiedot(j))~=sign(Jiedot(j-1))
        if ~qei && ~qie
            j=j-round(Jieback/dJie)-1; % come back round(Jieback/dJie) steps backwards (Jie=Jie-1)
            dt=dt/10;
            qei=1; % the cue switches from zero to one and vice versa
            qie=1; % the cue switches from zero to one and vice versa
        else
            if abs(Jie-JbarD^2/Jei)<0.01 % are we on the bifurcation line
                NullclineJei=[Jei JbarD^2/Jei];
                NullclineJie=[Jei JbarD^2/Jei];
            else % or on the nullcline itself
                NullclineJei=[Jei Jie-dJie-Jeidot(j-1)/((Jeidot(j)-Jeidot(j-1))/dJie)];
                NullclineJie=[Jei Jie-dJie-Jiedot(j-1)/((Jiedot(j)-Jiedot(j-1))/dJie)];
            end
            %dt=dt*10;
            %indei=indei+1;
            %indie=indie+1;
            break
        end
    elseif sign(Jeidot(j))~=sign(Jeidot(j-1)) % if we find a nullcine, then we go a Jie=Jie-1, and increase accuracy
        if ~qei % when q==0 we lower the accuracy in dt
            if qie==1 && ~doneie
                if abs(Jie-JbarD^2/Jei)<0.01 % are we on the bifurcation line
                    NullclineJei=[Jei JbarD^2/Jei];
                else % or on the nullcline itself
                    NullclineJei=[Jei Jie-dJie-Jeidot(j-1)/((Jeidot(j)-Jeidot(j-1))/dJie)];
                end
                doneei=1;
                j=j+1;
                continue
            end
            j=j-round(Jieback/dJie)-1; % come back round(Jieback/dJie) steps backwards (Jie=Jie-1)
            dt=dt/10;
            qei=1; % the cue switches from zero to one and vice versa
        elseif qei==1 && ~doneei % otherwise we take the point of the nullcline
            if abs(Jie-JbarD^2/Jei)<0.01 % are we on the bifurcation line
                NullclineJei=[Jei JbarD^2/Jei];
            else % or on the nullcline itself
                NullclineJei=[Jei Jie-dJie-Jeidot(j-1)/((Jeidot(j)-Jeidot(j-1))/dJie)];
            end
            dt=dt*10;
            doneei=1;
            if qie==1
                if doneie
                    break
                end
                dt=dt/10;
            end
        end
    elseif sign(Jiedot(j))~=sign(Jiedot(j-1))
        if ~qie % when q==0 we lower the accuracy in dt
            if qei==1 && ~doneei
                if abs(Jie-JbarD^2/Jei)<0.01 % are we on the bifurcation line
                    NullclineJie=[Jei JbarD^2/Jei];
                else % or on the nullcline itself
                    NullclineJie=[Jei Jie-dJie-Jiedot(j-1)/((Jiedot(j)-Jiedot(j-1))/dJie)];
                end
                dt=dt*10;
                %indie=indie+1;
                doneie=1;
                j=j+1;
                continue
            end
            j=j-round(Jieback/dJie)-1; % come back round(Jieback/dJie) steps backwards (Jie=Jie-1)
            dt=dt/10;
            qie=1; % the cue switches from zero to one and vice versa
        elseif qie==1 && ~doneie % otherwise we take the point of the nullcline
            if abs(Jie-JbarD^2/Jei)<0.01 % are we on the bifurcation line
                NullclineJie=[Jei JbarD^2/Jei];
            else % or on the nullcline itself
                NullclineJie=[Jei Jie-dJie-Jiedot(j-1)/((Jiedot(j)-Jiedot(j-1))/dJie)];
            end
            dt=dt*10;
            %indie=indie+1;
            doneie=1;
            if qei==1
                if doneei
                    break
                end
                dt=dt/10;
            end
        end
    end
    j=j+1;
end
end

