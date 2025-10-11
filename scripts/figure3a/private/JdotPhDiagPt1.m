%% Vector Flow in Jei-Jie phase diagram Pt.1 %%
function [Jeidot,Jiedot,Jei_arr,Jie_arr] = JdotPhDiagPt1(JbarD,lambda_i,lambda_e,alpha,tau_mE,tau_mI,tau_pE,tau_pI,H_E,H_I,mu,J_ie_max,Jee,Jii,dt,tf)
% This function computes the vector flow of Jei and Jie (pt.1), belongs to:
% "Functional stability and recurrent STDP in Rhythmogenesis"

%   Description:

%       The function parallely (parfor) for each Jei goes to
%       the function NullclinesParallelPt2 which computes for different Jie
%       the nullcline (checks if there is a nullcline and finds it)

%   Inputs: 

%       JbarD   -   Jbar on bifurcation
%       lambda_i, lambda_e   -  learning rate for synapse Jei and Jie
%       respectively
%       alpha    -    relative strength of depression
%       tau_pE, tau_mE, tau_pI, tau_mI - typical times synapse Jei (I) or
%       H_E,H_I - Hebbianity of Jei (I) and Jie (E)
%       KpE, KmE, KpI, KmI - K of potentiation (p) or depression (m) for
%       each synapse
%       mu   -   non-linearity
%       J_ie_max  -  JIEmax
%       Jee          -   Jee order parameter
%       Jii          -   Jii order parameter
%       dt            -   time bin
%       tf            -   final time of simulation

%   Outputs: 

%       Jeidot          - array of time derivative of Jei
%       Jiedot          - array of time derivative of Jie
%       Jei_arr         -   array of Jei
%       Jie_arr         -   array of Jie

%   Dependencies: 

%      Correlations_2D_full_diff  -  Comutes the correlations - we use it in
%   order to have the temporal kernels defined for times -Delta_extended to
%   Delta_extended

%   JdotPhDiagPt2.m   -   computes for different Jie the time detivative of Jei and Jie
%
% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the grid and create temporal kernels %%
Jei_init=0.05;
dJei=0.1;
Jei_final=1+Jii;
Jei_arr=Jei_init:dJei:(Jei_final);

Jie_init=0.5;
dJie=1;
Jie_final=20;
Jie_arr=Jie_init:dJie:(Jie_final);


Jeidot=zeros(length(Jei_arr),length(Jie_arr));
Jiedot=zeros(length(Jei_arr),length(Jie_arr));

[~,~,~,~,~,~,~,~,Delta_extended]=proj.common.Correlations_2D_full_diff(Jee,0.5,10,Jii,dt,tf); % calculate Delta_extended - the time on which the kernels and correlations are integrated together
KpE=1/tau_pE*exp(-Delta_extended*H_E/tau_pE).*heaviside(H_E*Delta_extended); % K plus of synapses Jie
KmE=1/tau_mE*exp(Delta_extended*H_E/tau_mE).*heaviside(-H_E*Delta_extended); % K minus of synapses Jie
KpI=1/tau_pI*exp(-Delta_extended*H_I/tau_pI).*heaviside(H_I*Delta_extended); % K plus of synapses Jei
KmI=1/tau_mI*exp(Delta_extended*H_I/tau_mI).*heaviside(-H_I*Delta_extended); % K minus of synapses Jei
%% Compute parallely Jeidot and Jiedot for different Jei %% 
parfor i=1:length(Jei_arr) 
    Jei=Jei_arr(i);
    [vei,vie]=JdotPhDiagPt2(Jei,Jie_arr,Jee,Jii,JbarD,lambda_i,lambda_e,alpha,KmE,KmI,KpE,KpI,mu,J_ie_max,dt,tf)
    Jeidot(i,:)=vei;
    Jiedot(i,:)=vie;
end