%% Nullclines Pt.1 %%
function [NullclineJei,NullclineJie,dJei,dJie,Jei_arr,Jie_arr] = NullclinesParallelPt1(JbarD,lambda_i,lambda_e,alpha,tau_pE,tau_mE,tau_pI,tau_mI,H_E,H_I,KpE,KmE,KpI,KmI,mu,J_ie_max,Jee,Jii,dt,tf)
% This function computes the nullclines of Jei and Jie, belonges to:
% "Functional stability and recurrent STDP in Rhythmogenesis"

%   Description:

%       The function parallely (parfor) for each Jei goes to
%       the function NullclinesParallelPt2 which computes for different Jie
%       the nullcline (checks if there is a nullcline and finds it)

%   Inputs: 

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

%       NullclineJei    -   Nullcline of Jei
% ,     NullclineJie    -   Nullcline of Jie
%       dJei            -   resolution of Jei
%       dJie            -   resolution of Jie
%       Jei_arr         -   array of Jei
%       Jie_arr         -   array of Jie

%   Dependencies:

%   NullclinesParallelPt2.m          -   computes for different Jie the nullcline (checks if there is a nullcline and finds it)
%
% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grid of Jei and Jie %%
Jei_init=0.01;
dJei=0.005;
Jei_fin=1+Jii-dJei;
Jie_init=0.05;
dJie=0.1;
Jie_fin=J_ie_max-dJie;

Jei_arr=Jei_init:dJei:Jei_fin;
Jie_arr=Jie_init:dJie:Jie_fin;
%% compute nullcline %%
NullclineJei=zeros(length(Jei_arr),2);
NullclineJie=zeros(length(Jie_arr),2);

parfor i=1:length(Jei_arr) % computes parallely the nullcline for each Jei
    Jei=Jei_arr(i);
    [NullclineJei(i,:),NullclineJie(i,:)]=NullclinesParallelPt2(JbarD,lambda_i,lambda_e,alpha,tau_pE,tau_mE,tau_pI,tau_mI,H_E,H_I,KpE,KmE,KpI,KmI,mu,J_ie_max,dt,tf,Jee,Jii,Jei,Jie_arr,Jie_init,dJie,Jie_fin)
end


end
