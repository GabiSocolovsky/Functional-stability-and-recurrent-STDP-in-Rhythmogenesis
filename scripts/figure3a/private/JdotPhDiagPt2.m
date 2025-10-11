%% Vector Flow in Jei-Jie phase diagram Pt.2 %%
function [vei,vie] = JdotPhDiagPt2(Jei,Jie_arr,Jee,Jii,JbarD,lambda_i,lambda_e,alpha,KmE,KmI,KpE,KpI,mu,J_ie_max,dt,tf)
% This function computes the vector flow of Jei and Jie (pt.2), belongs to:
% "Functional stability and recurrent STDP in Rhythmogenesis"

%   Description:

%       The function recieves runs on different values of Jie and computes
%       the vector flow - Jeidot and Jiedot

%   Inputs: Jei,Jie_arr,Jee,Jii,JbarD,lambda_i,lambda_e,alpha,KmE,KmI,KpE,KpI,mu,J_ie_max,dt,tf

%       Jei    -    Jei order parameter
%       Jie_arr  -   Jie array on which the loop will run
%       Jee          -   Jee order parameter
%       Jii          -   Jii order parameter
%       JbarD   -   Jbar on bifurcation
%       lambda_i, lambda_e   -  learning rate for synapse Jei and Jie
%       respectively
%       alpha    -    relative strength of depression
%       KpE, KmE, KpI, KmI - K of potentiation (p) or depression (m) for
%       each synapse
%       mu   -   non-linearity
%       J_ie_max  -  JIEmax
%       dt            -   time bin
%       tf            -   final time of simulation

%   Outputs: vei,vie

%       vei, vie -   vector of Jeidot (vei) and Jiedot (vie) for different Jie

%   Dependencies: 

%       Correlations_2D_full_diff  -  Computes the correlations 

%
% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vei=zeros(1,length(Jie_arr)); % vector of Jeidot for different Jie (empty)
vie=zeros(1,length(Jie_arr)); % vector of Jiedot for different Jie (empty)

for j=1:length(Jie_arr)
    Jie=Jie_arr(j);
    f=(1-Jie/J_ie_max)^mu; % synptic dependent function
    if Jei*Jie<JbarD^2 % FP Regime
        me=(1+Jii-Jei)/((1-Jee)*(1+Jii)+Jei*Jie);
        mi=(1-Jee+Jie)/((1-Jee)*(1+Jii)+Jei*Jie);
        vei(j)=lambda_i*(1-alpha)*me*mi;
        vie(j)=lambda_e*(f-alpha)*me*mi;
    else % R Regime
        [Corr_ie,Corr_ei,~,~,~,~,~,~,~]=proj.common.Correlations_2D_full_diff(Jee,Jei,Jie,Jii,dt,tf);
        vei(j)=lambda_i*(KpI-alpha*KmI)*squeeze(Corr_ie)*dt;
        vie(j)=lambda_e*(f*KpE-alpha*KmE)*squeeze(Corr_ei)*dt;
    end
end

end