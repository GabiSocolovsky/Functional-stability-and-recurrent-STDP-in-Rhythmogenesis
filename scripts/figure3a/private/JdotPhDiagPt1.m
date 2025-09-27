function [Jeidot,Jiedot,Jei_arr,Jie_arr] = JdotPhDiagPt1(JbarD,lambda_i,lambda_e,alpha,tau_mE,tau_mI,tau_pE,tau_pI,H_E,H_I,mu,J_ie_max,Jee,Jii,dt,tf)
% This function computes Jeidot and Jiedot (only order parameters - 2D)
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

[~,~,~,~,~,~,~,~,Delta_extended]=proj.common.Correlations_2D_full_diff(Jee,0.5,10,Jii,dt,tf);
KpE=1/tau_pE*exp(-Delta_extended*H_E/tau_pE).*heaviside(H_E*Delta_extended);
KmE=1/tau_mE*exp(Delta_extended*H_E/tau_mE).*heaviside(-H_E*Delta_extended);
KpI=1/tau_pI*exp(-Delta_extended*H_I/tau_pI).*heaviside(H_I*Delta_extended);
KmI=1/tau_mI*exp(Delta_extended*H_I/tau_mI).*heaviside(-H_I*Delta_extended);

parfor i=1:length(Jei_arr)
    Jei=Jei_arr(i);
    [vei,vie]=JdotPhDiagPt2(Jei,Jie_arr,Jee,Jii,JbarD,lambda_i,lambda_e,alpha,KmE,KmI,KpE,KpI,mu,J_ie_max,dt,tf)
    Jeidot(i,:)=vei;
    Jiedot(i,:)=vie;
end