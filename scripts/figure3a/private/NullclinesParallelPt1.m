function [NullclineJei,NullclineJie,dJei,dJie,Jei_arr,Jie_arr] = NullclinesParallelPt1(JbarD,lambda_i,lambda_e,alpha,tau_pE,tau_mE,tau_pI,tau_mI,H_E,H_I,KpE,KmE,KpI,KmI,mu,J_ie_max,Jee,Jii,dt,tf)
% Finding the Nullclines of the phase diagram %
% These nullclines are approximations in the case of very small DeltaJ

Jei_init=0.01;
dJei=0.005;
Jei_fin=1+Jii-dJei;
Jie_init=0.05;
dJie=0.1;
Jie_fin=J_ie_max-dJie;

Jei_arr=Jei_init:dJei:Jei_fin;
Jie_arr=Jie_init:dJie:Jie_fin;




%dt=dtfirst;

NullclineJei=zeros(length(Jei_arr),2);
NullclineJie=zeros(length(Jie_arr),2);

parfor i=1:length(Jei_arr)
    Jei=Jei_arr(i);
    [NullclineJei(i,:),NullclineJie(i,:)]=NullclinesParallelPt2(JbarD,lambda_i,lambda_e,alpha,tau_pE,tau_mE,tau_pI,tau_mI,H_E,H_I,KpE,KmE,KpI,KmI,mu,J_ie_max,dt,tf,Jee,Jii,Jei,Jie_arr,Jie_init,dJie,Jie_fin)
    %dt=dtfirst;
    %qei=0;
    %qie=0;
    %doneei=0;
    %doneie=0;
    %j=1;
end


end
