function [vei,vie] = JdotPhDiagPt2(Jei,Jie_arr,Jee,Jii,JbarD,lambda_i,lambda_e,alpha,KmE,KmI,KpE,KpI,mu,J_ie_max,dt,tf)
% This function computes Jeidot and Jiedot (only order parameters - 2D)
vei=zeros(1,length(Jie_arr));
vie=zeros(1,length(Jie_arr));

for j=1:length(Jie_arr)
    Jie=Jie_arr(j);
    f=(1-Jie/J_ie_max)^mu;
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