function f = FrequencyContourJbarJee(JbarD,Jbar_arrlast,Jee,Jii,dt,tf)
f=nan(1,length(Jbar_arrlast));
for j=1:length(Jbar_arrlast)
    Jbar=Jbar_arrlast(j);
    if isnan(JbarD)
        if Jee>2*Jbar-Jii
            continue
        end
        [~,~,Tperiod,~,~]=proj.common.Two_populations_full_rate_model_history(0.2,0.5,Jee,0.5,Jbar^2/0.5,Jii,dt,tf);
        f(j)=1/Tperiod;
    elseif Jbar<JbarD
        continue
    else
        [~,~,Tperiod,~,~]=proj.common.Two_populations_full_rate_model_history(0.2,0.5,Jee,0.5,Jbar^2/0.5,Jii,dt,tf);
        f(j)=1/Tperiod;
    end
end