%% Frequency contour in Jbar-Jee phase diagram %%
function f = FrequencyContourJbarJee(JbarD,Jbar_arrlast,Jee,Jii,dt,tf)
% This function computes the freuqnecy contour in the Jbar-Jee phase
% diagram from "Functional stability and recurrent STDP in Rhythmogenesis"

%   Description:

%  Computes the frequency for a grid of Jbar-Jee phase diagram.
%  It returns nans outside of the R region.

%   Inputs:

%       JbarD         -   Jbar on the bif. line
%       Jbar_arrlast  -   Jbar_arr (for which the frequency is calculated)
%       Jee           -   E-to-E order parameter synapse
%       Jii           -   I-to-I order parameter synapse
%       dt            -   time bin
%       tf            -   final time of simulation

%   Outputs:

%       f             -   grid of freuqnecies for Jbar-Jee

%   Dependencies:

%   Two_populations_full_rate_model_history.  - computes neural dynamics,
%   required for computing the time period

% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-09-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=nan(1,length(Jbar_arrlast)); % frequency grid (empty)
for j=1:length(Jbar_arrlast)
    Jbar=Jbar_arrlast(j);
    if isnan(JbarD)
        if Jee>2*Jbar-Jii % outside the R region
            continue
        end
        [~,~,Tperiod,~,~]=proj.common.Two_populations_full_rate_model_history(0.2,0.5,Jee,0.5,Jbar^2/0.5,Jii,dt,tf); % compute the time period
        f(j)=1/Tperiod; % freuqnecy
    elseif Jbar<JbarD
        continue
    else
        [~,~,Tperiod,~,~]=proj.common.Two_populations_full_rate_model_history(0.2,0.5,Jee,0.5,Jbar^2/0.5,Jii,dt,tf); % compute the time period
        f(j)=1/Tperiod; % freuqnecy
    end
end