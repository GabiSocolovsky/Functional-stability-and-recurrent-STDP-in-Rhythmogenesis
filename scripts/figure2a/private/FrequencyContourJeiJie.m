function f = FrequencyContourJeiJie(JbarD, Jee, Jei, Jie_arr, Jii, dt, tf)
% FrequencyContourJeiJie Compute oscillation frequencies across Jie values.
%
%   f = FREQUENCYCONTOURJEIJIE(JbarD, Jee, Jei, Jie_arr, Jii, dt, tf)
%   computes the oscillation frequency of a two-population network in the 
%   rhythmic region (R) for different values of the synaptic weight parameter 
%   Jie. 
%
%   If the condition Jei * Jie < JbarD^2 is not satisfied, the network is 
%   outside the R region, and the function returns NaN for that value. 
%
%   This function is intended to be used inside a PARFOR loop over values of Jei.
%
%   INPUTS:
%       JbarD   - Critical Jbar value on the bifurcation line
%       
%       Order parameters:
%           Jee     - Excitatory-to-excitatory synaptic weight 
%           Jei     - Inhibitory-to-excitatory synaptic weight 
%           Jie_arr - Array of excitatory-to-inhibitory synaptic weights 
%           Jii     - Inhibitory-to-inhibitory synaptic weight 
%       
%       dt      - Time step size
%       tf      - Simulation duration
%
%   OUTPUT:
%       f       - Array of oscillation frequencies corresponding to Jie_arr.
%                 NaN indicates non-rhythmic dynamics.
%
%   DEPENDENCIES:
%       Two_populations_full_rate_model_history.m
%
%

    % Initialize frequency vector
    f = nan(1, length(Jie_arr));

    % Loop over Jie values
    for j = 1:length(Jie_arr)
        Jie = Jie_arr(j);

        % Check if network is in rhythmic region
        if Jei * Jie < JbarD^2
            continue
        end

        % Run simulation and extract oscillation period
        [~, ~, Tperiod, ~, ~] = proj.common.Two_populations_full_rate_model_history( ...
            0.2, 0.5, Jee, Jei, Jie, Jii, dt, tf);

        % Store frequency
        f(j) = 1 / Tperiod;
    end
end
