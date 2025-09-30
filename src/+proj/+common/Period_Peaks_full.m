%% Finding the time Period of a variable %%
function [mean_period,std_period] = Period_Peaks_full(m,dt) 
%   This function computes the time period of a periodic function

%   Description:

%       Takes the variable m and finds its peaks with "findpeaks" after a
%       steady state time. By the time difference of this peaks it computes
%       the average time period.
    
%   Notes:  

%    -  This function is more suitable for a signal that is already known to be
%       periodic and that has one local maximum per period.

%   Inputs:

%       m   -   The firing rates
%       dt  -   Time bin 

%   Outputs:

%       mean_period           -   the average time period of m
%       std_period            -   the standard deviation across different
%       cycles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Output variables and others %%%%%%%%%
std_period=zeros(1,size(m,1)); % the standard deviation across different (empty array)
mean_period=zeros(1,size(m,1)); % the average time period of m (empty array)
t_int=round(7*size(m,2)/10)+1; % an arbitrary time in which steady state is assumed
for j=1:size(m,1)
[vals,locs] = findpeaks(m(j,t_int:end)); % finds the indices of the peaks    
    if abs(diff(vals))<10^-13 % very small amplitude differences - which indicate convergence to a constant value
    std_period(j)=nan; % no time period when activity is constant
    mean_period(j)=nan; % no time period when activity is constant
    continue
    end
std_period(j)=std(dt*diff(locs)); % the standard deviation across different
mean_period(j)=dt*mean(diff(locs)); % the average time period of m
end

end

