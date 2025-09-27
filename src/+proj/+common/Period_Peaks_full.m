function [mean_period,std_period] = Period_Peaks_full(m,dt)
% Finding the Period of a function

% This function is more suitable for a signal that is already known to be
% periodic.

std_period=zeros(1,size(m,1));
mean_period=zeros(1,size(m,1));
t_int=round(7*size(m,2)/10)+1;
for j=1:size(m,1)
[vals,locs] = findpeaks(m(j,t_int:end)); % finds the indices of the peaks    
    if abs(diff(vals))<10^-13
    std_period(j)=nan;
    mean_period(j)=nan;
    continue
    end
std_period(j)=std(dt*diff(locs));
mean_period(j)=dt*mean(diff(locs));
end

end

