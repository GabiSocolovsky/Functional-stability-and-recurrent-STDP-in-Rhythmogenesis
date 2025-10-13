%% Save Parallely Data %%
function parsave(fname,Vnvol,Vevol,Vavol,Vdvol,stdmebar,stdmibar,stdmetil,stdmitil,Jeimeandyn,Jiemeandyn,dJ_arr_t)
% This function saves the data relevant to figures 6a-6d in the paper:
% "Functional stability and recurrent STDP in Rhythmogenesis"

%   Description:

%       1. Puts the data relevant for figure 6a-6d in a strcture variable -
%       DynVar
%       2. Saves this variable in a another variable - fname# (fname is
%       defined outside the function, and # is the number of trial

%   Inputs:

%       fname - filename
%       Vnvol,Vevol,Vavol,Vdvol - Volume of families I, II, III, IV respectively (not shown in paper
%       stdmebar,stdmibar - SD of me and mi (bar) across population
%       stdmetil,stdmitil - SD of me and mi (tilde) across population
%       Jeimeandyn,Jiemeandyn - mean synaptic weights
%       dJ_arr_t - synaptic fluctuations in a vector first dJei and then dJie

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DynVar=struct('Vnvol',Vnvol,'Vevol',Vevol,'Vavol',Vavol,'Vdvol',Vdvol,'stdmebar',stdmebar,'stdmibar',stdmibar,'stdmetil',stdmetil,'stdmitil',stdmitil,'Jeimeandyn',Jeimeandyn,'Jiemeandyn',Jiemeandyn,'dJ_arr_t',dJ_arr_t);
save(fullfile('data', fname), 'DynVar'); % save DynVar in fname#
end


