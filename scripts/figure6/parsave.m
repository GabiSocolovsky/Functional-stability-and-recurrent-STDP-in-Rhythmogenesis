function parsave(fname,Vnvol,Vevol,Vavol,Vdvol,stdmebar,stdmibar,stdmetil,stdmitil,Jeimeandyn,Jiemeandyn,dJ_arr_t)
DynVar=struct('Vnvol',Vnvol,'Vevol',Vevol,'Vavol',Vavol,'Vdvol',Vdvol,'stdmebar',stdmebar,'stdmibar',stdmibar,'stdmetil',stdmetil,'stdmitil',stdmitil,'Jeimeandyn',Jeimeandyn,'Jiemeandyn',Jiemeandyn,'dJ_arr_t',dJ_arr_t);
save(fullfile('data', fname), 'DynVar');
end


