function ds = dynamic_set_auxiliary_series(ds, params)
%
% Computes auxiliary variables of the dynamic model
%
ds.AUX_ENDO_LAG_4_1=ds.x(-1);
ds.AUX_ENDO_LAG_30_1=ds.z(-1);
end
