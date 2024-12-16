


% load estimated parameters
if isfile([EP_options.estimation_file_name '_mean.mat']) && strcmp(EP_options.load_theta,'mean')
	load([EP_options.estimation_file_name '_mean.mat'])
elseif isfile([EP_options.estimation_file_name '_mean.mat']) && strcmp(EP_options.load_theta,'MHmode')
	load([EP_options.estimation_file_name '_mean.mat'])
    theta=mode_metropolis;
elseif isfile([EP_options.estimation_file_name '_mode.mat'])  && strcmp(EP_options.load_theta,'mode')
	load([EP_options.estimation_file_name '_mode.mat'])
else
    load([EP_options.estimation_file_name '_mode.mat'])
	load([EP_options.estimation_file_name '_mle_temp.mat'])
    filtered_errors = e_obj.spfm_exo_simul;
    filtered_errors((1:size(best_exo,1))+1,:) = best_exo;
    theta=params;
    filtered_data=best_filtered;
end


% update order in case parameter order changed
for i1=(1+length(e_obj.theta.is_sd)):length(e_obj.theta.id)
    if ~isempty(strmatch(deblank(e_obj.theta.names(i1,:)),M_.param_names,'exact'))
        set_param_value(deblank(e_obj.theta.names(i1,:)), theta(i1));
        eval([ deblank(e_obj.theta.names(i1,:)) ' = ' num2str(theta(i1)) ]);
     end
end 

%% load smoothed data
load([EP_options.estimation_file_name '_filtered_data'],'filtered_errors','filtered_data','filtered_errors_ts','variable_names');
Tdata           = size(filtered_errors,1);
exo_filtered    = zeros(Tdata,M_.exo_nbr);
endo_filtered   = repmat(oo_.steady_state,[1 size(filtered_data,2)]);

% shocks
M_.Sigma_e = zeros(M_.exo_nbr);
for i1=1:size(filtered_errors_ts.name,1)
    % update covariance matrix
	idx = strmatch(deblank(filtered_errors_ts.name{i1,:}),M_.exo_names,'exact');
    idtheta = strmatch(deblank(filtered_errors_ts.name{i1,:}),e_obj.theta.names,'exact');
    if ~isempty(idtheta)
        M_.Sigma_e(idx,idx)=theta(idtheta).^2;
    end
    % update filtered shock matrix
	idx = strmatch(filtered_errors_ts.name{i1},M_.exo_names,'exact');
    if ~isempty(idx)
        exo_filtered(:,idx) = filtered_errors(:,i1);
    end
end;

% endogenous variable_names
endo_filtered = repmat(oo_.steady_state,[1 size( filtered_data,2)]);
for i1 = 1:size(variable_names,1)
	idx = strmatch(variable_names{i1},M_.endo_names,'exact');
    
    if ~isempty(idx)
        endo_filtered(idx,:) =  filtered_data(i1,:);
    end
end

% set them in adequate ss form
steady;

