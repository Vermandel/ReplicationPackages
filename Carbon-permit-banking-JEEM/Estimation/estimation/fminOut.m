function stop = fminOut(x,optimvalues,state,e_obj)
    
    global M_;

    stop = false;
	updte= 0;
	
	if size(x,2) > 1
		params=x';
	else
		params=x;
	end
	
	%if  optimvalues.iteration <20 ~isempty(getCurrentTask()) % if parallel
        % if parallel, initialize the persistent variables
	%	lnprior = priordens(theta_mode(e_obj.thet_ids),bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4,1);
	%	% one at the end starts the initialization
    %end		
    
	if optimvalues.iteration == 0
		% just starting iteration
		posterior = optimvalues.fval;
		% initialize
		datavec = posterior;
		fstory = [ posterior ];
		xstory = [ x ];
		% save
		save([M_.fname '_mle_estimates_temp'], 'posterior', 'datavec', 'fstory', 'xstory', 'params');
	else
		% get new posterior
		load([M_.fname '_mle_estimates_temp.mat'],'posterior','datavec','fstory','xstory');
		newposterior = optimvalues.fval;

		%	% if improvement
		if newposterior < posterior
		
			fprintf('\n');
			 
			if size(xstory,1) == 1
                xstory = xstory';
            end
			% adopt posterior
			posterior = newposterior;
			% append in history matrices for hessian calculation
			datavec = [datavec; optimvalues.fval];
			fstory  = [ fstory; posterior];
			xstory  = [ xstory params ];
			fprintf('%4s \t posterior: %.8f \t ',['#' num2str(optimvalues.iteration)],-optimvalues.fval)
			fprintf('\n');
			if mod(optimvalues.iteration+1,2) == 0
				% then shows estimated parameters
				load ve_names;
				for i1=1:length(params)
					fprintf('%10s : %.8f \n ',ve_names{i1},params(i1))
				end
				fprintf('\n');
            end
            
            try
    			save([M_.fname '_mle_estimates_temp.mat'], 'posterior', 'datavec', 'fstory', 'xstory', 'params');
            catch
            end
            updte=1;
		else
			fprintf('.')
		end
		
	end

	
	
end