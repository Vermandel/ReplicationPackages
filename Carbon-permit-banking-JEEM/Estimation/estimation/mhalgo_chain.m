
function [] = mhalgo_chain(posterior_func,theta_mode,P,c,options,toss,toss_normal,chain_nb,ndraws,load_switch,obs,e_obj,oo_,M_,options_,bayestopt_)

	save_every_n_draws = 50;
	filename = [M_.fname '_mh_chain_' num2str(chain_nb) '.mat'];
	
	% if parallel pool, update some folder
	%if ~isempty(getCurrentTask())
	%	gg =which('dynare');
	%	addpath([gg(1:end-8) 'qz\']);
	%	addpath('estimation');
	%	which('dsge_llk')
	%end
	
	nparams = size(theta_mode,1);

	if load_switch
		load(filename)
		accept			= accept_ratio(end)*length(accept_ratio);
		theta_history 	= [theta_history zeros(nparams,ndraws)];
		fval_history 	= [fval_history zeros(1,ndraws)];
		accept_ratio 	= [accept_ratio, zeros(1,ndraws)];
		fval 			= fval_history(:,indx);
		theta 			= theta_history(:,indx);
		toss			= [toss rand(1,ndraws)];
		toss_normal		= [toss_normal randn(nparams,ndraws)];
		indx_init 		= indx;
  		[~,ee] 			= feval(posterior_func,theta,obs,e_obj,oo_,M_,options_,bayestopt_);		
	else
		if exist(filename,'file')
			delete(filename);
        end
        
 		theta = theta_mode;
  		[fval,ee] = feval(posterior_func,theta,obs,e_obj,oo_,M_,options_,bayestopt_);
		fval = -fval;

		theta_history = zeros(nparams,ndraws);
		fval_history = zeros(1,ndraws);
		accept_ratio = zeros(1,ndraws);
		success		 = nan(1,ndraws);

		indx_init = 0;
		accept = 0;
		
	end
	
	% using the same sequence of shocks to speed up
	%ee = zeros(size(obs));

	for indx = indx_init+1:indx_init+ndraws

		
		theta_new = theta + c^2*P'*toss_normal(:,indx);
		[fval_new,ee_new] = feval(posterior_func,theta_new,obs,e_obj,oo_,M_,options_,bayestopt_,ee);
		fval_new = -fval_new;
		
		fprintf('chain(#%d) %.2f%%',chain_nb,indx/(indx_init+ndraws)*100);
		
            
		r = min(exp(fval_new-fval),1);     % probability of acceptance
		% check if new draw is accepted or rejected
			
		if (toss(indx)<r && r >= 0)
			fprintf('\t %.3f > %.3f (llk:%8.2f > %8.2f) ',r,toss(indx),fval_new,fval);
			fprintf(' (ACCEPT) ');
			fval = fval_new;
			theta = theta_new;
			accept = accept+1;
			success(indx)=1;
			ee = ee_new;
		else
			fprintf('\t %.3f < %.3f (llk:%8.2f < %8.2f) ',r,toss(indx),fval_new,fval);
			fprintf(' (reject) ');
			success(indx)=0;
		end
		accept_ratio(indx)=accept/indx;
	 %       disp([ 'Current ar = ' num2str(accept_ratio(indx),2) ])
		fprintf(' (AR=%.2f) \n',accept_ratio(indx));
		theta_history(:,indx)=theta;
		fval_history(indx)=fval;
		%else
		%	save hey;
		%	error('hey')
		%end
		
		if mod(indx,save_every_n_draws)==0
			save(filename,'fval_history','theta_history','c','accept_ratio','accept','indx', 'toss', 'toss_normal')
		end
		
	end
	save(filename,'fval_history','theta_history','c','accept_ratio','accept','indx', 'toss', 'toss_normal')
	accept_average = accept/ndraws;
end



