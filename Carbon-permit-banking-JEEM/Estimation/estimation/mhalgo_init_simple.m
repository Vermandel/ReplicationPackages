
function [] = mhalgo_init_simple(posterior_func,theta_mode,Hessian,nb_chains,ndraws,load_switch,jump,obs,e_obj,oo_,M_,options_,bayestopt_)
	
	if ~exist('nb_chains')
		nb_chains = 1;
	end
	
	if size(theta_mode,1)<size(theta_mode,2)
		theta_mode = theta_mode';
	end
	
	if load_switch
		disp('Loading previous results...')
		idmax = [];
		disp('Number of iterations per chains:')
		% compare the chains, and get the same number of iterations
		for i1=1:nb_chains
			load([M_.fname '_mh_chain_' num2str(i1) '.mat'],'fval_history');
			id2 = find(fval_history~=0);
			disp(['Chain(' num2str(i1) ') : ' num2str(max(id2)) ]);
			idmax = [idmax max(id2)];
		end
		% among all the chains, selects the lowest number of iterations
		idmax = min(idmax);
		disp('Removing extra iterations...')	
		% Time to clean chains
		for i1=1:nb_chains
			load([M_.fname '_mh_chain_' num2str(i1) '.mat']);
			fval_history 	= fval_history(1:idmax);
			theta_history	= theta_history(:,1:idmax);
			toss_normal		= toss_normal(:,1:idmax);
			toss			= toss(1:idmax);
			indx			= idmax;
			accept_ratio	= accept_ratio(1:idmax);
			c				= jump;
			save([M_.fname '_mh_chain_' num2str(i1) '.mat'],'fval_history','theta_history','toss_normal','toss','indx','accept_ratio','c','accept');
		end
		
	end
	
	% generate random numbers
	strc = 's1';
	if nb_chains > 1
		for i1 = 2:nb_chains
			strc = [strc ',s' num2str(i1)];
		end
	end
	% generate random generator objects
	% and store them into s1,s2,...
	eval(['[' strc '] = RandStream.create(''mrg32k3a'',''NumStreams'',nb_chains);']) ;

	P = chol(inv(Hessian));
	% hit target acceptance rate.
	options = optimset('display','iter','TolFun',1e-2,'TolX',1e-2,'MaxIter',1000,'MaxFunEvals',10000);
	c=jump;
	
	nparams = size(theta_mode,1);
	for i1 = 1:nb_chains
		toss_mat{i1}		= eval(['rand(s' num2str(i1) ',1,ndraws)']);
		toss_normal_mat{i1}	= eval(['randn(s' num2str(i1) ',nparams,ndraws)']);
	end
	
	
	% execute chains
	for i0 = 1:nb_chains
		%% persistent variables are troublesome in parallel computing
		% initialize them
		%if  ~isempty(getCurrentTask()) % if parallel
			lnprior = priordens(theta_mode(e_obj.thet_ids),bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4,1);
			% one at the end starts the initialization
		%end		
		mhalgo_chain(posterior_func,theta_mode,P,c,options,toss_mat{i0},toss_normal_mat{i0},i0,ndraws,load_switch,obs,e_obj,oo_,M_,options_,bayestopt_);
	end
end



