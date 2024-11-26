% options for drawing
numC = 3;
numR = 5;
np   = length(theta);
do_plots=1;
ndrop = 500;

% if function exists, then export data
exporting = exist('matrix_to_txt.m')>0;

% get ids of metropolis chains
list_MC=dir(['*' M_.fname '_mh_chain_' '*.mat']);
list_MC=char(list_MC.name);
idMC   = list_MC(:,end-4);
nMC    = length(idMC);


title_chains = char(['Chain ' idMC(1)]);
load([M_.fname '_mh_chain_' idMC(1) '.mat']);
idx = find(fval_history~=0);
%drop first draws
idx2=ceil(length(idx)*options_.mh_drop);
theta_all = theta_history(:,idx(idx2:end));
idmax = idx(end);
if nMC > 1
	for i1=2:nMC
		load([M_.fname '_mh_chain_' idMC(i1) '.mat']);
		idx = find(fval_history~=0);
		theta_all = [theta_all theta_history(:,idx)];
		idmax = min(idmax,idx(end));
		title_chains = char(title_chains,['Chain ' num2str(i1)]);
	end
end


% get maximum mode
load([M_.fname '_mh_chain_' idMC(1) '.mat']);
fval_tot=fval_history(1:idmax);
if nMC > 1
	for i1 = 2:nMC
		load([M_.fname '_mh_chain_' idMC(i1) '.mat']);
		fval_tot=[fval_tot;fval_history(1:idmax)];
	end
end
amax=max(fval_tot(:));
[maxrow,maxcol]=find(fval_tot==amax);
load([M_.fname '_mh_chain_' idMC(maxrow(1)) '.mat']);
mode_metropolis = theta_history(:,maxcol(1));



if do_plots==1
    indx = size(theta_all,2);
	

	thePosteriorsX = [];	% for exporting
	thePosteriorsY = [];
	thePosteriorMnX = [];
	thePosteriorMnY = [];
	thePriorsX = [];		% for exporting
	thePriorsY = [];
	theplottot=1;
	for fig = 1:ceil(np/(numC*numR))
	theplot = 1;
		figure(100+fig)
		for ii=1:(numC*numR)
			if theplottot <= np
				subplot(numR,numC,theplot)
				
				ksdensity(theta_all(theplottot,2:indx)); 

				hold on;
							
				plab = e_obj.theta.names(theplottot,:);
				
				% get id in dynare
				idd = strmatch(deblank(e_obj.theta.names(theplottot,:)),bayestopt_.name,'exact');
				
				if plab(1:2)=='ST'; axis tight; end        

				[xp,fp] = myplot_priors2(bayestopt_.pshape(idd),bayestopt_.p3(idd),bayestopt_.p4(idd),bayestopt_.p6(idd),bayestopt_.p7(idd),bayestopt_.name(idd));
				plot_lines_green(mode_metropolis(theplottot));
				
				if exporting && exist('mean_metropolis')
					[fi,xi] = ksdensity(theta_all(theplottot,2:indx));
					thePosteriorsX = [thePosteriorsX  xi'];
					thePosteriorsY = [thePosteriorsY  fi'];
					thePriorsX = [thePriorsX  xp'];
					thePriorsY = [thePriorsY  fp'];
					thePosteriorMnX = [thePosteriorMnX [mean_metropolis(theplottot);mean_metropolis(theplottot)]];
					thePosteriorMnY = [thePosteriorMnY [min(fi);max(fi)+.1*abs(max(fi)-min(fi))]];
				end
				
				v = axis;
				if ~isempty(find(e_obj.theta.is_sd==theplottot)) % std shocks
                    idd = strmatch(deblank(e_obj.theta.names(theplottot,:)),M_.exo_names,'exact');
					if isfield(M_,'exo_names_tex') && ~isempty(M_.exo_names_tex2(idd,:)) && ~isempty(M_.exo_names_long2(idd,:))
						title([  M_.exo_names_long2(idd,:) ' ($' deblank(M_.exo_names_tex2(idd,:)) '$)'],'interpreter','latex')	
					else
						title(e_obj.theta.names(theplottot,:))					
					end
				else % structural parameter
                    idd = strmatch(deblank(e_obj.theta.names(theplottot,:)),M_.param_names,'exact');
					if isfield(M_,'param_names_tex') && ~isempty(M_.param_names_tex2(idd,:)) && ~isempty(M_.param_names_long2(idd,:))
						title([  M_.param_names_long2(idd,:) ' ($' deblank(M_.param_names_tex2(idd,:)) '$)'],'interpreter','latex')	
					else
						title(e_obj.theta.names(theplottot,:))					
					end
				end
				xlim([v(1) max(0.02,v(2))])
			end
			theplot = theplot+1;
			theplottot=theplottot+1;
		end
		legend('Posterior','Prior','mode')
	end
	
	% exporing this to a text file
	if exporting
		matrix_to_txt([thePriorsX thePriorsY],'priors.txt',char(strcat(char(e_obj.theta.names),'_x'),strcat(char(e_obj.theta.names),'_y') ))	
		matrix_to_txt([thePosteriorsX thePosteriorsY],'posteriors.txt',char(strcat(char(e_obj.theta.names),'_x'),strcat(char(e_obj.theta.names),'_y') ))
		matrix_to_txt([thePosteriorMnX thePosteriorMnY],'posteriors_mean.txt',char(strcat(char(e_obj.theta.names),'_x'),strcat(char(e_obj.theta.names),'_y') ))
	end
	
	figure(202)
	% plot log likelihood
	subplot(2,1,1)
	load([M_.fname '_mh_chain_' idMC(1) '.mat']);
	idx = find(fval_history~=0);
    plot(fval_history(1:idmax),'LineWidth',.5)
	fval_tot=fval_history(1:idmax);
    title([ 'The posterior draws ' num2str(max(fval_tot),'%0.4f')])
    if nMC > 1
		hold on;
		for i1 = 2:nMC
			load([M_.fname '_mh_chain_' idMC(i1) '.mat']);
			plot(fval_history(1:idmax),'LineWidth',.5)
			fval_tot=[fval_tot;fval_history(1:idmax)];
		end
		hold off;
		%iddx = find(fval_tot>100000 & fval_tot<100000); % remove extremes
		ylim([mean(fval_tot,'all')-3.5*std(fval_tot,[],'all') mean(fval_tot,'all')+3.5*std(fval_tot,[],'all')])
	end
	% plot acceptance rate
	subplot(2,1,2)
	load([M_.fname '_mh_chain_' idMC(1) '.mat']);
	idx = find(fval_history~=0);
    plot(accept_ratio(1:idmax),'LineWidth',.5)
    title(['The acceptance rate'])
    if nMC > 1
		hold on;
		for i1 = 2:nMC
			load([M_.fname '_mh_chain_' idMC(i1) '.mat']);
			plot(accept_ratio(1:idmax),'LineWidth',.5)
			fval_tot=[fval_tot;fval_history(1:idmax)];
		end
		hold off;
	end
	
	theplottot=1;
	numC = 1;
	numR = 4;
	for fig = 1:ceil(np/(numC*numR))
	theplot = 1;
		figure(202+fig)
		for ii=1:(numC*numR)
			if theplottot <= np

				subplot(numR,numC,theplot)
				load([M_.fname '_mh_chain_' idMC(1) '.mat']);
				plot(theta_history(theplottot,1:idmax),'LineWidth',.5); 
				hold on; axis tight
				plab = e_obj.theta.names(theplottot,:);
				if plab(1:2)=='ST'; axis tight; end
				title(e_obj.theta.names(theplottot,:));
				if nMC > 1
					for i1=2:nMC
						load([M_.fname '_mh_chain_' idMC(i1) '.mat']);
						plot(theta_history(theplottot,1:idmax),'LineWidth',.5); 
					end
				end
				hold off;
			end
			legend(title_chains)
			theplot = theplot+1;
			theplottot=theplottot+1;
		end
    end
    
end

cut_percent_obs=0.2;

% get truncated parameters
load([M_.fname '_mh_chain_' idMC(1) '.mat']);
theta_history_cut = theta_history(:,(cut_percent_obs*idmax+1:idmax));
fval_history_cut = fval_history(cut_percent_obs*idmax+1:idmax);
if nMC > 1
	for i1 = 2:nMC
		load([M_.fname '_mh_chain_' idMC(i1) '.mat']);
		theta_history_cut = [theta_history_cut theta_history(:,(cut_percent_obs*idmax+1:idmax))];
		fval_history_cut  = [fval_history_cut fval_tot(cut_percent_obs*idmax+1:idmax)];
	end
end


trunc=0.1:0.1:0.9;

md = mhm_marginal_density(theta_history_cut,fval_history_cut,trunc);








ndraws_cut = size(theta_history_cut,2);
disp('-----------')
disp(['Number of runs is ' num2str(idmax*nMC)])
disp(['Objective @ mode  ' num2str(amax,6)])
disp(['Log Marginal data density is ' num2str(mean(md))])
disp(' ')
fprintf('%10s \t %10s \t %10s \t %10s \t %10s \t %10s %12s \n','PARAMETER','DISTRIB','MEAN','STD','MODE','MEAN','[0.050;0.950]')
jj = 0;
mean_metropolis = nan(size(mode_metropolis));
for ii=1:length(e_obj.thet_ids)
	
    idd=e_obj.thet_ids(ii);
	fprintf('%10s \t ',bayestopt_.name{idd})
	
	theprior=bayestopt_.p1(idd);
	switch bayestopt_.pshape(idd)
		case 6
			theprior = 'InvGam2';
		case 1
			theprior = 'Beta     ';
		case 2
			theprior = 'Gamma    ';
		case 3
			theprior = 'Normal   ';
		case 4
			theprior = 'InvGam1   ';
	end
	fprintf('%10s \t ',lower(theprior))
	
	fprintf('%10.2f \t ',bayestopt_.p1(idd))
	fprintf('%10.2f \t ',bayestopt_.p2(idd))
	
	fprintf('%10.3f \t ',mode_metropolis(ii))
	
	% draw the CDF and get percentiles
	[y,x]=ecdf(theta_history_cut(ii,1:indx));
	id_lb = find(y>=0.049999 & y<0.08);
	id_ub = find(y>=0.949999 & y<0.98);
	id_mb = find(y>=0.499999 & y<0.53);
	fprintf('%10.3f [%.3f;%.3f]\n',x(id_mb(1)),x(id_lb(1)),x(id_ub(1)))
    mean_metropolis(ii) = x(id_mb(1));
end
[llk,ee,res,yy] = dsge_llk(mean_metropolis,obs,e_obj,oo_,M_,options_,bayestopt_);
eval(['save ' M_.fname '_mean.mat mean_metropolis mode_metropolis  e_obj ee;'])
