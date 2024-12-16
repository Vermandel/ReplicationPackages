
% this file simulate climate 
%----------------------------------------------------------------
% 1. Options
%----------------------------------------------------------------

EP_options.load_theta 			= 'MHmode';				% 
EP_options.estimation_file_name = 'run_estimation';		% name of the file to load estimated parameter values
options_.simul.maxit            = 220;					% number of iterations for perfect foresight solver

%----------------------------------------------------------------
% 2. Load Files
%----------------------------------------------------------------

% add additional variables to simulations :
@#define 		SIMULATIONS = 1
% load core model file:
@#include 		"model_file.mod"
% load additional matlab functions
@#includepath 	"estimation_v6"
% call a routine that loads estimated parameters
@#include 		"load_parameters.mod"


%----------------------------------------------------------------
% 3. Computation
%----------------------------------------------------------------


spfm_exo_simul 	= IF_EP_init_shocks(size(exo_filtered,1),M_,options_);
% check observable 
[endo_simul,endo_simul0] = simult_EP(exo_filtered(2:end,:),oo_,M_,options_,e_obj.id.shocks,endo_filtered,155);
endo_simulx 						= endo_simul0;
endo_simulx(:,1:size(endo_simul,2)) = endo_simul;
endo_simul_ts  						= dseries(endo_simulx',filtered_errors_ts.dates(1)-1,M_.endo_names);
exo_filtered_ts						= dseries(exo_filtered,filtered_errors_ts.dates(1)-1,M_.exo_names);


%% options
% then solve extended path
Sim_dates = dates('2023Q3'):dates('2123Q3'); 
dummy_surprises  					= zeros(1,M_.exo_nbr);
dummy_surprises(surprise_shocks) 	= 1;
exo_simul							= exo_filtered_ts(Sim_dates(1:end)).data;


% end of sample date is:
T_ = size(exo_filtered,1);

% filtered data are:
y_ = endo_simul_ts(Sim_dates(1:end)).data';

%% <-------------------------------- aa -------------------------------------->

T2	=	1985:1/4:2200;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute IRFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Twindow     = options_.expectation_window;				% time horizon of IRF
nirfs 		= 500;										% number of draws
t 			= 2;										% date of IRF

% find 'surprise' shocks
err_ids 	= find(sum(M_.Sigma_e,1)~=0);
Sigma2  	= M_.Sigma_e(err_ids,err_ids);

% saving results in this
GIRFs		= nan([size(y_,1) Twindow+2 nirfs length(err_ids)]);
GIRFs_pp	= nan([size(y_,1) Twindow+2 nirfs length(err_ids)]);

f 		= waitbar(0,'Computing...');
for it_ = 1:nirfs

	% avoid same seeds across burns
    rng(it_^2*t)
	waitbar(it_/nirfs,f,['Compute IRFs ']);
	

	% load shocks
	oo_.endo_simul	 			= y_(:,(t-1):(t+options_.expectation_window));
	options_.periods			= size(oo_.endo_simul,2)-2;
	
	%% add contemporaneous shocks
	oo_.exo_simul			= zeros(size(oo_.endo_simul,2),M_.exo_nbr);
	%oo_.exo_simul(:,idx) 	= randn(size(oo_.endo_simul,2),length(idx))*(chol(M_.Sigma_e(idx,idx))');
	baseline_exo 			= randn(1,length(err_ids))*chol(Sigma2)';
	
		
	% get baseline 
	oo_.exo_simul(2,err_ids)	= baseline_exo;
	%oo_ 						= perfect_foresight_solver_core(M_,options_,oo_);
	perfect_foresight_solver;
	if ~oo_.deterministic_simulation.status; continue; end
	baseline_endo	 = oo_.endo_simul;
	
	% add shock
	for i0 = 1:length(err_ids)
		% reboot path
		oo_.endo_simul = baseline_endo;
		% size of impulse
		ex_ 	= zeros(1,length(err_ids));
		ex_(i0) = 1;    % 1 sd shock
		% feed with current shocks
		oo_.exo_simul(2,err_ids)	= baseline_exo + ex_*chol(Sigma2)';
		% solve
		%oo_ 						= perfect_foresight_solver_core(M_,options_,oo_);
		perfect_foresight_solver;
		if oo_.deterministic_simulation.status
			% compute distance
			GIRFs(:,:,it_,i0)    = oo_.endo_simul-baseline_endo;
			GIRFs_pp(:,:,it_,i0)  = GIRFs(:,:,it_,i0);
			
			% if non zero initial state -> express in pp from initial state 
			nozeros = find(baseline_endo(:,1)~=0);
			for i1 = 1:length(nozeros) 
				GIRFs_pp(nozeros(i1),:,it_,i0) = GIRFs_pp(nozeros(i1),:,it_,i0)/baseline_endo(nozeros(i1))*100;
			end
		end
	end
end
close(f)


TT			= 2:(size(GIRFs,2)-1);
yall    	= (1:Twindow)';
yall_pp 	= (1:Twindow)';
colnames 	= char('time');
IC 			= [5 50 95];

for i1 = 1:length(err_ids)
    
    notnans = find(~isnan(GIRFs(1,1,:,i1)));

    % get mean of level IRFs
	y1 = quantile(GIRFs(:,:,notnans,i1),IC/100,3);
	y1 = y1(:,TT,:);
    % get mean of relative IRFs
	y2 = quantile(GIRFs_pp(:,:,notnans,i1),IC/100,3);
	y2 = y2(:,TT,:);
	y2 = (abs(y2)>options_.dynatol.f).*y2;
    % save endo names
    endo_names = M_.endo_names;
    exo_names  = M_.exo_names;
    endo_names_long = M_.endo_names_long;
    exo_names_long  = M_.exo_names_long;
    endo_names_tex  = M_.endo_names_tex;
    exo_names_tex   = M_.exo_names_tex;
    %
	for ics = 1:length(IC)
		yall    = [yall y1(:,:,ics)'];
		yall_pp = [yall_pp y2(:,:,ics)'];
		colnames = char(colnames,strcat(char(endo_names),['_' char(M_.exo_names{err_ids(i1)}) '_' int2str(IC(ics))] ));
	end
    
 
end



varnames=char('lny','pi100','rr100');
nn=size(varnames,1);
nt = size(GIRFs_pp,2);
T2 = 1:(nt-1)';
nx = size(GIRFs_pp,4);
figure('Name','Figure 3: Decomposition of detrended output during the transition');
ni = 1;
for i2 = 1:nx
	for i1 = 1:nn
		subplot(nx,nn,ni)
			ni = ni + 1;
			idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
			y1 = squeeze(quantile(GIRFs(idx,:,:,i2),[0.025 0.5 0.975],3));
			
			plot(T2,y1(T2,2),'linewidth',2,'Color',[50, 168, 82]/255)
			
			title(['$' M_.endo_names_tex{idx} '$ ' M_.endo_names_long{idx}],'Interpreter','latex');
			hold on;
			plot(T2,y1(T2,1),'linewidth',1,'Color',[140, 135, 135]/255)
			plot(T2,y1(T2,3),'linewidth',1,'Color',[140, 135, 135]/255)
			hold off;
			grid on
			if i1==1
				ylabel(M_.exo_names_long{i2})
			end
	end
end

