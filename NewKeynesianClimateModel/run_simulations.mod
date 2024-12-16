
% this file simulate climate 
%----------------------------------------------------------------
% 1. Options
%----------------------------------------------------------------

%

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

steady;
%----------------------------------------------------------------
% 3. Computation
%----------------------------------------------------------------

% identify what is a sur
surprise_shocks 					= logical(zeros(M_.exo_nbr,1));
surprise_shocks(e_obj.id.shocks) 	= logical(1);

% simulate the model under estimated scenario
[endo_simul,endo_simul0] 			= simult_EP(exo_filtered(2:end,:),oo_,M_,options_,e_obj.id.shocks,endo_filtered,215);
endo_simulx 						= endo_simul0;
endo_simulx(:,1:size(endo_simul,2)) = endo_simul;
endo_simul_ts  						= dseries(endo_simulx',filtered_errors_ts.dates(1)-1,M_.endo_names);
exo_filtered_ts						= dseries(exo_filtered,filtered_errors_ts.dates(1)-1,M_.exo_names);


% then solve extended path
Sim_dates = dates('2023Q4'):dates('2773Q3'); 
dummy_surprises  					= zeros(1,M_.exo_nbr);
dummy_surprises(surprise_shocks) 	= 1;
exo_simul							= exo_filtered_ts(Sim_dates(2:end)).data;
SampleSize 							= 400; % length of the window of simulation to re-update
% Save initial endogenous variable path before making change in carbon tax path
eval('y0 = endo_simul_ts(exo_filtered_ts.dates(1):(Sim_dates(1)-1)).data;')

% make a backup of parameter before adjusting expectations
M_bak=M_;


%% BUSINESS AS USUAL SCENARIO : PHI*tau = 0*tau
M_bau 				= M_bak;
M_bau.params(strmatch('Et',M_.param_names,'exact')) = 0;		% psi*tau = 0*tau : no carbon tax
endo_simul_bau		= endo_simul_ts(Sim_dates).data';		
M_ 					= M_bau;
oo_.endo_simul      = endo_simul_ts(Sim_dates).data';			% use as inital conditions the path under baseline 2023Q4:2777Q3
oo_.exo_simul		= exo_filtered_ts(Sim_dates(1:end)).data;	% take the shocks/path of carbon tax
options_.periods 	= size(oo_.endo_simul,2)-2;					% rescale the update of the extended simulation over period 2023Q4:2777Q3
oo_.endo_simul(:,end) = feval([M_.fname '.steadystate'],oo_.steady_state, zeros(M_.exo_nbr,1), M_.params);
perfect_foresight_solver;										% solve deterministic path
endo_simul_bau		= oo_.endo_simul;
% add shocks:
endo_simul_bau 		= [y0' simult_EP_solve_path(exo_simul,endo_simul_bau,options_.expectation_window,oo_,M_bau,options_,dummy_surprises,SampleSize)];
endo_simul_bau_ts  	= dseries(endo_simul_bau',filtered_errors_ts.dates(1)-1,M_.endo_names);


%% NET ZERO SCENARIO: PHI*tau = 1*tau
M_nz0 				= M_bak;
M_nz0.params(strmatch('Et',M_.param_names,'exact')) = 1;
endo_simul_nz0		= endo_simul_ts(Sim_dates).data';
M_ 					= M_nz0;
oo_.endo_simul      = endo_simul_ts(Sim_dates).data';
oo_.exo_simul		= exo_filtered_ts(Sim_dates(1:end)).data;
options_.periods 	= size(oo_.endo_simul,2)-2;
perfect_foresight_solver;
endo_simul_nz0		= oo_.endo_simul;
endo_simul_nz0 		= [y0' simult_EP_solve_path(exo_simul,endo_simul_nz0,options_.expectation_window,oo_,M_nz0,options_,dummy_surprises,SampleSize)];
endo_simul_nz0_ts  	= dseries(endo_simul_nz0',filtered_errors_ts.dates(1)-1,M_.endo_names);


%% BASELINE SCENARIO : PHI*tau = 0.53*tau
M_bsl 				= M_bak;
endo_simul_bsl		= endo_simul_ts(Sim_dates).data';
M_ 					= M_bsl;
oo_.endo_simul      = endo_simul_ts(Sim_dates).data';
oo_.exo_simul		= exo_filtered_ts(Sim_dates(1:end)).data;
options_.periods 	= size(oo_.endo_simul,2)-2;
oo_.endo_simul(:,end) = feval([M_.fname '.steadystate'],oo_.steady_state, zeros(M_.exo_nbr,1), M_.params);
perfect_foresight_solver;
endo_simul_bsl		= oo_.endo_simul;
endo_simul_bsl 		= [y0' simult_EP_solve_path(exo_simul,endo_simul_bsl,options_.expectation_window,oo_,M_bsl,options_,dummy_surprises,SampleSize)];
endo_simul_bsl_ts  	= dseries(endo_simul_bsl',filtered_errors_ts.dates(1)-1,M_.endo_names);

% PLOT FIGURE 3
nx = 3;ny = 3;varnames=char('lny','pi100','rr100','lny_n','output_gap','rrn100','tau_USD','E');

nn=size(varnames,1);
T2 = dseries_to_num(endo_simul_nz0_ts)';
figure('Name','Figure 2: Model-implied projections based on alternative control rates of emissions');
for i1 =1:nn
	subplot(nx,ny,i1)
	
		idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
		y1 = endo_simul_nz0(idx,:);
		y2 = endo_simul_bau(idx,:);
		y3 = endo_simul_bsl(idx,:);
		
		plot(T2,y1,':','linewidth',2,'Color',[50, 168, 82]/255)
		
		title(['$' M_.endo_names_tex{idx} '$ ' M_.endo_names_long{idx}],'Interpreter','latex')

	
		hold on;
		xlim([2023.5 2100])
		plot(T2,y2,'linewidth',2,'Color',[194, 68, 43]/255)
		plot(T2,y3,'--','linewidth',2,'Color',[50, 96, 168]/255)
			hold off;
		grid on

end
legend('Paris-Agreement','Laissez-faire','Estimated tax path')


nx = 1;ny = 4;varnames=char('lny','IS','rotemberg_to_gdp','abatement_to_gdp');
nn=size(varnames,1);
T2 = dseries_to_num(endo_simul_nz0_ts)';
figure('Name','Figure 3: Decomposition of detrended output during the transition');
for i1 =1:nn
	subplot(nx,ny,i1)
	
		idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
		y1 = endo_simul_nz0(idx,:);
		y2 = endo_simul_bau(idx,:);
		y3 = endo_simul_bsl(idx,:);
		
		plot(T2,y1,':','linewidth',2,'Color',[50, 168, 82]/255)
		
		title(['$' M_.endo_names_tex{idx} '$ ' M_.endo_names_long{idx}],'Interpreter','latex')

	
		hold on;
		xlim([2023.5 2100])
		plot(T2,y2,'linewidth',2,'Color',[194, 68, 43]/255)
		plot(T2,y3,'--','linewidth',2,'Color',[50, 96, 168]/255)
			hold off;
		grid on

end
legend('Paris-Agreement','Laissez-faire','Estimated tax path')



nx = 3;ny = 2;varnames=char('pi_tot','pi_hat_s','pi_hat_m','pi_hat_g','pi_hat_x');
nn=size(varnames,1);
T2 = dseries_to_num(endo_simul_nz0_ts)';
figure('Name','Figure 4: Decomposition of inflation during the green transition');
for i1 =1:nn
	subplot(nx,ny,i1)
	
		idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
		y1 = endo_simul_nz0(idx,:);%- endo_simul_bsl(idx,:);
		y2 = endo_simul_bau(idx,:);%- endo_simul_bsl(idx,:);
		y3 = endo_simul_bsl(idx,:);
		
		plot(T2,y1,':','linewidth',2,'Color',[50, 168, 82]/255)
		
		title(['$' M_.endo_names_tex{idx} '$ ' M_.endo_names_long{idx}],'Interpreter','latex')

	
		hold on;
		xlim([2023.5 2100])
		plot(T2,y2,'linewidth',2,'Color',[194, 68, 43]/255)
		plot(T2,y3,'--','linewidth',2,'Color',[50, 96, 168]/255)
			hold off;
		grid on

end
legend('Paris-Agreement','Laissez-faire','Estimated tax path')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot TRENDS
nx = 2;ny = 3;varnames=char('Z','L','SIG','THETA1','pi_bar','tau_USD');
nn=size(varnames,1);
T2 = dseries_to_num(endo_simul_nz0_ts)';
figure;
for i1 =1:nn
	subplot(nx,ny,i1)
		idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
		y1 = endo_simul_nz0(idx,:);
		plot(T2,y1,':','linewidth',2,'Color',[50, 168, 82]/255)
		title([ M_.endo_names_long{idx}],'Interpreter','latex')
		xlim([1985 2100])
        grid;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE: DETERMINISTIC VERSUS STOCHASTIC PATH
endo_simul0_ts  					= dseries(endo_simul0',filtered_errors_ts.dates(1)-1,M_.endo_names);
T0                                  = dseries_to_num(endo_simul0_ts)';
figure;
subplot(1,4,1)
    plot(T2,endo_simul_nz0_ts.dy.data,T0,endo_simul0_ts.dy.data);
    xlim([1971 2023.25]);grid on;
    title('Output growth','Interpreter','latex')
subplot(1,4,2)
    plot(T2,endo_simul_nz0_ts.pi.data,T0,endo_simul0_ts.pi.data);
    xlim([1971 2023.25]);grid on;
    title('Inflation','Interpreter','latex')
subplot(1,4,3)
    plot(T2,endo_simul_nz0_ts.r.data,T0,endo_simul0_ts.r.data);
    xlim([1971 2023.25]);grid on;
    title('Interest rate','Interpreter','latex')
subplot(1,4,4)
    plot(T2,endo_simul_nz0_ts.rr_real.data,T0,endo_simul0_ts.rr_real.data);
    xlim([1971 2023.25]);grid on;
    title('Real rate','Interpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% CHANGING RULES : NATURAL RATE MONETARY POLICY RULE
% r(t)/(pi_bar(t)*rn(t)) Nat Rule
M_nrule					= M_bak;
M_nrule.params(strmatch('Et',M_.param_names,'exact')) = 1;
M_nrule.params(strmatch('MP_r_n',M_.param_names,'exact')) = 1;
endo_simul_nrule_nz0	= endo_simul_ts(Sim_dates).data';
M_ 						= M_nrule;
oo_.endo_simul      	= endo_simul_ts(Sim_dates).data';
oo_.exo_simul			= exo_filtered_ts(Sim_dates(1:end)).data;
options_.periods 		= size(oo_.endo_simul,2)-2;
[oo_]=perfect_foresight_solver(M_, options_, oo_);
endo_simul_nrule_nz0		= oo_.endo_simul;
endo_simul_nrule_nz0 		= [y0' simult_EP_solve_path(exo_simul,endo_simul_nrule_nz0,options_.expectation_window,oo_,M_nrule,options_,dummy_surprises,SampleSize)];
endo_simul_nrule_nz0_ts  	= dseries(endo_simul_nrule_nz0',filtered_errors_ts.dates(1)-1,M_.endo_names);
%% BAU
M_nrule.params(strmatch('Et',M_.param_names,'exact')) = 0;
endo_simul_nrule_bau	= endo_simul_ts(Sim_dates).data';
M_ 						= M_nrule;
oo_.endo_simul      	= endo_simul_ts(Sim_dates).data';
oo_.exo_simul			= exo_filtered_ts(Sim_dates(1:end)).data;
options_.periods 		= size(oo_.endo_simul,2)-2;
[oo_]=perfect_foresight_solver(M_, options_, oo_);
endo_simul_nrule_bau		= oo_.endo_simul;
endo_simul_nrule_bau 		= [y0' simult_EP_solve_path(exo_simul,endo_simul_nrule_bau,options_.expectation_window,oo_,M_nrule,options_,dummy_surprises,SampleSize)];
endo_simul_nrule_bau_ts  	= dseries(endo_simul_nrule_bau',filtered_errors_ts.dates(1)-1,M_.endo_names);




%%%%%%%%%% CHANGING RULES
% y(t)/(Ybar*L(t)*Z(t)) STEADY STATE RULE
M_yrule					= M_bak;
M_yrule.params(strmatch('Et',M_.param_names,'exact')) = 1;
M_yrule.params(strmatch('MP_y_n',M_.param_names,'exact')) = 0;
endo_simul_yrule_nz0	= endo_simul_ts(Sim_dates).data';
M_ 						= M_yrule;
oo_.endo_simul      	= endo_simul_ts(Sim_dates).data';
oo_.exo_simul			= exo_filtered_ts(Sim_dates(1:end)).data;
options_.periods 		= size(oo_.endo_simul,2)-2;
[oo_]=perfect_foresight_solver(M_, options_, oo_);
endo_simul_yrule_nz0		= oo_.endo_simul;
endo_simul_yrule_nz0 		= [y0' simult_EP_solve_path(exo_simul,endo_simul_yrule_nz0,options_.expectation_window,oo_,M_yrule,options_,dummy_surprises,SampleSize)];
endo_simul_yrule_nz0_ts  	= dseries(endo_simul_yrule_nz0',filtered_errors_ts.dates(1)-1,M_.endo_names);
%%% BAU
M_yrule.params(strmatch('MP_y_n',M_.param_names,'exact')) = 0;
M_yrule.params(strmatch('Et',M_.param_names,'exact')) = 0;
endo_simul_yrule_bau	= endo_simul_ts(Sim_dates).data';
M_ 						= M_yrule;
oo_.endo_simul      	= endo_simul_ts(Sim_dates).data';
oo_.exo_simul			= exo_filtered_ts(Sim_dates(1:end)).data;
options_.periods 		= size(oo_.endo_simul,2)-2;
[oo_]=perfect_foresight_solver(M_, options_, oo_);
endo_simul_yrule_bau		= oo_.endo_simul;
endo_simul_yrule_bau 		= [y0' simult_EP_solve_path(exo_simul,endo_simul_yrule_bau,options_.expectation_window,oo_,M_yrule,options_,dummy_surprises,SampleSize)];
endo_simul_yrule_bau_ts  	= dseries(endo_simul_yrule_bau',filtered_errors_ts.dates(1)-1,M_.endo_names);





varnames	= char('lny','pi100','rr100');
nn			= size(varnames,1); nx = 2; ny = 3;
figure('Name','FIGURE 6. The impact of climate transition under alternative monetary policy rules');
for i1 =1:(2*nn)
	subplot(nx,ny,i1)
		
		if i1 < nn+1
		idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
		y1 = endo_simul_nz0(idx,:);
		y2 = endo_simul_nrule_nz0(idx,:);
		y3 = endo_simul_yrule_nz0(idx,:);
		else
		idx = strmatch(deblank(varnames(i1-nn,:)),M_.endo_names,'exact');
		y1 = endo_simul_bau(idx,:);
		y2 = endo_simul_nrule_bau(idx,:);
		y3 = endo_simul_yrule_bau(idx,:);
		end

		plot(T2,y1,'linewidth',2,'Color',[32, 161, 73]/255)
		title(['$' M_.endo_names_tex{idx} '$ ' M_.endo_names_long{idx}],'Interpreter','latex')
		xlim([2023 2100])
		hold on;
		plot(T2(1:length(y2)),y2,'--','linewidth',2,'Color',[224, 182, 11]/255)
		plot(T2(1:length(y3)),y3,':','linewidth',2,'Color',[33, 114, 176]/255)
		hold off;
		grid on;
		if i1==1
			ylabel('Paris-Agreement')
		elseif i1==nn+1
			ylabel('Laissez-faire')
		end
end


legend('Baseline rule','Natural rule','Steady-state rule')




