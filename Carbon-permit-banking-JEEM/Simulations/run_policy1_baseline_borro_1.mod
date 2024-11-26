//*************************************************************
//
// A General Equilibrium Approach to Carbon Permit Banking
//
// L. Dubois, J.-G. Sahuc, G. Vermandel  - August 2024
//
// Policy Exercise - Baseline Scenario
//
//*************************************************************


@#include "core_model.mod"


%-------------------------------
% 1. Simulations
%-------------------------------  

% simulation options
@#define t_simul = 1525
idex = strmatch('eta_e',M_.exo_names,'exact');
t0   = 1;
options_.ep.verbosity=1;

steady;
        %*******************************************************
        % 4.1 Previous situation:
        %   -Start in January 2023 (940 Million allowances)
        %   -cap decrease by 2.2% each year until the end  
        %*******************************************************


        initval;
            eta_e = 0;
        end;
        [y0, M_.params] = eval([ M_.fname  '.steadystate(oo_.steady_state, zeros(M_.exo_nbr,1), M_.params);']);


		@#if  variable_phi1 
			y0(strmatch('deltaphi',M_.endo_names,'exact')) = 1;
		@#endif 


        endval;
            eta_e = -0.970241988227596;
        end;
        yT = eval([ M_.fname  '.steadystate(oo_.steady_state, oo_.exo_steady_state, M_.params);']);

        perfect_foresight_setup(periods=@{t_simul});

		disp('DONE 1')
		
		if exist('run_policy1_baseline.mat') > 0
			load('run_policy1_baseline.mat')
			tsim0 = size(v1.endo_simul,2);
			for ix=1:v1_M_.endo_nbr
				id1 = strmatch(v1_M_.endo_names{ix},M_.endo_names,'exact');
				if ~isempty(id1)
					oo_.endo_simul(id1,2:(tsim0-1)) = v1.endo_simul(ix,2:(tsim0-1));
				end
			end
		end

		
%        oo_.exo_simul(:,idex) = linspace(0, -0.819660102439086,length(oo_.exo_simul(:,idex)));

        j=1/12;
        % Baisse de 2.2% par an
        for i=1:414
            oo_.exo_simul(t0+i,idex)  = -(43/1529)*j;
            j=j+1/12;
        end
        oo_.exo_simul(t0+415:end,idex) = -0.970241988227596;

        % set initial and terminal states
        oo_.endo_simul(:,1) = y0;
        oo_.endo_simul(strmatch('b',M_.endo_names,'exact'),1) = 8.5;  % (1087*(12/1529), pour démarrer à 1135)
        oo_.endo_simul(:,end) = yT;


        perfect_foresight_solver(lmmcp, maxit=5);
        v1 = oo_;

		disp('DONE 2')
       %*******************************************************
        % 4.2 Baseline:
        %   -start in december 2021 (940 millions of allowances)
        %   -LRF at 2.2% in 2023  
        %   -LRF at 4.3% in 2024-27
        %   -one shot reduction of 90Mt in 2024
        %   -one shot reduction of 27Mt in 2026
        %   -LRF at 4.4% from 2028 on
        %*******************************************************

        t1 = t0+5 ; % Annoucement of the new policy in May 2023

     	oo_.exo_simul = zeros(size(v1.exo_simul(t1:end,:)));

        endval;
            eta_e = -0.970241988227596;
        end;
        yT = eval([ M_.fname  '.steadystate(oo_.steady_state, oo_.exo_steady_state, M_.params);']);
      

        j=5/12;
        % Baisse de 2.2% jusqu'à la fin 2023
        for i=1:8
            oo_.exo_simul(i,idex) = -(43/1529)*j;
            j=j+1/12;
        end
     
        j = 1/12;
        % Reduction unique de 90Mt en 2024 et baisse de 4.3% par an en 2024 et 2025 
        for i=9:32 
            oo_.exo_simul(i,idex) = -90/1529 -43/1529 -(84/1529)*j;
            j=j+1/12;
        end
        % Reduction unique de 27Mt en 2026 et baisse de 4.3% par an en 2026 et 2027
        for i=33:56 
            oo_.exo_simul(i,idex) = -90/1529 -43/1529 -27/1529 -(84/1529)*j ;
            j=j+1/12;
        end
        j=1/12;
        % Baisse de 4.4% par an a partir de 2028
        for i=57:193
            oo_.exo_simul(i,idex) = -90/1529 -43/1529 -27/1529 -(84/1529)*4 -(86/1529)*j;
            j=j+1/12;
        end

        oo_.exo_simul(194:end,idex) = -0.970241988227596;

        oo_.endo_simul				= v1.endo_simul(:,t1:end);
        oo_.endo_simul(:,end) 		= yT;
        options_.periods			= size(oo_.endo_simul,2)-2;
        perfect_foresight_solver;
        v2 = oo_;
        v2.endo_simul = [v1.endo_simul(:,1:(t1-1)) v2.endo_simul(:,:)] ;

save run_policy1_baseline_borro_1 v2 M_;

