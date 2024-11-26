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




%-------------------------------
% 2. Figure
%-------------------------------

vars = {'thet' 'b' 'lambda' 'e' 'pe' 'mu_cost' 'y' 'c' 'i'};
var_names = {'Permit supply' 'Permit bank' 'Shadow permit value' 'Carbon emissions'...
 'Carbon price' 'Abatement' 'Output' 'Consumption' 'Investment'} ;
date = [2023:1/12:2060] ;
thresh = 1081-600+240-180-60-2+2; 

y_v1 = v1.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley1 = y_v1./y_v1(:,1); 
yy_v1 = 100*(scaley1-1);
y_v2 = v2.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley2 = y_v2./y_v2(:,1); 
yy_v2 = 100*(scaley2-1);

c_v1 = v1.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v2 = v2.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
scalec1 =  c_v1./ c_v1(:,1); scalec2 =  c_v2./ c_v2(:,1); 
cc_v1 = 100*(scalec1-1); cc_v2 = 100*(scalec2-1);

x_v1 = v1.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v2 = v2.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
scalex1 =  x_v1./ x_v1(:,1); scalex2 =  x_v2./ x_v2(:,1); 
xx_v1 = 100*(scalex1-1); xx_v2 = 100*(scalex2-1);


cnorm = (phi_1*(MU)^phi_2)*100;
muc_v1 = v1.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v1 = (cnorm./muc_v1(1))*muc_v1 ;
muc_v2 = v2.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v2 = (cnorm./muc_v2(1))*muc_v2 ;


p_s1 = (1529/12)*v1.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;
p_s2 = (1529/12)*v2.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;

b1 = (1529/12)*v1.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;
b2 = (1529/12)*v2.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;

e1 = (1529/12)*v1.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;
e2 = (1529/12)*v2.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;

pe1 = v1.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe1 = (80./pe1(1))*pe1 ;
pe2 = v2.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe2 = (80./pe2(1))*pe2 ;




figure
subplot(3,3,1)
  plot(date,p_s1(1:length(date)),'linewidth',1.5,'LineStyle',':','Color',[0.57, 0.64, 0.69])
  hold on
  plot(date,p_s2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold off 
  title({'Permit supply','(in million permits)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,2)
  plot(date,b1(1:length(date)),'linewidth',1.5,'LineStyle',':','Color',[0.57, 0.64, 0.69])
  hold on
  plot(date,b2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold off
  title({'Permit bank', '(in million permits)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,3)
  plot(date,pe1(1:length(date)),'linewidth',1.5,'LineStyle',':','Color',[0.57, 0.64, 0.69])
  hold on
  plot(date,pe2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold off
  title({'Carbon price', '(in euros)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,4)
  plot(date,e1(1:length(date)),'linewidth',1.5,'LineStyle',':','Color',[0.57, 0.64, 0.69])
  hold on
  plot(date,e2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold off
  ylim([0 80])
  title({'Carbon emissions', '(in million tons of CO2)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on
  ylim([0 100])

subplot(3,3,5)
  plot(date,muc_v1(1:length(date)),'linewidth',1.5,'LineStyle',':','Color',[0.57, 0.64, 0.69])
  hold on
  plot(date,muc_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold off
  title({'Abatement costs','(in percentage of output)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal'; ax.YAxis.Exponent = 0;
  grid on
  ylim([0 2])

subplot(3,3,6)
  plot(date, yy_v1(1:length(date)),'linewidth',1.5,'LineStyle',':','Color',[0.57, 0.64, 0.69])
  hold on
  plot(date, yy_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold off
  title({'Output', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,7)
  plot(date,cc_v1(1:length(date)),'linewidth',1.5,'LineStyle',':','Color',[0.57, 0.64, 0.69])
  hold on
  plot(date,cc_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold off
  title({'Consumption', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,8)
  plot(date, xx_v1(1:length(date)),'linewidth',1.5,'LineStyle',':','Color',[0.57, 0.64, 0.69])
  hold on
  plot(date, xx_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold off
  title({'Investment', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

set(gcf,'PaperType','<custom>','PaperSize',[20.4 15.25],'PaperPosition',[0.1 0.02158 20.3046 15.2284]);
print -dpdf Figure_policy1_baseline.pdf

v1_M_ = M_;
save run_policy1_baseline v1 v2 v1_M_;
