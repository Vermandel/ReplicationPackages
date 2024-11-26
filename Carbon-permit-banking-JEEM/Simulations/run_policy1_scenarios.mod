//*************************************************************
//
// A General Equilibrium Approach to Carbon Permits Banking
//
// L. Dubois, J.-G. Sahuc, G. Vermandel  - August 2024
//
// Cap scenarios from 2031
//
//*************************************************************


@#include "core_model.mod"


%-------------------------------
% 1. Simulations
%-------------------------------  

% simulation options
@#define t_simul = 1525

@#define t_simul1 = 401
idex = strmatch('eta_e',M_.exo_names,'exact');
t0   = 1;

steady;

        %*******************************************************
        % 4.1 Previous & Baseline:
        % Load from run_policy1_baseline.mod
        %*******************************************************

        load('run_policy1_baseline.mat')
 

        %*******************************************************
        % 4.3 Baseline with return to 2.2% in Phase V:
        %   -start in december 2021 (940 millions of allowances)
        %   -LRF at 2.2% in 2023  
        %   -LRF at 4.3% in 2024-27
        %   -one shot reduction of 90Mt in 2024
        %   -one shot reduction of 27Mt in 2026
        %   -LRF at 4.4% in 2028-30
        %   -LRF at 2.2% from 2031 on
        %*******************************************************

        t1 = t0+5 ; % Annoucement of the new policy in May 2023

     	oo_.exo_simul = zeros(size(v1.exo_simul(t1:end,:)));
        %perfect_foresight_setup(periods=@{t_simul1});

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
        for i=57:92
            oo_.exo_simul(i,idex) = -90/1529 -43/1529 -27/1529 -(84/1529)*4 -(86/1529)*j;
            j=j+1/12;
        end
        j=1/12;
        % Baisse de 2.2% par an a partir de 2031
        for i=93:295
            oo_.exo_simul(i,idex) = -90/1529 -43/1529 -27/1529 -(84/1529)*4 -(86/1529)*3 -(43/1529)*j;
            j=j+1/12;
        end

        oo_.exo_simul(296:end,idex) = -0.970241988227596;


        oo_.endo_simul				= v1.endo_simul(:,t1:end);
        oo_.endo_simul(:,end) 		= yT;
        options_.periods			= size(oo_.endo_simul,2)-2;
        perfect_foresight_solver(lmmcp, maxit=5);
        v3 = oo_;
        v3.endo_simul = [v1.endo_simul(:,1:(t1-1)) v3.endo_simul(:,:)] ;



        %*******************************************************
        % 4.3 Baseline with increase to 10% in Phase V:
        %   -start in december 2021 (940 millions of allowances)
        %   -LRF at 2.2% in 2023  
        %   -LRF at 4.3% in 2024-27
        %   -one shot reduction of 90Mt in 2024
        %   -one shot reduction of 27Mt in 2026
        %   -LRF at 4.4% in 2028-30
        %   -LRF at 10% from 2031 on
        %*******************************************************

        t1 = t0+5 ; % Annoucement of the new policy in May 2023

     	oo_.exo_simul = zeros(size(v1.exo_simul(t1:end,:)));
        %perfect_foresight_setup(periods=@{t_simul1});

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
        for i=57:92
            oo_.exo_simul(i,idex) = -90/1529 -43/1529 -27/1529 -(84/1529)*4 -(86/1529)*j;
            j=j+1/12;
        end
        j=1/12;
        % Baisse de 10% par an a partir de 2031
        for i=93:136
            oo_.exo_simul(i,idex) = -90/1529 -43/1529 -27/1529 -(84/1529)*4 -(86/1529)*3 -(195.4/1529)*j;
            j=j+1/12;
        end


        oo_.exo_simul(137:end,idex) = -0.970241988227596;


        oo_.endo_simul				= v1.endo_simul(:,t1:end);
        oo_.endo_simul(:,end) 		= yT;
        options_.periods			= size(oo_.endo_simul,2)-2;
        perfect_foresight_solver(lmmcp, maxit=5);
        v4 = oo_;
        v4.endo_simul = [v1.endo_simul(:,1:(t1-1)) v4.endo_simul(:,:)] ;


%-------------------------------
% 5. Figure
%-------------------------------

vars = {'thet' 'b' 'lambda' 'e' 'pe' 'mu_cost' 'y' 'c' 'i'};
var_names = {'Permit supply' 'Permit bank' 'Shadow permit value' 'Carbon emissions'...
 'Carbon price' 'Abatement' 'Output' 'Consumption' 'Investment'} ;
date = [2023:1/12:2060] ;
thresh = 1081-600+240-180-60; 



y_v1 = v1.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley1 = y_v1./y_v1(:,1); 
yy_v1 = 100*(scaley1-1);
y_v2 = v2.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley2 = y_v2./y_v2(:,1); 
yy_v2 = 100*(scaley2-1);
y_v3 = v3.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley3 = y_v3./y_v3(:,1); 
yy_v3 = 100*(scaley3-1);
y_v4 = v4.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley4 = y_v4./y_v4(:,1); 
yy_v4 = 100*(scaley4-1);

c_v1 = v1.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v2 = v2.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v3 = v3.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v4 = v4.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
scalec1 =  c_v1./ c_v1(:,1); scalec2 =  c_v2./ c_v2(:,1); scalec3 =  c_v3./ c_v3(:,1); scalec4 =  c_v4./ c_v4(:,1); 
cc_v1 = 100*(scalec1-1); cc_v2 = 100*(scalec2-1); cc_v3 = 100*(scalec3-1); cc_v4 = 100*(scalec4-1);

x_v1 = v1.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v2 = v2.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v3 = v3.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v4 = v4.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
scalex1 =  x_v1./ x_v1(:,1); scalex2 =  x_v2./ x_v2(:,1); scalex3 =  x_v3./ x_v3(:,1); scalex4 =  x_v4./ x_v4(:,1); 
xx_v1 = 100*(scalex1-1); xx_v2 = 100*(scalex2-1); xx_v3 = 100*(scalex3-1); xx_v4 = 100*(scalex4-1);

coef = (phi_1*(MU)^phi_2)*100;
muc_v1 = v1.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v1 = (coef./muc_v1(1))*muc_v1 ;
muc_v2 = v2.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v2 = (coef./muc_v2(1))*muc_v2 ;
muc_v3 = v3.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v3 = (coef./muc_v3(1))*muc_v3 ;
muc_v4 = v4.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v4 = (coef./muc_v4(1))*muc_v4 ;


p_s1 = (1529/12)*v1.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;
p_s2 = (1529/12)*v2.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;
p_s3 = (1529/12)*v3.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;
p_s4 = (1529/12)*v4.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;

b1 = (1529/12)*v1.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;
b2 = (1529/12)*v2.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;
b3 = (1529/12)*v3.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;
b4 = (1529/12)*v4.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;

e1 = (1529/12)*v1.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;
e2 = (1529/12)*v2.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;
e3 = (1529/12)*v3.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;
e4 = (1529/12)*v4.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;


pe1 = v1.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe1 = (80./pe1(1))*pe1 ;
pe2 = v2.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe2 = (80./pe2(1))*pe2 ;
pe3 = v3.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe3 = (80./pe3(1))*pe3 ;
pe4 = v4.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe4 = (80./pe4(1))*pe4 ;


figure
subplot(3,3,1)
  hold on
  plot(date,p_s2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,p_s3(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.4, 0.69, 0.2])
  plot(date,p_s4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0., 0.42, 0.24])
  hold off 
  title({'Permit supply','(in million permits)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,2)
  hold on
  plot(date,b2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,b3(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.4, 0.69, 0.2])
  plot(date,b4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0., 0.42, 0.24])
  hold off
  title({'Permit bank', '(in million permits)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,3)
  hold on
  plot(date,pe2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,pe3(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.4, 0.69, 0.2])
  plot(date,pe4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0., 0.42, 0.24])
  hold off
  title({'Carbon price', '(in euros)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,4)
  hold on
  plot(date,e2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,e3(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.4, 0.69, 0.2])
  plot(date,e4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0., 0.42, 0.24])
  hold off
  title({'Carbon emissions', '(in million tons of CO2)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,5)
  hold on
  plot(date,muc_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,muc_v3(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.4, 0.69, 0.2])
  plot(date,muc_v4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0., 0.42, 0.24])
  hold off
  title({'Abatement costs','(in percentage of output)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal'; ax.YAxis.Exponent = 0;
  grid on


subplot(3,3,6)
  hold on
  plot(date, yy_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date, yy_v3(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.4, 0.69, 0.2])
  plot(date, yy_v4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0., 0.42, 0.24])
  hold off
  title({'Output', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,7)
  hold on
  plot(date,cc_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,cc_v3(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.4, 0.69, 0.2])
  plot(date,cc_v4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0., 0.42, 0.24])
  hold off
  title({'Consumption', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,8)
  hold on
  plot(date, xx_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date, xx_v3(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.4, 0.69, 0.2])
  plot(date, xx_v4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0., 0.42, 0.24])
  hold off
  title({'Investment', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on
%legend('Baseline','Alternative I','Alternative II','Orientation','horizontal','Position',[0.42 0.011 0.1804 0.0710],'box','off','FontSize',10)

legend({'Baseline','Alternative I','Alternative II'},'Position',[0.75 0.15 0.1 0.1], 'Interpreter','Latex', FontSize=10)
legend boxoff

set(gcf,'PaperType','<custom>','PaperSize',[20.4 15.25],'PaperPosition',[0.1 0.02158 20.3046 15.2284]);
print -dpdf Figure_policy1_scenarios.pdf

save run_policy1_scenarios
