//*************************************************************
//
// A General Equilibrium Approach to Carbon Permits Banking
//
// L. Dubois, J.-G. Sahuc, G. Vermandel  - August 2024
//
// The Market Stability Reserve
//
//*************************************************************


@#define  msr = 1

@#include "core_model.mod"




%-------------------------------
% 4. Simulations
%-------------------------------  


% simulation options
@#define t_simul = 1525
@#define t_simul1 = 401
idex = strmatch('eta_e',M_.exo_names,'exact');
t0   = 1;



        %*******************************************************
        % 4.1 Previous situation:
        %   -Start in January 2023 (Cap at 1529 million allowances)
        %   -LRF at 2.2%  
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

		
		if exist('run_policy3_msr.mat') > 0
			load('run_policy3_msr.mat')
			tsim0 = size(v3.endo_simul,2);
			for ix=1:v3_M_.endo_nbr
				id3 = strmatch(v3_M_.endo_names{ix},M_.endo_names,'exact');
				if ~isempty(id3)
					oo_.endo_simul(id3,2:(tsim0-1)) = v3.endo_simul(ix,2:(tsim0-1));
				end
			end
		end

        oo_.exo_simul(1:t0,idex) = 0;

        j=1/12;
        % Baisse de 2.2% par an
        for i=1:416
            oo_.exo_simul(t0+i,idex)  = -(43/1529)*j;
            j=j+1/12;
        end

        oo_.exo_simul(t0+415:end,idex) = -0.970241988227596;

	
        % set initial and terminal states
        oo_.endo_simul(:,1) = y0;
        %oo_.endo_simul(strmatch('msr',M_.endo_names,'exact'),1) = 4.8;
        oo_.endo_simul(strmatch('b',M_.endo_names,'exact'),1) = 8.5;  % (1087*(12/1529), pour démarrer à 1135)
        oo_.endo_simul(:,end) = yT;

        perfect_foresight_solver(lmmcp, maxit=5);
        v3 = oo_;



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

     	oo_.exo_simul = zeros(size(v3.exo_simul(t1:end,:)));
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
        for i=57:193
            oo_.exo_simul(i,idex) = -90/1529 -43/1529 -27/1529 -(84/1529)*4 -(86/1529)*j;
            j=j+1/12;
        end


        oo_.exo_simul(194:end,idex) = -0.970241988227596;


        oo_.endo_simul				= v3.endo_simul(:,t1:end);
        oo_.endo_simul(:,end) 		= yT;
        options_.periods			= size(oo_.endo_simul,2)-2;
        perfect_foresight_solver;
        v4 = oo_;
        v4.endo_simul = [v3.endo_simul(:,1:(t1-1)) v4.endo_simul(:,:)] ;


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

     	oo_.exo_simul = zeros(size(v3.exo_simul(t1:end,:)));
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
        % Reduction unique de 90Mt en 2024 et baisse de 10% par an en 2024 et 2025 
        for i=9:32 
            oo_.exo_simul(i,idex) = -90/1529 -43/1529 -(195.4/1529)*j;
            j=j+1/12;
        end
        % Reduction unique de 27Mt en 2026 et baisse de 10% par an
        for i=33:89 
            oo_.exo_simul(i,idex) = -90/1529 -43/1529 -27/1529 -(195.4/1529)*j ;
            j=j+1/12;
        end


        oo_.exo_simul(90:end,idex) = -0.970241988227596;


        oo_.endo_simul				= v3.endo_simul(:,t1:end);
        oo_.endo_simul(:,end) 		= yT;
        options_.periods			= size(oo_.endo_simul,2)-2;
        perfect_foresight_solver;
        v5 = oo_;
        v5.endo_simul = [v3.endo_simul(:,1:(t1-1)) v5.endo_simul(:,:)] ;




%-------------------------------
% 5. Figure
%-------------------------------
load('run_policy1_baseline.mat')

vars = {'thet' 'b' 'lambda' 'e' 'pe' 'mu_cost' 'y' 'c' 'i'};
var_names = {'Permit supply' 'Permit bank' 'Shadow permit value' 'Carbon emissions'...
 'Carbon price' 'Abatement' 'Output' 'Consumption' 'Investment'} ;
date = [2022+11/12:1/12:2060] ;
thresh = 1081-600+240-180-60; 

y_v1 = v1.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley1 = y_v1./y_v1(:,1); 
yy_v1 = 100*(scaley1-1);
y_v2 = v2.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley2 = y_v2./y_v2(:,1); 
yy_v2 = 100*(scaley2-1);
y_v3 = v3.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley3 = y_v3./y_v2(:,1); 
yy_v3 = 100*(scaley3-1);
y_v4 = v4.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley4 = y_v4./y_v2(:,1); 
yy_v4 = 100*(scaley4-1);
y_v5 = v5.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley5 = y_v5./y_v2(:,1); 
yy_v5 = 100*(scaley5-1);

c_v1 = v1.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v2 = v2.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v3 = v3.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v4 = v4.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v5 = v5.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
scalec1 =  c_v1./ c_v1(:,1); scalec2 =  c_v2./ c_v2(:,1); scalec3 =  c_v3./ c_v2(:,1); scalec4 =  c_v4./ c_v2(:,1); scalec5 =  c_v5./ c_v2(:,1); 
cc_v1 = 100*(scalec1-1); cc_v2 = 100*(scalec2-1); cc_v3 = 100*(scalec3-1); cc_v4 = 100*(scalec4-1); cc_v5 = 100*(scalec5-1);


x_v1 = v1.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v2 = v2.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v3 = v3.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v4 = v4.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v5 = v5.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
scalex1 =  x_v1./ x_v1(:,1); scalex2 =  x_v2./ x_v2(:,1); scalex3 =  x_v3./ x_v2(:,1); scalex4 =  x_v4./ x_v2(:,1); scalex5 =  x_v5./ x_v2(:,1); 
xx_v1 = 100*(scalex1-1); xx_v2 = 100*(scalex2-1); xx_v3 = 100*(scalex3-1); xx_v4 = 100*(scalex4-1); xx_v5 = 100*(scalex5-1);


coef = (phi_1*(MU)^phi_2)*100;
muc_v1 = v1.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v1 = (coef./muc_v1(1))*muc_v1 ;
muc_v2 = v2.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v2 = (coef./muc_v2(1))*muc_v2 ;
muc_v3 = v3.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v3 = (coef./v1.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1))*muc_v3 ;
muc_v4 = v4.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v4 = (coef./v2.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1))*muc_v4 ;
muc_v5 = v5.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v5 = (coef./v2.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1))*muc_v5 ;


p_s1 = (1529/12)*v1.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;
p_s2 = (1529/12)*v2.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;
p_s3 = (1529/12)*v3.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;
p_s4 = (1529/12)*v4.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;
p_s5 = (1529/12)*v5.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;


b1 = (1529/12)*v1.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;
b2 = (1529/12)*v2.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;
b3 = (1529/12)*v3.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;
b4 = (1529/12)*v4.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;
b5 = (1529/12)*v5.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;


e1 = (1529/12)*v1.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;
e2 = (1529/12)*v2.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;
e3 = (1529/12)*v3.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;
e4 = (1529/12)*v4.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;
e5 = (1529/12)*v5.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;



pe1 = v1.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe1 = (80./pe1(1))*pe1 ;
pe2 = v2.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe2 = (80./pe2(1))*pe2 ;
pe3 = v3.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe3 = (80./v1.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1))*pe3 ;
pe4 = v4.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe4 = (80./v2.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1))*pe4 ;
pe5 = v5.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe5 = (80./v2.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1))*pe5 ;


b2 = [b2(1),b2];
e2 = [e2(1),e2];
pe2 = [pe2(1),pe2];
muc_v2 = [muc_v2(1),muc_v2];
yy_v2 = [yy_v2(1),yy_v2];
cc_v2 = [cc_v2(1),cc_v2];
xx_v2 = [xx_v2(1),xx_v2];
p_s4 = [p_s2(1),p_s4];
b4 = [b2(1),b4];
e4 = [e2(1),e4];
pe4 = [pe2(1),pe4];
muc_v4 = [muc_v2(1),muc_v4];
yy_v4 = [yy_v2(1),yy_v4];
cc_v4 = [cc_v2(1),cc_v4];
xx_v4 = [xx_v2(1),xx_v4];
p_s5 = [p_s2(1),p_s5];
b5 = [b2(1),b5];
e5 = [e2(1),e5];
pe5 = [pe2(1),pe5];
muc_v5 = [muc_v2(1),muc_v5];
yy_v5 = [yy_v2(1),yy_v5];
cc_v5 = [cc_v2(1),cc_v5];
xx_v5 = [xx_v2(1),xx_v5];

bbar = 833*ones(size(date));
btilde = 1096*ones(size(date));



figure
subplot(3,3,1)
  hold on
  plot(date,p_s2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,p_s4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0.4, 0.69, 0.2])
  hold off 
  title({'Permit supply','(in million permits)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,2)
  hold on
  plot(date,b2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,b4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0.4, 0.69, 0.2])
  plot(date,bbar(1:length(date)),'linewidth',1,'LineStyle',':','Color',[0, 0, 0])
  plot(date,btilde(1:length(date)),'linewidth',1,'LineStyle',':','Color',[0, 0, 0])
  hold off
  title({'Permit bank', '(in million permits)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,3)
  hold on
  plot(date,pe2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,pe4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0.4, 0.69, 0.2])
  hold off
  title({'Carbon price', '(in euros)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,4)
  hold on
  plot(date,e2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,e4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0.4, 0.69, 0.2])
  ylim([3 80])
  hold off
  title({'Carbon emissions', '(in million tons of CO2)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,5)
  hold on
  plot(date,muc_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,muc_v4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0.4, 0.69, 0.2])
  hold off
  title({'Abatement costs','(in percentage of output)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal'; ax.YAxis.Exponent = 0;
  grid on

subplot(3,3,6)
  hold on
  plot(date, yy_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date, yy_v4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0.4, 0.69, 0.2])
  hold off
  title({'Output', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,7)
  hold on
  plot(date,cc_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date,cc_v4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0.4, 0.69, 0.2])
  hold off
  title({'Consumption', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,8)
  hold on
  plot(date, xx_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  plot(date, xx_v4(1:length(date)),'linewidth',1.5,'LineStyle','-.','Color',[0.4, 0.69, 0.2])
  hold off
  title({'Investment', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on
%legend('Baseline','Baseline with MSR','Orientation','horizontal','Position',[0.42 0.011 0.1804 0.0710],'box','off','FontSize',10)

legend({'Baseline','Baseline with MSR'},'Position',[0.75 0.15 0.1 0.1], 'Interpreter','Latex', 'FontSize', 10)
legend boxoff

set(gcf,'PaperType','<custom>','PaperSize',[20.4 15.25],'PaperPosition',[0.1 0.02158 20.3046 15.2284]);
print -dpdf Figure_policy3_msr.pdf

v3_M_ = M_;
save run_policy3_msr v3 v4 v5 v3_M_;
