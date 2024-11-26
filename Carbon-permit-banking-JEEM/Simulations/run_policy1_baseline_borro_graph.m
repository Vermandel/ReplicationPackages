%*************************************************************
%
% A General Equilibrium Approach to Carbon Permit Banking
%
% L. Dubois, J.-G. Sahuc, G. Vermandel  - August 2024
%
% Policy Exercise - Baseline Scenario with Borrowing
%
%*************************************************************


load run_policy1_baseline_borro_1;
load run_policy1_baseline_borro_2;


%-------------------------------
% 2. Figure
%-------------------------------

vars = {'thet' 'b' 'lambda' 'e' 'pe' 'mu_cost' 'y' 'c' 'i' 'lambda'};
var_names = {'Permit supply' 'Permit bank' 'Shadow permit value' 'Carbon emissions'...
 'Carbon price' 'Abatement' 'Output' 'Consumption' 'Investment'} ;
date = [2023:1/12:2060] ;
thresh = 1081-600+240-180-60-2+2; 
t0   = 1;


y_v2 = v2.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley2 = y_v2./y_v2(:,1); 
yy_v2 = 100*(scaley2-1);
y_v22 = v22.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley22 = y_v22./y_v22(:,1); 
yy_v22 = 100*(scaley22-1);

c_v22 = v22.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v2 = v2.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
scalec22 =  c_v22./ c_v22(:,1); scalec2 =  c_v2./ c_v2(:,1); 
cc_v22 = 100*(scalec22-1); cc_v2 = 100*(scalec2-1);

x_v22 = v22.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v2 = v2.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
scalex22 =  x_v22./ x_v22(:,1); scalex2 =  x_v2./ x_v2(:,1); 
xx_v22 = 100*(scalex22-1); xx_v2 = 100*(scalex2-1);


cnorm = (phi_1*(MU)^phi_2)*100;
muc_v22 = v22.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v22 = (cnorm./muc_v22(1))*muc_v22 ;
muc_v2 = v2.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v2 = (cnorm./muc_v2(1))*muc_v2 ;

p_s22 = (1529/12)*v22.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;
p_s2 = (1529/12)*v2.endo_simul(strmatch(vars{1},M_.endo_names,'exact'),t0+1:end-thresh) ;

b22 = (1529/12)*v22.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;
b2 = (1529/12)*v2.endo_simul(strmatch(vars{2},M_.endo_names,'exact'),t0+1:end-thresh) ;

e22 = (1529/12)*v22.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;
e2 = (1529/12)*v2.endo_simul(strmatch(vars{4},M_.endo_names,'exact'),t0+1:end-thresh) ;

pe22 = v22.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe22 = (80./pe22(1))*pe22 ;
pe2 = v2.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe2 = (80./pe2(1))*pe2 ;




figure
subplot(3,3,1)
  plot(date,p_s2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold on
    plot(date,p_s22(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.57, 0.64, 0.69])
  hold off 
  title({'Permit supply','(in million permits)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,2)
  plot(date,b2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold on
  plot(date,b22(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.57, 0.64, 0.69])
  hold off
  title({'Permit bank', '(in million permits)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,3)
  plot(date,pe2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold on
    plot(date,pe22(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.57, 0.64, 0.69])
  hold off
  title({'Carbon price', '(in euros)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,4)
  plot(date,e2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold on
  plot(date,e22(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.57, 0.64, 0.69])
  hold off
  ylim([0 80])
  title({'Carbon emissions', '(in million tons of CO2)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on
  ylim([0 100])

subplot(3,3,5)
  plot(date,muc_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold on 
  plot(date,muc_v22(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.57, 0.64, 0.69])
  hold off
  title({'Abatement costs','(in percentage of output)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal'; ax.YAxis.Exponent = 0;
  grid on
  ylim([0 2])

subplot(3,3,6)
  plot(date, yy_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold on
  plot(date, yy_v22(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.57, 0.64, 0.69])
  hold off
  title({'Output', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,7)
  plot(date,cc_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
 hold on 
  plot(date,cc_v22(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.57, 0.64, 0.69])
  hold off
  title({'Consumption', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,3,8)
  plot(date, xx_v2(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold on
  plot(date, xx_v22(1:length(date)),'linewidth',1.5,'LineStyle','--','Color',[0.57, 0.64, 0.69])
  hold off
  title({'Investment', '(in percentage deviation)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',9)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

  legend({'Baseline','Baseline with limited borrowing'},'Position',[0.75 0.15 0.1 0.1], 'Interpreter','Latex', FontSize=10)
  legend boxoff

set(gcf,'PaperType','<custom>','PaperSize',[20.4 15.25],'PaperPosition',[0.1 0.02158 20.3046 15.2284]);
print -dpdf Figure_policy1_baseline_borro.pdf;

close all;


  plot(date,e2(1:length(date)),'linewidth',2,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  hold on
  plot(date,e22(1:length(date)),'linewidth',2,'LineStyle','--','Color',[0.57, 0.64, 0.69])
  hold off
  ylim([0 80])
  %title({'Carbon emissions', '(in million tons of CO2)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',14)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on
  ylim([0 100])
  legend({'Baseline','Baseline with limited borrowing'}, 'Interpreter','Latex', FontSize=14)

set(gcf,'PaperType','<custom>','PaperSize',[20.4 15.25],'PaperPosition',[0.1 0.02158 20.3046 15.2284]);
print -dpdf Figure_policy1_baseline_borro2.pdf;

%'Position',[0.75 0.15 0.1 0.1],