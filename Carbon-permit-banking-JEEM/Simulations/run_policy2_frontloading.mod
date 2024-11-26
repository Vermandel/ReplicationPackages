//*************************************************************
//
// A General Equilibrium Approach to Carbon Permit Banking
//
// L. Dubois, J.-G. Sahuc, G. Vermandel  - August 2024
//
// Frontloading policy
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

        %*******************************************************
        % 4.1 Previous & Baseline
        % Retrieved from run_policy1_baseline.mod
        %*******************************************************

        load('run_policy1_baseline.mat');
        
        
        %*******************************************************
        % 4.2 Frontloading: 
        % (first step, before the announcements of may 2023)
        %*******************************************************
        
        
        initval;
            eta_e = 0;
        end;
        [y0, M_.params] = eval([ M_.fname  '.steadystate(oo_.steady_state, zeros(M_.exo_nbr,1), M_.params);']);


        endval;
            eta_e = -0.970241988227596;
        end;
        yT = eval([ M_.fname  '.steadystate(oo_.steady_state, oo_.exo_steady_state, M_.params);']);

        perfect_foresight_setup(periods=@{t_simul});
        
        
        oo_.exo_simul = v1.exo_simul ;


        oo_.exo_simul(2:49,idex) = oo_.exo_simul(2:49,idex) + 20/1529;
        oo_.exo_simul(50:97,idex) = oo_.exo_simul(50:97,idex) - 20/1529;
        oo_.endo_simul = v1.endo_simul; % Just to speed things up

        perfect_foresight_solver(lmmcp, maxit=200);
     	v3 = oo_;


    	%*******************************************************
        % 4.3 Frontloading: 
        % (second step, after the announcements of may 2023)
        %*******************************************************

        t1 = t0+5 ; % Annoucement of the new policy in May 2023

     	oo_.exo_simul = v2.exo_simul;

        oo_.exo_simul(1:44,idex) = oo_.exo_simul(1:44,idex) + 20/1529;
        oo_.exo_simul(45:92,idex) = oo_.exo_simul(45:92,idex) - 20/1529;


        oo_.endo_simul = v2.endo_simul(:,6:end); % Just to speed things up
        oo_.endo_simul				= v3.endo_simul(:,t1:end);
        oo_.endo_simul(:,end) 		= v2.endo_simul(:,end);

        options_.periods			= size(oo_.endo_simul,2)-2;
        perfect_foresight_solver;
        v4 = oo_;
        v4.endo_simul = [v3.endo_simul(:,1:(t1-1)) v4.endo_simul(:,:)] ;
        


%-------------------------------
% 5. Figure
%-------------------------------

vars = {'thet' 'b' 'lambda' 'e' 'pe' 'mu_cost' 'y' 'c' 'i'};
var_names = {'Permit supply' 'Permit bank' 'Shadow permit value' 'Carbon emissions'...
 'Carbon price' 'Abatement' 'Output' 'Consumption' 'Investment'} ;
date = [2022+11/12:1/12:2032] ;
thresh = 1081-600+240-180-60; 


sd1= 2; ed1=49 ;sd2=50; ed2 =97 ;


y_v1 = v1.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley = y_v1./y_v1(:,1); 
yy_v1 = scaley*1080;          %billion 
y_v2 = v2.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley = y_v2./y_v2(:,1); 
yy_v2 = scaley*1080;
y_v3 = v3.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley = y_v3./v1.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1); 
yy_v3 = scaley*1080;
y_v4 = v4.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1:end-thresh);
scaley = y_v4./v2.endo_simul(strmatch(vars{7},M_.endo_names,'exact'),t0+1); 
yy_v4 = scaley*1080;

c_v1 = v1.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v2 = v2.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v3 = v3.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
c_v4 = v4.endo_simul(strmatch(vars{8},M_.endo_names,'exact'),t0+1:end-thresh);
scalec1 =  c_v1./ y_v1; scalec2 =  c_v2./ y_v2; scalec3 =  c_v3./ y_v3; scalec4 =  c_v4./ y_v4;
cc_v1 = scalec1.*yy_v1; cc_v2 = scalec2.*yy_v2; cc_v3 = scalec3.*yy_v3; cc_v4 = scalec4.*yy_v4;


x_v1 = v1.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v2 = v2.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v3 = v3.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
x_v4 = v4.endo_simul(strmatch(vars{9},M_.endo_names,'exact'),t0+1:end-thresh);
scalex1 =  x_v1./ y_v1; scalex2 =  x_v2./ y_v2; scalex3 =  x_v3./ y_v3; scalex4 =  x_v4./ y_v4;
xx_v1 = scalex1.*yy_v1; xx_v2 = scalex2.*yy_v2; xx_v3 = scalex3.*yy_v3; xx_v4 = scalex4.*yy_v4;


coef = (phi_1*(MU)^phi_2)*1080;
muc_v1 = v1.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v1 = (coef./muc_v1(1))*muc_v1 ;
muc_v2 = v2.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v2 = (coef./muc_v2(1))*muc_v2 ;
muc_v3 = v3.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v3 = (coef./v1.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1))*muc_v3 ;
muc_v4 = v4.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1:end-thresh) ;
muc_v4 = (coef./v2.endo_simul(strmatch(vars{6},M_.endo_names,'exact'),t0+1))*muc_v4 ;


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
pe3 = (80./v1.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1))*pe3 ;
pe4 = v4.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1:end-thresh) ;
pe4 = (80./v2.endo_simul(strmatch(vars{5},M_.endo_names,'exact'),t0+1))*pe4 ;

diff_permits = [0,p_s4(1:109)-p_s2(1:109)];
diff_bank = [0,b4(1:109)-b2(1:109)];
diff_emissions = e4-e2;
diff_price = pe4-pe2;
diff_abatement = muc_v4-muc_v2;
diff_y =  yy_v4 - yy_v2;
diff_c =  cc_v4 - cc_v2;
diff_x =  xx_v4 - xx_v2;


figure
subplot(3,2,1)
plot(date,diff_permits(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88])
  title({'Permit supply difference','(in million permits)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',10)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on

subplot(3,2,2)
plot(date,diff_bank(1:length(date)),'linewidth',1.5,'LineStyle','-','Color',[0.25, 0.41, 0.88]) 
title({'Permit bank difference','(in million permits)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',10)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
  grid on


subplot(3,2,3)
X = categorical({'2023-2026','2027-2030'});
X = reordercats(X,{'2023-2026','2027-2030'});
cum_e1=cumsum(diff_emissions(sd1:ed1));
cum_e2=cumsum(diff_emissions(sd2:ed2));
bar(X,[cum_e1(end),cum_e2(end)],'BarWidth',0.5,'FaceColor',[0.25, 0.41, 0.88])
title({'Cumulated carbon emission difference','(in million tons of CO2)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',10)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
grid on

subplot(3,2,4)
X = categorical({'2023-2026','2027-2030'});
X = reordercats(X,{'2023-2026','2027-2030'});
bar(X,[mean(diff_price(sd1:ed1)),mean(diff_price(sd2:ed2))],'BarWidth',0.5,'FaceColor',[0.25, 0.41, 0.88])
title({'Mean carbon price difference','(in euros)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',10)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
grid on


subplot(3,2,5)
X = categorical({'2023-2026','2027-2030'});
X = reordercats(X,{'2023-2026','2027-2030'});
cum_y1=cumsum(diff_y(sd1:ed1));
cum_y2=cumsum(diff_y(sd2:ed2));
bar(X,[cum_y1(end),cum_y2(end)],'BarWidth',0.5,'FaceColor',[0.25, 0.41, 0.88])
title({'Cumulated output difference', '(in billion euros)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',10)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
grid on

subplot(3,2,6)
X = categorical({'2023-2026','2027-2030'});
X = reordercats(X,{'2023-2026','2027-2030'});
cum_c1=cumsum(diff_c(sd1:ed1));
cum_c2=cumsum(diff_c(sd2:ed2));
bar(X,[cum_c1(end),cum_c2(end)],'BarWidth',0.5,'FaceColor',[0.25, 0.41, 0.88])
title({'Cumulated consumption difference', '(in billion euros)'},'Interpreter','Latex');
  set(gca,'fontname','times','FontSize',10)
  ax=gca; c=ax.FontWeight; ax.TitleFontWeight='normal';
grid on


set(gcf,'PaperType','<custom>','PaperSize',[20.4 15.25],'PaperPosition',[0.1 0.02158 20.3046 15.2284]);
print -dpdf Figure_policy2_frontloading.pdf  
