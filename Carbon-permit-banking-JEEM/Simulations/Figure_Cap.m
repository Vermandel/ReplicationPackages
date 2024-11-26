
% Figure Cap


[previous,baseline,alt1,alt2] = deal(nan(1,43)) ;

%% Construct Previous
previous(5) = 1529;
previous(4) = 1529 + 43;
previous(3) = previous(4) +240;
for i = 2:-1:1
    previous(i) = previous(i+1) + 34 ;
end


for i = 5:size(previous,2)
   previous(i) = previous(i-1) - 43 ;
end
previous = max(previous,0) ;

%% Construct Baseline
baseline(1:6) = previous(1:6) ;
%baseline(7) = previous(6) - 90 - 84 ;
baseline(7) = baseline(6) - 90 - 84 ;
baseline(8) = baseline(7) - 84 ;
baseline(9) = baseline(8) - 27 - 84 ;
baseline(10) = baseline(9) - 84 ;

for i = 11:size(baseline,2)
    baseline(i) = baseline(i-1) - 86 ;
end
baseline = max(baseline,0) ;

%% Construct Alt1

alt1(14) = baseline(14);

for i = 15:size(alt1,2)
    alt1(i) = alt1(i-1) - 43 ;
end
alt1(14:end) = max(alt1(14:end),0) ;

%% Construct Alt2

alt2(14) = baseline(14);

for i = 15:size(alt2,2)
    alt2(i) = alt2(i-1) - 105.5 ;
end
alt2(14:end) = max(alt2(14:end),0) ;

%% Graph

date = [2018:2060] ;

figure
plot(date,previous,'linewidth',2,'LineStyle','--','Color',[0.57, 0.64, 0.69])
hold on
plot(date,baseline,'linewidth',2,'LineStyle','-','Color',[0.25, 0.41, 0.88])
%plot(date,alt1,'linewidth',2,'LineStyle','-.','Color',[0.4, 0.69, 0.2])
%plot(date,alt2,'linewidth',1.5,'LineStyle','-.','Color',[0., 0.42, 0.24])
xlim([2018 2060])
line([2023 2023],[0 2000],'linewidth',1,'LineStyle',':','Color',[0, 0, 0])
%line([2031 2031],[0 2000],'linewidth',1,'LineStyle',':','Color',[0, 0, 0])
%legend('Cap with LRF at 2.2%','Cap reform (baseline)','Cap reform (with LRF at 2.2% beyond 2030)')
legend('Cap with LRF at 2.2% (pre-reform)','Cap with 2023 reform (baseline)')
set(gca,'fontname','times','FontSize',15)
grid on
ylabel('Mt CO2 e')
%title({'Cap scenarios','(in million permits)'},'Interpreter','Latex');


set(gcf,'PaperType','<custom>','PaperSize',[20.4 15.25],'PaperPosition',[0.1 0.02158 20.3046 15.2284]);
print -dpdf Figure_cap.pdf
