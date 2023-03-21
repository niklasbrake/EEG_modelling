function temp_plot(asynch_high_res,params1)


figureNB;
[full_model,synFun,apFun] = fittingmodel('avalanches');
% f_high_res = 0.1:0.1:1e3;
f_high_res = 0.5:0.5:1e3;
f_high_res = f_high_res(f_high_res<50);
plot(f_high_res,asynch_high_res,'k','LineWidth',1);
hold on;
plot(f_high_res,10.^apFun(f_high_res,params1),'color','g','LineWidth',1);
plot(f_high_res,10.^full_model(f_high_res,params1),'color','r','LineWidth',0.5);
plot(f_high_res,10.^synFun(f_high_res,params1),'color','b','LineWidth',1);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([0.1,100])
xticks(10.^[-1:2]);
xticklabels([0.1,1,10,100]);
xlabel('Frequency (Hz)');
ylabel('log powere');
yticks([]);