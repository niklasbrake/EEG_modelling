lambda_E = @(p) pdf('logn',p,log(0.5),1);
lambda_I = @(p) pdf('logn',p,log(2.5),1);
tauE = @(p) pdf('unif',p,1,3.5);
tauI = @(p) pdf('unif',p,5,20);
erev = @(p) pdf('unif',p,-75,-45);
gleak = @(p) pdf('unif',log10(p),-5,log10(0.005));


t1 = 10.^linspace(-3,2,1e3);
lex = lambda_E(t1); lex = lex/max(lex);
lix = lambda_I(t1); lix = lix/max(lix);

t2 = linspace(0,25,1e3);
tex = tauE(t2); tex = tex/max(tex)*0.6;
tix = tauI(t2); tix = tix/max(tix)*0.6;

t3 = linspace(-75,-40,1e3);
elx = erev(t3); elx = elx/max(elx)*0.6;

t4 = 10.^linspace(-5,-2,1e3);
glx = gleak(t4); glx = glx/max(glx)*0.6;


figureNB(3.3,2.8);
axes('Position',[0.25,0.8475,0.63,0.11]);
    fill([t1,fliplr(t1)],[lex,0*t1]+1.3,'k');
    hold on;
    fill([t1,fliplr(t1)],[lix,0*t1],'k');
    set(get(gca,'YAxis'),'Visible','off')
    set(gca,'xscale','log');
    xlim([1e-2,1e2])
    xticks([0.1,1,10]);
    xticklabels([0.1,1,10]);
    xax = get(gca,'xaxis');
    set(gca,'XMinorTick','on')
    xax.MinorTickValues = [0.1,10];    
    ylim([0,2.3]);
    gcaformat;
    set(gca,'FontSize',6)
    text(16e-4,2.25,'\lambda_E','FontSize',7);
    text(16e-4,0.25,'\lambda_I','FontSize',7);
    txt = text(10,-1.6,'   Hz','FontSize',6);
axes('Position',[0.25,0.585,0.63,0.11]);
    fill([t2,fliplr(t2)],[tex,0*t2]+1.3,'k');
    hold on;
    fill([t2,fliplr(t2)],[tix,0*t2],'k');
    set(get(gca,'YAxis'),'Visible','off')
    xlim([0,25])
    ylim([0,2.3]);
    xax = get(gca,'xaxis');
    set(gca,'XMinorTick','on')
    xax.MinorTickValues = [5:10:25];    
    gcaformat;
    set(gca,'FontSize',6)
    text(-5,2.25,'\tau_E','FontSize',7);
    text(-5,0.25,'\tau_I','FontSize',7);
    txt = text(20,-1.6,'   ms','FontSize',6);
axes('Position',[0.25,0.373,0.63,0.06]);
    fill([t3,fliplr(t3)],[glx,0*t3],'k');
    set(get(gca,'YAxis'),'Visible','off')
    xlim([-80,-40])
    xticks([-70,-50]);
    xax = get(gca,'xaxis');
    set(gca,'XMinorTick','on')
    xax.MinorTickValues = [-80,-60,-40];    
    gcaformat;
    set(gca,'FontSize',6)
    ylim([0,1]);
    text(-80,0.5,'E_L  ','FontSize',7,'HorizontalAlignment','right');
    txt = text(-49,-1.25,'   mV','FontSize',6);
axes('Position',[0.25,0.16,0.63,0.06]);
    fill([t4,fliplr(t4)]*1e3,[glx,0*t4],'k');
    set(get(gca,'YAxis'),'Visible','off')
    set(gca,'xscale','log');
    xax = get(gca,'xaxis');
    xax.Limits = [0.001,10];
    set(gca,'XMinorTick','on')
    xax.MinorTickValues = [0.001,0.01,0.1,1,10]      
    xticks([0.01,1])
    xticklabels([0.01,1])
    gcaformat;
    set(gca,'FontSize',6)
    ylim([0,1]);
    text(1e-3,0.5,'g_L  ','FontSize',7,'HorizontalAlignment','right');
    txt = text(1,-1.125,'  mS/cm^2','FontSize',6);