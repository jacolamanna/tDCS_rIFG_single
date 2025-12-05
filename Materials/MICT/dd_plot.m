function [k_low, k_hi, rr_low, rr_hi] = dd_plot(dell, imm_low, imm_hi, fign, low_err, hi_err)

fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[-Inf],...
               'Upper',[Inf],...
               'StartPoint',[0.5]);
ft_hi = fittype('10000 / (1+(k*x))','options',fo);
ft_low = fittype('500 / (1+(k*x))','options',fo);

% k [1/months] = k / 30 [1/days]
[c_low,g1] = fit(dell,imm_low,ft_low);
rr_low = g1.adjrsquare;
[c_hi,g2] = fit(dell,imm_hi,ft_hi);
rr_hi = g2.adjrsquare;
k_low = c_low.k;
k_hi = c_hi.k;

tx_low = ['k_{500€} = ',num2str(c_low.k,3),' [1/months] | R^{2}_{adj} = ',num2str(g1.adjrsquare,3),''];
tx_hi = ['k_{10k€} = ',num2str(c_hi.k,3),' [1/months] | R^{2}_{adj} = ',num2str(g2.adjrsquare,3),''];

dell2 = [dell(1):0.001:dell(end)];

% Fonts eccetera 
set(0, 'DefaultFigureVisible', 'on')
set(0, 'DefaultAxesFontName', 'Helvetica','DefaultAxesFontSize',14);
set(0,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
set(0,'DefaultTextFontSize',14);
set(0,'DefaultTextFontName','Helvetica');
set(0,'DefaultAxesLineWidth',1);

c1 = 'k';
c2 = [0 114 189]/255+0.1;
c3 = [59 114 46]/255+0.1;
c4 = [237 206 0]/255;

f = figure('visible','off');
if isempty(low_err)
    plot(dell,imm_low,'sb','MarkerSize',10);
else
    errorbar(dell,imm_low,low_err,'sb','MarkerSize',10,'linewidth',1);
end
hold on, plot(dell2, c_low(dell2),'-b','linewidth',1)

title(tx_low);
xlabel('Reward delay [months]');
ylabel('Subjective value [€]');

box off
set(gca,'TickDir','out');

saveas(f,['',fign,'_LOW.fig']);

f = figure('visible','off');
if isempty(hi_err)
    plot(dell,imm_hi,'sr','MarkerSize',10);
else
    errorbar(dell,imm_hi,hi_err,'sr','MarkerSize',10,'linewidth',1);
end
hold on, plot(dell2, c_hi(dell2),'-r','linewidth',1)

title(tx_hi);
xlabel('Reward delay [months]');
ylabel('Subjective value [€]');

box off
set(gca,'TickDir','out');

saveas(f,['',fign,'_HI.fig']);

salvaFigura(['',fign,'_LOW'],5,4);
salvaFigura(['',fign,'_HI'],5,4);
