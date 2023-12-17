function Ch2fig13_pdf()
% Purpose:  This function will make the probability density function figure
%           similar to figure 13 from Tobias Hasse's dissertation
%
%           Dependancies:  input files created by compute_dists.m
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     edited October 2021, written January 2019

% Load files and add PDF variables to an array
load('TRH 205k ALL storage time cdf pdf.mat')
pdf(1,:)=PDFall;
load('TRH 205k upstream storage time cdf pdf.mat')
pdf(2,:)=PDFp;
pdf(3,:)=PDFv;
load('TRH 205k downstream storage time cdf pdf.mat')
pdf(4,:)=PDFp;
pdf(5,:)=PDFv;

fx=px(2:end); % grab the x coordinates for plotting

%% Make a scatter plot
hpdf = figure(103);
% set(h(1),'marker','none','linestyle','-','color','k')
% set(h(2),'marker','none','linestyle','-','color','b')
% set(h(3),'marker','none','linestyle','-','color','g')
% set(h(4),'marker','none','linestyle','-','color','y')
% set(h(5),'marker','none','linestyle','-','color','r')

h1=scatter(fx,pdf(5,:),15);
    set(h1,'marker','o','MarkerFaceColor',[255 156 0]./255,...
        'MarkerEdgeColor','none')
%     h.MarkerEdgeAlpha = .2;
%     alpha(.2)
hold on
h2=scatter(fx,pdf(4,:),15);
    set(h2,'marker','o','MarkerFaceColor','r','MarkerEdgeColor','none')
%     h.MarkerEdgeAlpha = .2;
%     alpha(.2)
h3=scatter(fx,pdf(3,:),15);
    set(h3,'marker','o','MarkerEdgeColor','b')   
%     h.MarkerEdgeAlpha = .2;
%     alpha(.2)
h4=scatter(fx,pdf(2,:),15);
    set(h4,'marker','o','MarkerEdgeColor','g')
%     h.MarkerEdgeAlpha = .2;
%     alpha(.2)
h5=scatter(fx,pdf(1,:),15);
    set(h5,'marker','o','MarkerFaceColor','k','markerEdgeColor','none')
%     h.MarkerFaceAlpha = .02;
hold off
%     alpha(0) % makes points invisible

% junk = 61
box on
set(gca,'yscale','log','xscale','log')
ylim([ 10^-6 10^-1])

% xlim([ 20 2000])
% ylim([ 10^-3 10^-1])
% legend entries in default order
% legend('Downstream flood','Downstream point bar','Upstream flood',...
%     'Upstream point bar','All sediment','location','southwest')
% reordered legend entries
legend([h5, h4, h3, h2, h1],...
    ['All sediment        ';...
    'Upstream point bar  ';...
    'Upstream flood      ';...
    'Downstream point bar';...
    'Downstream flood    '],'location','northeast')
ylabel('Fraction / 30 yr')
xlabel('Storage time (years)')
set(hpdf,'color','w','position',[40 40 500 500],'PaperPositionMode','auto')
set(gca,'position',[.1 .1 .8 .8])
% get(gca,'position')

drawnow
    outfile=sprintf('Probability density functions');
    print(hpdf,'-painters', '-dpng', '-r600',outfile)

end
