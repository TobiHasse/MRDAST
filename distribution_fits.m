function distribution_fits()
% Purpose:  This function will make the survivor function (CCDF) figures 8,
%           9 & 10 including the fits to the distribution.  It is reordered
%           and cleaned up code from the files curve_fits.m and
%           curv_fits_array.m both originally written in December 2018
%
%           Dependancies:  input files created by compute_dists.m
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     edited October 2021, written December 2018

% Load files and add CDF variables to an array of CCDF's aka survior
% functions or Complementary Cumulative Distribution Functions
load('TRH 205k upstream storage time cdf pdf.mat')
ccd(2,:)=1-CDFp;
ccd(3,:)=1-CDFv;
load('TRH 205k downstream storage time cdf pdf.mat')
ccd(4,:)=1-CDFp;
ccd(5,:)=1-CDFv;
% load this one last to use other variables in file for plots
load('TRH 205k ALL storage time cdf pdf.mat')
ccd(1,:)=1-CDFall;

fx=px(2:end);   % grab the x coordinates for plotting

% make figure with exponential and Pareto fit options
Ch2fig8_candidate_functions(fx, ccd(1,:), aveall, m70, mdall);

% obtain fit objects for all the survivor functions
[fits,~] = survivor_fit_eqn_8(fx,ccd,1);

% show the data and model fit all sediment
Ch2fig9_model(fits{1}, fx, ccd(1,:), aveall, m70, mdall);

% show the data and models for each facies and reach
Ch2fig10_all_four(fx,fits,ccd)
end

function Ch2fig10_all_four(fx,fits,ccd)
% This function plots the surviovor functions and fits for each reach and
% facies

% pull out variables for plotting
% second character o = fit;       d = data
% third  character p = point bar; v = vertically acreted floodplain
% fourth character u = upstream;  d = downstream
fopu = fits{2};
fdpu = ccd(2,:);
fovu = fits{3};
fdvu = ccd(3,:);
fopd = fits{4};
fdpd = ccd(4,:);
fovd = fits{5};
fdvd = ccd(5,:);

hall4=figure(96);

% plot and format first data set and fit
h1=plot(fovu,fx,fdvu);
% set(h(1),'marker','+','markeredgecolor','b','markersize',2)
set(h1(1),'marker','none','linestyle','-','color','b')
set(h1(2),'color','k','linestyle',':')

% plot and format additional sets
hold on
h2=plot(fopu,fx,fdpu);
% set(h(1),'marker','x','markeredgecolor','g','markersize',2)
set(h2(1),'marker','none','linestyle','-','color','g')
set(h2(2),'color','k','linestyle','-.')

h3=plot(fovd,fx,fdvd);
% set(h(1),'marker','s','markeredgecolor','y','markersize',2)
set(h3(1),'marker','none','linestyle','-','color','y')
set(h3(2),'color','k','linestyle','--')

h4=plot(fopd,fx,fdpd);
% set(h(1),'marker','o','markeredgecolor','r','markersize',2)
set(h4(1),'marker','none','linestyle','-','color','r')
set(h4(2),'color','k','linestyle','-')

hold off

% axis scale and labels
ylim([10^-5 10^0])
xlim([10 10^6])
set(gca,'yscale','log',...
    'xscale','log')
set(hall4,'color','w')
ylabel('Fraction older')
xlabel('Storage time (years)')
title(strcat('Curve fits for upstream and downstream storage time',...
    ' survivor functions') )
% axes equal

% reorder legend using vectors of legend handles and label strings
legend([h2(1),h2(2),h1(1),h1(2),h4(1),h4(2),h3(1),h3(2)],...
    ['Upstream point bar  ';'Fit                 ';...
    'Upstream flood      ';'Fit                 ';...
    'Downstream point bar';'Fit                 ';...
    'Downstream flood    ';'Fit                 '],'location','southwest')

% Set position and size of figure
pp = [354 374 510 489];
set(hall4,'position',pp,'PaperPositionMode','auto')

drawnow
outfile=sprintf('Upstream and Downstream survivor function fits');
print(hall4,'-painters', '-dpng', '-r600',outfile)

end


function Ch2fig8_candidate_functions(fx, CCDFall, aveall, m70, mdall)
%% This function will take a sruvivor function and attempt to fit it with
% Pareto
% Truncated Pareto
% Exponential
% Exponential on first 900 years of data
% Truncated Pareto on data over 900 years old
% 900 years is determined implicitly by using the 30th element of an array
% with a 30 year time step

% Pareto
ft= fittype('(a/x)^c')
sp = [50 .5];
[fobjp,gofp,outp] = fit(fx',CCDFall',ft,'startpoint',sp)

% Truncated Pareto
ft = fittype('c*a^c*(x^-c-b^-c)./( (1-(a/b)^c)) ')
sp = [138 170000 .9 ];
[fobjpt,gofpt,outpt] = fit(fx',CCDFall',ft,'startpoint',sp)
% Truncated Pareto over 900 years
[fobjpt9,gofpt9,outpt9] = fit(fx(30:end)',CCDFall(30:end)',ft,...
    'startpoint',sp)

% Exponential
ft = fittype('a*exp(b*x)')
sp = [1 -.0001];
[fobje, gofe, oute ] = fit(fx',CCDFall',ft,'startpoint',sp)
% Exponential under 900 years
[fobjet, gofet, outet ] = fit(fx(1:30)',CCDFall(1:30)',ft,'startpoint',sp)

mksz=8;
hfuncs = figure(97);

% plot and format first data set and fit
h=plot(fobjp,fx,CCDFall);
set(h(2),'color','k','linestyle','-')
set(h(1),'markeredgecolor','k')
set(h(1),'linestyle','none','marker','o','color','g','markersize',...
    mksz/2,'markerfacecolor',[.6 .6 .6],'markeredgecolor','none')

set(gca,'yscale','log',...
    'xscale','log')
set(hfuncs,'color','w')

% plot additional sets
hold on
h=plot(fobjpt,fx,CCDFall);
set(h(1),'markeredgecolor','none')
set(h(2),'color','k','linestyle',':')
h=plot(fobje,fx,CCDFall);
set(h(1),'markeredgecolor','none')
set(h(2),'color','k','linestyle','-.')
h=plot(fobjet,fx,CCDFall);
set(h(1),'markeredgecolor','none')
set(h(2),'color','b','linestyle','-')
h=plot(fobjpt9,fx,CCDFall);
set(h(1),'markeredgecolor','none')
set(h(2),'color','r','linestyle','-')

hold off

% optional add median, 70th percentile and residence time (average) markers
add_markers(fx,CCDFall,mdall,m70,aveall,mksz)


xlim([10^1 10^6])
ylim([.0000001 10])
% axis normal

title('Comparison of functional fits to survivor function')
xlabel('Storage time (years)')
ylabel('Fraction older')

legend('Fraction Older','Pareto fit','','Truncated Pareto fit','',...
    'Exponential fit','','Exponential fit < 900','',...
    'Truncated Pareto fit >900','Median','70^t^h Percentile','Mean',...
    'location','southwest')

% Set position and size of figure
pp= [ 354 156 510 784] %figure position for tall svg output with axes equal
set(hfuncs,'position',pp,'PaperPositionMode','auto')

drawnow
outfile=sprintf('Candidate function fits to All stored sediment');
print(hfuncs,'-painters', '-dpng', '-r600',outfile)

end % function

function Ch2fig9_model(fobj, fx, CCDFall, aveall, m70, mdall)
% Show the best fit model (equation 8) for the all sediment survior
% function

mksz=8;
hmod = figure(98);
clf
h=plot(fobj,fx,CCDFall);     % for best fit model
set(h(2),'color','k','linestyle','-')
set(h(1),'markeredgecolor','k')
set(h(1),'linestyle','none','marker','o','color','g','markersize',mksz,...
    'markerfacecolor',[.6 .6 .6],'markeredgecolor','none')

% optional add median, 70th percentile and residence time (average) markers
add_markers(fx,CCDFall,mdall,m70,aveall,mksz)

set(gca,'yscale','log','xscale','log')
set(hmod,'color','w')

xlim([10^1 10^6])
ylim([.0000001 10])
% axis normal
title('Equation 8 model functional fit to survivor function')
xlabel('Storage time (years)')
ylabel('Fraction older')
legend('Fraction Older','Equation 8','Median','70^t^h percentile',...
    'Mean','location','southwest')

% set size and aspect ratio of figure for 5.8 inch wide print area
pp= [ 354 156 510 784] %figure position for tall svg output with axes equal
set(hmod,'position',pp,'PaperPositionMode','auto')

drawnow
outfile=sprintf('Equation 8 model fit to All stored sediment');
print(hmod,'-painters', '-dpng', '-r600',outfile)

end % function Ch2fig8

function add_markers(fx,CCDFall,mdall,m70,aveall,mksz)
% this function will add markers for the median avearge and 70th percentile
% this code would be more robust if it used a figure handle

%plot median, find index of x coordinate
[~, idxd]=find(min(abs(fx-mdall))==abs(fx-mdall));
%plot 70th percentile
[~, idx7]=find(min(abs(fx-m70))==abs(fx-m70));
%plot mean
[~, idxa]=find(min(abs(fx-aveall))==abs(fx-aveall));

% add them to the current figure
hold on
h=plot(fx(idxd),CCDFall(idxd),...
    fx(idx7),CCDFall(idx7),...
    fx(idxa),CCDFall(idxa));

hold off

% format markers
set(h(:),'linestyle','none','markersize',mksz*2)
set(h(1),'marker','d','color','k')
set(h(2),'marker','s','color','k')
set(h(3),'marker','o','color','k')
end

