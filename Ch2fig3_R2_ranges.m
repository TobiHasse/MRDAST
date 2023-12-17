function Ch2fig3_R2_ranges(RS)
% Purpose:  This function will display the the changes in R^2 as the scale 
%           break transition parameter xo (called T in the dissertation) 
%           changes
%           Chapter 2 figure 3 from Tobias Hasse's dissertation
%           It is excerpted code from the files curv_fits_array.m 
%           originally written in December 2018.
%           
%           Dependancies:   input files created by distr_r_squared_ranges.m
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     edited October 2021, written December 2018
rs=RS{1};
hrng = figure(71);
h=plot(rs(6,:),rs(1:5,:)');
ylim([.99 1])
xlim([.999*min(rs(6,:)) 1.001*max(rs(6,:))])
hl=legend('All eroding sediment','Upstream Point Bar',...
    'Upstream Flood','Downstream Point Bar','Downstream Flood');
xlabel(strcat('Transition from exponential to Pareto distribution',...
    ' (Storage time in years)'))
ylabel('R^2 of best fit model given T') % T == xo
set(hrng,'color','w')
set(hl,'position',[.3 .35 .336 .233]) % move legend near center of fig
set(h(1),'color',[.6 .6 .6])
set(h(2),'color','k','linestyle',':')
set(h(3),'color','k','linestyle','-.')
set(h(4),'color','k','linestyle','--')
set(h(5),'color','k','linestyle','-')

set(hrng,'PaperPositionMode','auto')

drawnow
outfile=sprintf('R squared variation and range with transition parameter');
print(hrng,'-painters', '-dpng', '-r600',outfile)

end % Ch2fig3_R2_ranges
