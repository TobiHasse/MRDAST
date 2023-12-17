function Ch2fig11_fig12_ranges(ranges5, coef)
% Purpose:  This function will display the the central estimate and lower
%           and upper confidence bounds of fit parameters from both matlab
%           fit and the R^2 based method
%           It is excerpted code from the file curv_fits_array.m 
%           originally written in December 2018.
%           
%           Dependancies:   input files created by distr_r_squared_ranges.m
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     edited October 2021, written December 2018
hrng4 = figure(101);
set(hrng4,'color','w')
order=[1 2 4 3 5];              % vector to reorder variables
for j = 1:5
    i = order(j);               % work in my order
    try; subplot(2,2,j); end    % make a subplot 1-4
    if i == 5;                  % make separate fig
        hrngT = figure(102);
        set(hrngT,'color','w'); 
    end    

    % plot the matlab fit ranges
    h=plot(2:2:10,[coef.intervalLo(:,i),coef.vals(:,i),...
        coef.intervalHi(:,i)],...
        'marker','x','linestyle','none','markersize',5,'color','k');
    set(h(1:2:3),'marker','.')

    hold on
    % plot the R^2 ranges
    h=plot(1:2:9,ranges5(:,:,i),...
        'marker','o','linestyle','none','markersize',5,'color','k');
    set(h(1:2:3),'marker','.')
    % plot error bars for matlab fit
    for k = 1:5
        px = 2*k;
        plot([px,px],[coef.intervalLo(k,i),coef.intervalHi(k,i)],...
            'color','k');
    end
    % plot error bars for R^2 ranges
    for k = 1:5
        px = 2*k-1;
        plot([px,px,px],ranges5(k,:,i),'color','k');
    end
    hold off
    
    % axis tight
    if i == 1;    ylim([60000 250000]);     end
    if i == 2;    ylim([.5 1.1]);           end
    if i == 3;    ylim([.8 1.05]);          end
    if i == 4;    ylim([-.00179 -.00074]);  end
    if i == 5;
        ylim([0 2000]);
        legend('max','Matlab fit range','min','max','R^2 based range',...
            'location','northeast')
    end
    xlim([0 11])
    ylabel(sprintf(' Coefficient: %s ', coef.names{i} ))
    xlabel('All | up: bar, flood | down: bar, flood','fontsize',8)
    % waitforbuttonpress
    

end
    set(hrng4,'PaperPositionMode','auto')
    set(hrngT,'PaperPositionMode','auto')

    drawnow
    outfile=sprintf('R squared based parameter ranges');
    print(hrng4,'-painters', '-dpng', '-r600',outfile)
    
    outfile=sprintf('R squared based transition parameter T ranges');
    print(hrngT,'-painters', '-dpng', '-r600',outfile)


end %Ch2fig11_fig12_ranges
