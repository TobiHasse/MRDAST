function Ch2fig15_meander_belt(closesolid,mbsolid,mb_start_step,i )
% Purpose:  This function will generate figures similar to figure 15
%           showing the age distribution of the floodplain nearby to the
%           channel
% Author:   Tobias Hasse    tobiack@udel.edu
% Date:     2017, edited June 2021

% old input files from when this file was a script
% cd 'C:\Users\User\Documents\MATLAB\MeanderBelt'
% load(TRH 205ky meander belt mb solidity and maps.mat)
% load('TRH 205ky meander belt mb solidity and maps July.mat')
% load(['TRH 205ky meander belt mb solidity and maps ',...
%     'Jan 2020 Complete 2.mat'])
%    
dt_yrs = 30;
confidence_interval = 95;
keyboard
htan = figure(22);
clf   % mb_start_step should replace references to 4000
%     csa=closesolid(4000:end,:); % tangent meander belt (complete data)
    csa=closesolid(mb_start_step:i,:); % tangent meander belt→partial data
    [ppx, ppy, px, m] = sub_plotSolid(csa, dt_yrs, confidence_interval); 
        patch(ppx,ppy,[.7 .7 .7],'edgecolor','none')

    hold on
    plot(px, 1-m,'color','k','linewidth',2)
    hold off

    set(gcf,'color','w')
    % tangent meander belt tick marks
    set(gca,'xtick',[0,1000,1500,5000,10000],... 
            'ytick',[0 .05 .1 1],...
            'xgrid','on',...
            'ygrid','on')
    set(gcf,'position',[   100   100   853   358 ])
    legend('95% Confidence interval','Median','delete')
    ylabel('Fraction of nearby floodplain older','fontsize',10)
    xlabel('Age of floodplain near channel','fontsize',10)
    title(['Solidity statistic for nearby floodplain computed ',...
        'with closesolid'],'fontsize',10)
drawnow
outfile=sprintf('Tangent meander belt age distribution. Time %d.png',i);
print(htan,'-painters', '-dpng', '-r600',outfile)

hcamp = figure(21);
clf
% csa=mbsolid(4000:end,:);  % Camporeale 2005 meanderbelt figure (complete)
csa=mbsolid(mb_start_step:i,:);  %Camporeale 2005 meanderbelt fig (partial)
    [ppx, ppy, px, m] = sub_plotSolid(csa, dt_yrs, confidence_interval); 
        patch(ppx,ppy,[.7 .7 .7],'edgecolor','none')

    hold on
    plot(px, 1-m,'color','k','linewidth',2)
    hold off

    set(gcf,'color','w')
    set(gca,'xtick',[0,1000,2000,5000,10000],...
            'ytick',[0 .05 .1 .25 1],...
            'xgrid','on',...
            'ygrid','on')
    set(gcf,'position',[   100   600   853   358 ])
    legend('95% Confidence interval','Median','delete')
    ylabel('Fraction of nearby floodplain older','fontsize',10)
    xlabel('Age of floodplain near channel','fontsize',10)
    title(['Solidity statistic for floodplain width of 3.4 lambda ',...
        'after Camporeale 2005'],'fontsize',10)
drawnow
outfile=sprintf('Camporeale meander belt age distribution. Time %d.png',i);
print(hcamp,'-painters', '-dpng', '-r600',outfile)

%%
% outf = sprintf('persistance of proximal old material');
% fig2svg(sprintf('%s tangent meander belt Jan.svg',outf)
% fig2svg(sprintf('%s 3.4 lambca meander belt Jan.svg',outf)
% fig2svg(sprintf('%s both meander belts Jan.svg',outf)
%%

% figure(11)
% imagesc(bool_channel(:,345:2240)+closeMB(:,345:2240))
% axis equal
% title(['closing the channel centerlin image with a 200 ',...
%     structuring element'])
% print('-painters', '-dpng', '-r600','closed channel 200.png')
end % function

