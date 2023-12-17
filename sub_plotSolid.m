function [ppx ppy px m] = sub_plotSolid(csa,dt_yrs,cInterval)
px=[1:size(csa,2)]*dt_yrs; 
% figure(1) % plot all the data
%     plot(1-csa','marker','.','linestyle','none')
m=median(csa);
% mn = mean(csa);
% Y = prctile(X,p)
p2_5 = prctile(csa,50-cInterval/2);         % Dec 2023 revert to this
p97_5 = prctile(csa,50+cInterval/2);
% p2_5 = prctile_TRH(csa',50-cInterval/2)'; % Dec 2023 cancel this
% p97_5 = prctile_TRH(csa',50+cInterval/2)';
%
% figure(1) % plot the mean and median lines onto all the data
%     hold on
%     plot(1-m,'linewidth',5)
%     plot(1-mn,'linewidth',3,'color','r')
% 
%     hold off
%
% figure(2) % plot the confidence bounds
%     plot(px,1-m)
%     hold on
%     plot(px,1-p2_5,'color','r')
%     plot(px,1-p97_5,'color','g')
%     % plot(px,p75-p25,'color','k')
%     hold off
% 
%     ylabel('Fraction of nearby floodplain older than age')
%     xlabel('Age of floodplain near channel')
%     set(gcf,'color','w')
%     title('Solidity statistic for nearby floodplain computed with closesolid')

%
% figure(4)

    % h=area(px,[1-p95',1-p5'-(1-p95')],'linestyle','none');
    % set(h(1),'FaceColor','w')
    % set(h(2),'FaceColor',[.5 .5 .5])
    ppx=[px,fliplr(px)];
    ppy=[1-p2_5,fliplr(1-p97_5)];

end