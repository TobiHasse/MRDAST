function ticker_plot(ticker_times,ticker_variables)
% Purpose:  This function will create a figure showing how long it is
%           taking for each subsection of the storage time model.  The
%           figure is in 4 panels and includes the cumulative times in the
%           upper left and the sparcity of the stratigraphic arrays in the
%           lower left.  The sparcity is the fraction of non zero data
%           points in the array.  There is significant sparcity and
%           changing the algorithm could be beneficial.
%           In the upper right is the additional time required at each step
%           to compute the storage time distribution.  Early in the
%           simulation it is faster to index into part of the arrays, later
%           it is faster to compute on the entire array.  If the algorithm
%           switch is in the right space the blue dots should increase
%           linearly and then level off.
%           In the lower right is for time to delete the eroded 
%           chronostratigraphy, same as upper right. The red dots are the 
%           change in (change in) computational time.
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     edited June 2021, written 2015

% ticker check function name
% fprintf('Welcome to: "%s" \n',...
%         fullfile(fileparts(mfilename('fullpath')),mfilename))

% optional read in from file, can be helpful to check while simulation is
% running on a different machine
% cd(uigetdir)
% load('tickers.mat')
% ticker_variables
% ticker_times(end,:)

nd = find(ticker_times(:,1),1,'last'); % find end of data

% plot the change in time increment for various timers
tick_diff = diff( [ 0,0,0; ticker_times( :, [ 3,6,4 ] ) ] );
tick_2_diff = diff( [ 0,0,0; tick_diff ] ); 

% remove large numbers if program is restarted
first_idx = find( tick_diff > 0 , 1 ,'first' ) ;  

y2 = max( tick_diff( first_idx+1 : end,:) ) / .8; % scale the plots nicely
y2(y2==0)=.5; % if y2 is not greater than y-axis min ylim() will error

figure(6)
subplot(2,2,1) % cumulative times
    set(gcf,'color','w')
    plot(ticker_times(1:nd,8),sum(ticker_times(1:nd,1:7),2),...
        ticker_times(1:nd,8),ticker_times(1:nd,3),...
        ticker_times(1:nd,8),ticker_times(1:nd,4),...
        ticker_times(1:nd,8),ticker_times(1:nd,1),...
        ticker_times(1:nd,8),ticker_times(1:nd,2),...
        ticker_times(1:nd,8),ticker_times(1:nd,5),...
        ticker_times(1:nd,8),ticker_times(1:nd,6),...
        ticker_times(1:nd,8),ticker_times(1:nd,7))
    legend('total time','permute','clear','island','rep','deposit',...
        'reset','clean','location','northwest')
    xlabel('model completion (%)')
    ylabel('run time (s)')
    title('Computation time for storage time analysis')

subplot(2,2,2) % storage time calculations
    set(gcf,'color','w')
    plot( ticker_times(1:nd,8),tick_diff(1:nd,1),...
        ticker_times(1:nd,8),tick_2_diff(1:nd,1),...
        'marker','.','linestyle','none')
    legend('computation time incriment',...
        'change in computation time increment','location','northwest')
    xlabel('model completion (%)')
    ylabel('additional time (s)')
    title(sprintf(['1^s^t & 2^n^d derivatives of \n',...
        'storage time computation time']))
    yl = ylim;
    ylim( [ 0 , y2(1) ] )

subplot(2,2,3) % sparcity of chronostratigraphy matricies
    set(gcf,'color','w')
    plot( ticker_times(1:nd,8),ticker_times(1:nd,12),...
        ticker_times(1:nd,8),ticker_times(1:nd,13),...
        'marker','.','linestyle','none')
    legend('Overbank sparsity','Point Bar sparsity','location','northwest')
    xlabel('model completion (%)')
    ylabel('fraction of nonzero elements')
    title('Sparsity of ')
    set(gca,'yscale','log')

subplot(2,2,4) % clear eroded sediment
      set(gcf,'color','w')
        plot( ticker_times(1:nd,8),tick_diff(1:nd,3),...
            ticker_times(1:nd,8),tick_2_diff(1:nd,3),...
            'marker','.','linestyle','none')
    legend('computation time incriment',...
        'change in computation time increment','location','northwest')
    xlabel('model completion (%)')
    ylabel('additional time (s)')
    title(sprintf(['1^s^t & 2^n^d derivatives of \n'...
        'clear eroded sediment computation time']))
    yl = ylim;
    ylim( [ 0 , y2(3) ] )
    drawnow
%%
end

%%
% for i = 1:20
%     sp16(i)=prctile(ticker_times(i:20:end,12),16);
%     sp50(i)=prctile(ticker_times(i:20:end,12),50);
%     sp84(i)=prctile(ticker_times(i:20:end,12),84);
% end
% plot(1:20,[sp84;sp50;sp16])
% legend('84','50','16','location','northwest')
