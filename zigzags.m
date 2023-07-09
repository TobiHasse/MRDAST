function zigzags()
% Purpose:  Show swath profiles in stacked plots similar to wheeler
%           diagrams to show channel position in the valley through time
%           This code is excerpted from xsteps.m
%           Create figure 14 from Chapter 2
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     Spring 2016, edited October 2021

% private functions dependencies
%           cmapsPUB.m
%           animation.m
%           chan_occ.m
%           slices.m
%           params_storage.mat
%           Chan_occ files from the complete model _display_ run


% make pill bug movie: get array of channel_frequency at each time interval
working_dir = pwd;              % get the current directory
disp('please select the folder with channel occupation slices')
% EDIT chan_freq TO RETURN ALL chan_ %{i}
[val_w_pxl , ~            ,chan_freq] = chan_occ('spam.gif',1);
% [val_w_pxl              ,chan_freq] = chan_occ('spam.gif',1);

% return to the working directory
% cd 'C:\Users\thasse\Documents\MATLAB\test' %
cd( working_dir )
%% selected swaths, normalized
%  yl = ylim
% stts = 575:500:2200;

%% Load previously used or recommended slices
% This code is useful for deciding the change from upstream to downstream
% run this code to get the swaths from the analyzed slices.


load params_storage.mat         % contains bounds for reaches, etc
% could run this code to get x_start_rx and x_end_rx
% Hasse_211* files available at: https://doi.org/10.5281/ZENODO.5651841.
% load('Hasse_211ka_30yr_A3_Cfo24_2Eo_offset.mat')
% load('Hasse_211ka_30yr_A3_Cfo24_2Eo.mat')
% [ x_start_rx, x_start_buffered, x_end_rx, end_sim_rx ] = ...
%                             starts_and_ends( riv, riv2, pix_per_chan, B);

% show the slices
show_slice_figures = true;
array_size=450*10^6;            % dissertation setting 
[x_starts, x_ends, ymn, ymx] = slices( array_size, end_sim, ...
    x_start, x_end, show_slice_figures);

% retrun to the working directory
% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
cd( working_dir )

% DEBUG save OR reload variables to save time when rerunning code
% save('chan_freq and slices.mat','chan_freq','x_starts','x_ends','-v7.3')
% load('chan_freq and slices.mat')

%% these are the swaths for the slices in the dissertation.  This was used
% to decide to use slices
% 1-16 in the upstream reach and
% 17 - 36 in the downstream reach
swath_w = x_ends-x_starts;      % vector of swath widths
yl2=nan;                        % y limit for probability density fig
for i = 1:4:length(x_starts)    % put 4 swaths on each figure
    end_swath = min(i+3, length( x_starts ) );
    % get the mean probability densiy within the swath
    swath_count = get_swath_counts( chan_freq, ...
        x_starts( i : end_swath ), swath_w( i : end_swath ) );
    % make a set of swath profiles
    yl2 = zigzag_fig( swath_count, x_starts( i : end_swath ),...
        display_progress_years, yl2 );
end

%% constant width slices
swath_width = 100;              % scalar for constant width swath profiles

stts = x_start:swath_width:x_end;   % vector of starting x coordinates
yl2=nan;                        % y limit for probability density fig

for i = 1:4:length(stts)
    end_swath = min(i+3,length(stts));
    % get the mean probability densiy within the swath
    swath_count = get_swath_counts( chan_freq, stts( i : end_swath ),...
        swath_width );
    % make a set of swath profiles
    yl2 =  zigzag_fig( swath_count, stts( i : end_swath ),...
        display_progress_years, yl2 )
end

%% exemplary swaths used in Chapter 2 figure 14
stts = 900:1100:2200;
swath_count = get_swath_counts( chan_freq, stts, swath_width );
yl2 = zigzag_fig( swath_count, stts, display_progress_years, nan)

end % function zigzags

function [swath_count] = get_swath_counts(chan_freq, stts, swath_w)
% Purpose   This function computes the average amount of time pixels within
%           the swath were occupied during each time interval.  The units
%           are the same as the units in chan_freq wich should be in the
%           percentage of each slice interval

if isequal(length(swath_w),1)   % make swath_w a vector
    swath_w = repmat(swath_w,1,numel(stts));
end

% starting with the last time interval
for j = numel(chan_freq):-1:1
    if numel(chan_freq{j})>0
        for i = length(stts):-1:1
            % compile the average probability density (%) in the swath
            swath_count(:,i,j) = ...
                mean( chan_freq{j}(:,stts(i):stts(i) + swath_w(i) - 1),2 );
        end
    end
end

end % function swath_counts


function yl2 = zigzag_fig( swath_count, stts, dt, yl2 )
% Purpose   This function will make swath profile figures and save them to
%           file in the current folder
% Inputs    swath_counts    the average % of time pixels were occupied by
%                           the channel
%           stts            starting down valley coordinate for each swath
%           dt              years represented by each chan_freq slcie
%           yl2             y axis upper limit for probability density fig

% make a boolean shwing if any pixel in each swath was occupied
valley_widths = double( flipud( permute( swath_count, [3 1 2] ) > 0 ));
lnth = size( valley_widths, 1 );    % # of slices

% assign ages to the values in the boolean
for i = 1:lnth
    valley_widths(i,:,:) = valley_widths(i,:,:) * dt * (lnth+1- i) ;
end

% cmapsTRA      % original call for colorbar which is slightly nicer than 
                % the color bar in cmapsPUB
% cmap_grayer   % original colormap from cmapsTRA

cm = cmapsPUB('age',0);     % get colormap data
cm_step = size( cm, 1 ) / size( valley_widths, 1 );
cidx = floor( 2 : cm_step : size( cm, 1 ) ); % indices for short cmap
cidx( end-6 : end ) = cidx( end-6 : end ) + [6 6 8 10 10 5 1]; % adjust


hf = figure(55);
clf

% set figure height to fit number of subplots
fig_h = 2.6 * size( swath_count, 2 );   
hf.Units = 'inches';
hf.Position = [1 1 5.8 fig_h];

lnth = size( swath_count , 2 );     % # of subplot rows
for j = 1:lnth                      % go through the swaths
    hR = subplot( lnth, 2, 2*j );   % subplot on the right
        % quick plot for testing 
        % ha = imagesc(flipud(permute( swath_count(:,i,:),[3 1 2])>0));
    
        im_this = valley_widths(:,:,j);
        % find region not yet explored
        outside = flipud( cumsum( flipud( im_this ) ) );
        % find region not re-explored at end of simulation
        abandoned = cumsum( im_this );
        % make regions numeric so you can in color
        im_this( ~abandoned ) = -dt;    
        im_this( ~outside   ) = -2 * dt;
        ha = imagesc(im_this);

        % burgundy blue white then rest of colors
        % colormap( hr, [ .5,0,0; 0,0,1; 1,1,1; cm( cidx, : ) ] ) 
        % gray lt gray white then rest of colors
        colormap( hR, [ .5,.5,.5; .6,.6,.6; 1,1,1; cm( cidx, : ) ] ) 
    
    hL = subplot( lnth, 2, 2*j-1 );  % subplot on the left
        % normalizer= sum(sum(permute(swath_count(:,i,:),[1 3 2]),2))
        % normalizer= sum(sum(permute(swath_count(:,j,:),[1 3 2]),2))
        % normalizer = 6845 % when chan_freq is number of steps occupied
        % when swath_count is average % of swath
        normalizer = size( swath_count, 3 ); 

        % coordinate for bottom edge of first swath
        x_bot =zeros( size( swath_count, 1 ), 1 ); 
        hp = area( [ x_bot, permute( swath_count(:,j,:), [1 3 2])]./...
            normalizer, 'edgecolor', 'none', 'linewidth', 0.3 );
%         colormap(hL, [1,1,1; cm(cidx,:)]) % white end % R2015a
        colororder(hL, [1,1,1; cm(cidx,:)]) % white end % R2021a
        
        try hL.YLim(2) = yl2; end   % set YLim unless yl2 is nan
    
    if isequal(j,1)             % for the top row of subplots only
        yl=get(hL,'ylim');
        title(hL,'Fraction of channel visits','fontsize',10)
        title(hR,'Valley locations visited during 5 ky','fontsize',10)
    end
    if isequal (j, lnth)        % for the bottom row of subplots only
        xlabel(hL,'Cross valley distance (km)')
        xlabel(hR,'Cross valley distance (km)')
    else                        % for other subplots
        set(hL,'xticklabel',[])
        set(hR,'xticklabel',[])
    end
    
%% legend entries are HARD CODED for dissertation settings    
    set(hR,'ytick',[1 18 41],...
        'yticklabel',[205 120 5],...
        'xtick',[500 1000 1500],...
        'xticklabel',[7 14 21],...
        'linewidth', 2, 'fontsize', 8 )
    set(hL,'ylim',[ yl(1) yl(2) ],...
        'xtick',[500 1000 1500],...
        'xticklabel',[7 14 21],...
        'linewidth', 2, 'fontsize', 8 )
    hL.TickLength  =  hL.TickLength*2;
    hR.TickLength = hR.TickLength*2;
    ylabel(hL,'Probability Density (% Sim Time)')
    ylabel(hR,'Simulation time (ky)')

    % make subplots larger, and reduce space in margins
    p = get(hL, 'pos');
    p = p + [-.06 -.02 .07 .055];
    set(hL,'pos',p);
    p = get(hR, 'pos');
    p = p + [-.005 -.02 .07 .055];
    set(hR,'pos',p);
    
    % label the subplot with the starting down valley pixel coordinate
    text(0.050 * hL.XLim( 2 ), 0.95 * hL.YLim( 2 ),...
        sprintf('Swath %d', stts( j ) ) )
    
    % set PaperPosition so that *.png file preserves figure size
    set( hf, 'color', 'w', 'PaperPosition', hf.Position );
    drawnow
    yl2 = hL.YLim( 2 ); % return yl2 for next set of swaths
end % for j
%%
outfile = sprintf('6845 normalized zig zags %s.png',num2str( stts ))
print('-painters', '-dpng', '-r600', outfile)
%print('-painters','-dpng','-r600','6845 normilized zig zags 900 2000.png')
end % function zigzag_fig
%% Inspect chan_freq

% for i = 1:numel(chan_freq)
%     figure(5)
%         imagesc(chan_freq{i})
%         colorbar
%         waitforbuttonpress
% end
