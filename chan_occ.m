function [ valley_w_pxl, occ_bool, chan_freq_out ] = ...
                                chan_occ( outfile,   slice_interval )
% Purpose:  Create an animated .gif showings how frequently the channel is
%           present in each pixel on the floodplain after reading in the 
%           data from file.  The user can choose how many files of data are
%           combined by adjusting the slice interval.  

% Author:   Tobias Hasse tobiack@udel.edu

% Date:     Spring 2016, edited June 2021

% INPUT     outfile          output filename for animation
%           slice_ineterval  how many intervals to combine in one movie frame

% OUTPUT    valley_w_pxl     valley widht in pixels
%           animation.gif    saved to current folder

% private functions dependencies
%           cmapsPUB.m
%           animation.m

% Parsing the file names:
% 2d .mat files must be in the current folder with a filename that begins 
% 'channel_occ' and includes 
% "**r*##.**" where 'r' is the first 'r' in the filename, and '##' 
% represents numerical digits (representing the model time step).
% There must be one character between the r and first 

fprintf('Welcome to: "%s" \n\n',...
        fullfile(fileparts(mfilename('fullpath')),mfilename))

% load some parameters for the simulation
load params_storage.mat 
load params_meander B

% cd (['C:\Users\User\Documents\MATLAB\MyLib\PUBlication\Ch2\other\',...
%     'display_model complete run'])
cd('C:\Users\thasse\Documents\MATLAB\test')
fprintf(['Please select the folder containing the *.mat files with the',...
    ' number of times \n the channel occupied each floodplain pixel',...
    '\nThese files should be called "channel_occ****.mat"\n\n'])
cd(uigetdir)

% Axes are labeled assuming that channel width is in meters and w_px wide 
w_px = pix_per_chan;            % width of channel in pixels                                   
d_px = 2 * B / 10;              % one tenth the channel width

cdata = cmapsPUB('topo',0);     % get colormap data 
cdata(end,:)=[1 1 1];           % WHITE

files = dir;
for i = numel(files):-1:1       % remove extra files starting at the end
    if numel( files( i ).name) > 20 
        if strcmp( 'channel_occ', files(i).name(1:11) ) %19
            %do nothing
        else
            files(i) = [];      % remove from list
        end
    else
        files(i) = [];          % remove from list
    end
  
end

total_slices = numel(files);    % of time steps in interval
% in case the user selects more slices than available:
figure_interval = min( slice_interval , total_slices );
% denom scales the pixel counts to be a fraction of the interval
denom = figure_interval * display_progress_years / time_step_years;

[ ~, ids ] = sort( [ files(:).datenum ] );  % sort the files
first_time = true;

for i = 1:numel( ids )
    load( files( ids( i ) ).name )          % load file
    print_name = strrep( files( ids( i ) ).name, '_' , ' ' );   % fix name
    try
        % try adding current file to recent files
        % result is percentage of interval pixels were occupied by channel
        chan_occ_accum  = chan_occ_accum + chan_occupation / denom * 100 ;
    catch
        % result is percentage of interval pixels were occupied by channel
        chan_occ_accum  = chan_occupation / denom * 100 ; 
        occ_bool        = false( size( chan_occ_accum ) ) ;
    end
    % for the first frame, set up the dummy figure for initialization
    if first_time  
        hf =    figure( 1 );
        hf.Visible =  'off'; % hide the dummy figure
            imagesc( chan_occupation/denom*100 )
            set(gcf,'color','w','position',[1  41   1280   907])
            colorbar
            axis equal
        cx = caxis * figure_interval;        % for consistent axis limits
        cx = round( cx( 2 ) / 4 );
        nominal_k_age_start = 0;            % starting age of interval
    end
    
    % create a boolean of the floodplain visited within the time interval
    occ_bool( chan_occ_accum > 0 ) = true;  
    % occ_bool( chan_occupation > 0 )= occ_bool( chan_occupation > 0 ) + 1;

    figure(2)
    clf
        imagesc(occ_bool)
        title(sprintf('max %f, i is %d',max(chan_occupation(:)),i))
        xlabel('Pixels down valley')
        ylabel('Pixels accross valley')
        colormap(parula)
        colorbar
        axis equal
        drawnow
        
    % at the end of the interval: make figure, reset and store variables 
    if ~mod(i,figure_interval) | isequal( i , numel(files) )
        chan_{i}=chan_occ_accum;            % store the data
        chan_occ_accum = [];                % reset the variable

        % parse the file name for the model step portion, scale by 
        % model step, round
        nominal_k_age = round( str2num( print_name( ...
            find ( print_name == 'r' , 1 , 'first' ) + 2 : ...
            find ( print_name == '.' , 1 , 'first' ) ) ) / ...
            1000 * time_step_years ) ;
        % hard coded method
        %         nominal_k_age = nominal_k_age + 5;
        
        hf = figure( 1 ); % make the pretty figure
        clf
        hf.Visible = 'on';
            imagesc( chan_{i} )
            set(gcf,'color','w','position',[1   41   1000    700])
            h=colorbar;
            colormap(flipud(cdata))
            caxis([-cx/length(cdata) cx])
%             ylabel(h,'Number of timesteps pixel occupied by channel',...
%                                                     'fontsize',14 )
            ylabel(h,['Percent of simulation interval pixel occupied ',...
                'by channel']                        ,'fontsize',14 )
            xl = xlim;
            yl = ylim;
            % make axis ticks in kilometers but scaled to channel width
            set(gca,'fontsize',14,...
                'xtick',     [0 : (w_px*100) : xl(2)                ],...
                'XTickLabel',[0 : d_px       : xl(2)*d_px/(w_px*100)],...
                'ytick',     [0 : (w_px*100) : yl(2)                ],...
                'yTickLabel',[0 : d_px       : yl(2)*d_px/(w_px*100)])
            title(sprintf('Channel occupation %d to %d k.y.',...
                nominal_k_age_start , nominal_k_age) ,'fontsize',14 )
            xlabel('Down valley distance (km)'       ,'fontsize',14 )
            ylabel('Cross valley distance (km)'      ,'fontsize',14 )
            axis equal
            drawnow

        nominal_k_age_start = nominal_k_age;    % reset for next image
        animation_TRH(outfile,hf,first_time,2)    % send data to animation
        first_time = false;
%         if isequal(i,6) % DEBUG quit function early
%             return
%         end
    end % if slice_interval        
        
end
    
valley_w_pxl = yl(2)-yl(1);
chan_freq_out = chan_;%{i};


end % function chan_occ
