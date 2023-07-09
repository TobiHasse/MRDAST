function display_model(riv, riv2, B)
% Purpose:  This function will take the meandering river simulation and
%           show some characteristics to the modeler prior to running a 
%           full storage time analysis.  While this function is slow, it 
%           can save the modeler consederable time by helping you consider 
%           how and where to compute the storage time simulations.
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     edited June 2021, written 2015
fprintf('Welcome to: "%s" \n',...
        fullfile(fileparts(mfilename('fullpath')),mfilename))
time_now = datestr(clock)

load params_storage.mat
load params_deposition.mat
end_sim = min(end_sim,numel(riv))
% display model progress and make channel occupation files increment
disp_progress = round( display_progress_years / time_step_years )  
% disp_progress = 10            % DEBUG option display progress frequently
display_stops = [ disp_progress : disp_progress : end_sim,...
                    end_sim ]; % vector including last step of model

%% find x y limits of simulation area
limits = domain_limits(riv,riv2,pix_per_chan,B);
% limits for slice to extract stratigraphic cross section
ymn=400; ymx=1000; x_start = 1000; x_end = 1100; % strat xsec
% ymn=300; ymx=800; x_start = 1000; x_end = 1100; % strat xsec
% keyboard
%% make some arrays for storing variables
fp_area        = zeros( end_sim, 7 );   
chan_area      = zeros( end_sim, 1 );   
mass           = zeros( end_sim, 1 );
mig_bins       = 0:50;
mig_distance   = zeros( size( mig_bins , 1 ), size( mig_bins , 2 ) - 1 );
    
algorithm_start = clock
%% Make boolean array of channel location
% Initialize the array (needed outside of loop one time)
% doIntermp_TRH edited MATLAB function to return desired values
[ x, y ]     = doInterpm_TRH( riv( 1 ).Xcl, riv( 1 ).Ycl,...
                    pixel_sz/interp_nds_pix, 'lin' ); 
% bool_channel is a raster the size of the model domain, channel = true
bool_channel = node2grid( x, y, pixel_sz, limits); 
% update x & y limits, new since dissertation
ymn     = max( ymn,     floor( size( bool_channel, 1 ) / 2 ) );
ymx     = min( ymx,            size( bool_channel, 1 ) );
x_start = max( x_start, floor( size( bool_channel, 6 ) / 2 ) );
x_end   = min( x_end,          size( bool_channel, 2 ) );
%% Make some 2D arrays for tracking model perameters
prev_channel = false( size( bool_channel ) ); % true if visited by channel
prev_chan_analized = prev_channel;  % true if visited after t = 4000

chan_occ_one = zeros( size( bool_channel ) ); % counts # steps channel is 
chan_occ_all = chan_occ_one;                  % present, whole simulation

birth_step     = chan_occ_one; % time step of most recent pt bar deposition
old_elev       = chan_occ_one; % floodplain topo, previous step
new_elev       = chan_occ_one; % floodplain topo, current step
accum          = chan_occ_one; % accumulation of new sediment this step
bool_pt_bar    = prev_channel; % boolean of current point bar location
fp_old_and_new = prev_channel; % boolean of flood plain not eroded or chan

bool_xsec      = prev_channel; % to display cross section location
bool_xsec(499:751,1061-5:1061+5) = true; % for cross section location
bool_xsec(1795:1805,200:557)     = true; % for scale bar 357 pix is 5 km
bool_xsec = bool_xsec(ones(size(prev_channel))); % shrink to fit if needed

%% Make some 3D arrays for stratigraphy and a boolean to 'erode' data from
% the stratigraphic stack
stratigraphy   = zeros( ymx-ymn+1 , x_end-x_start+1 , end_sim ); 
bool_erosion   = false( size( stratigraphy ) ); 

%% Distance transform of channel offset so 0 is the edge of channel pixels
old_fp_distance = double( bwdist( bool_channel ) - ...
                    floor( pix_per_chan / 2) );  
old_fp_distance( old_fp_distance < 0 ) = 0; % make all channel pixels = 0

% create an initially full floodplain 'terrace'
old_elev( old_fp_distance > 0 ) = 1; 

% the down valley range 345:2240 is used throughout this code, it is hard
% coded to replicate the range of the upstream and downstream reach in
% Tobias Hasse's dissertation.
mass(1) = sum( sum( old_elev(:,345:2240) ) ); % full pixels on floodplain
loop_start = datestr(clock)

% Loop through the simulation
for i = 2 :1: end_sim % USUALLY start at i = 2
    % channel node spacing might be larger than pixels, interpolate more
    [ x, y ]        = doInterpm_TRH( riv( i ).Xcl, riv( i ).Ycl,...
                        pixel_sz/interp_nds_pix, 'lin' );
    % take the x y points and create pixel map with the channel = ture
    bool_channel    = node2grid( x, y, pixel_sz, limits);
    % Distance from channel edge offset so 0 is the edge of channel pixels
    new_fp_distance = double( bwdist( bool_channel ) - ...
                        floor( pix_per_chan / 2) );  
    new_fp_distance( new_fp_distance < 0 ) = 0; % make all chan pixels = 0

    % riv2 holds the river planform half way between riv(i-1) and riv(i)
    [ x, y ]        = doInterpm_TRH( riv2( i-1 ).Xcl, riv2( i-1 ).Ycl, ...
                        pixel_sz/interp_nds_pix, 'lin' );   % interpolate
    bool_channel    = node2grid( x, y, pixel_sz, limits);   % make raster
    mid_fp_distance = double( bwdist( bool_channel ) - ...
                        floor( pix_per_chan / 2) );  
    mid_fp_distance( mid_fp_distance < 0 ) = 0;    

    % At cutoff, point bar complexes appear as 'islands' on the floodplain.
    % When the channel migrates further than one channel width other
    % artifacts can appear as 'islands'.  These artifacts appear as
    % 'crescent' shapes, sometimes attatched to a point bar complex like a
    % string on a 'balloon'
    fp_old_and_new = rm_islands9( new_fp_distance, mid_fp_distance, ...
                        old_fp_distance, px_size_thresh, pix_per_chan );
    % make a 3D array of recently eroded or channel locations
    bool_erosion   = repmat( ~fp_old_and_new(ymn:ymx,x_start:x_end), ...
                        [1,1,end_sim]);
    % Erode the stored stratigraphy, two methods optimize computation speed
    if i < .55 * end_sim % previously 2700, (3x faster than permute method)
        stratigraphy( bool_erosion(:,:,1:i) ) = 0; % delete the eroded data
    else  % slow and steady( faster for large arrays )
        stratigraphy( bool_erosion          ) = 0; % delete the eroded data
    end
    old_elev( ~fp_old_and_new ) = 0;               % delete the eroded data
    bool_pt_bar( ~ (fp_old_and_new | new_fp_distance == 0 ) ) = 1;       

    % calculate the thickness of newly accumulated sediment: pseudocode:
    % newly formed point bars +
    % remaining acommodation space with exponential decay (vertical)
    % pelagic + diffusive sedimentation * expon decay (horizontal distance)
    accum = bool_pt_bar * pt_bar_elev + ...     
        ( 1-old_elev ) ./ exp( old_elev ) .* ...
        (Nu + Mu * exp( -new_fp_distance / ( 2*B/pixel_sz * lambda ) ) );           
    % compute new elevation and change in elevation
    accum( new_fp_distance == 0 ) = 0;    % cut new channel b.c. exp(0) = 1
    new_elev = old_elev .*(1-bool_pt_bar) + accum;  % find new elevation
    new_elev( new_fp_distance == 0 ) = 0; % cut new channel b.c. exp(0) = 1
    
    birth_step( bool_pt_bar ) = i;     % time step of most recet deposition
    birth_step( new_fp_distance == 0 ) = 0; % blank out the current channel
    
    % record the accumulated 'sediment' in the stratigraphy
    stratigraphy( :,:,i ) = accum( ymn:ymx, x_start:x_end );

    %% record channel and active floodplain area and volume
    % count how many times the channel is present at each pixel
    chan_occ_one = chan_occ_one + ...
                    single( ( ~fp_old_and_new - ~new_fp_distance ) | ...
                    ~old_fp_distance );
    % the upstream reach ranges from 345:1324
    % the downstream reach ranges from 1325:2240
    chan_area(i) = sum( sum( new_fp_distance( : , 345 : 2240 ) == 0) );
    mass(i)      = sum( sum( old_elev       ( : , 345 : 2240 )     ) );
    hmig = histogram( old_fp_distance( bwperim( ~fp_old_and_new ) ),...
                                                            mig_bins );
    mig_distance = mig_distance + hmig.Values;
    
    prev_channel(       ~fp_old_and_new ) = true;   % exclude the terrace
    prev_chan_analized( ~fp_old_and_new ) = true;   % exclude the terrace
    if isequal( i,4000 )  % chuck the spinup area
        prev_chan_analized = false(size(prev_channel) );
        fp_area(1:i,5:7)= 0;
    end
    pc = sum(prev_channel);             % area previously visited
    fp_area(i,1) = sum(pc);             % whole area
    fp_area(i,2) = sum(pc( 345:2240));  % analyzed area
    fp_area(i,3) = sum(pc( 345:1324));  % upstream reach
    fp_area(i,4) = sum(pc(1325:2240));  % downstream reach

    pc = sum(prev_chan_analized);       % area visited since reset
    fp_area(i,5) = sum(pc( 345:2240));  % analyzed area since reset
    fp_area(i,6) = sum(pc( 345:1324));  % upstream analyzed area
    fp_area(i,7) = sum(pc(1325:2240));  % downstream analyzed area

% reset arrays for next iteration
    bool_pt_bar(:,:)= 0; 
    old_fp_distance = new_fp_distance;
    old_elev        = new_elev;

% plot progress incrementally or at the end, save output
    if any( display_stops == i ) 
        pct = i / display_stops(end)*100;
        frac_complete = pct / 100 
        disp(sprintf('Starting figures, \nCurrent time:  %s',...
            datestr(datenum(clock))))
        
        % compile recent channel occupations & channel occ from whole sim
        chan_occ_all = chan_occ_all + chan_occ_one;
%       figure showing number of times pixels are occupied by the channel
        hocc = figure(1);  
            set(hocc,'color','w','position',[0 800 560 800]);
            subplot(2,1,1)
                imagesc( chan_occ_one )
                title(sprintf(['Channel occupation most recent %d ',...
                    'time steps'] ,disp_progress))
                xlabel('Down valley distance (pixels)')
                ylabel('Cross valley distance (pixels)')
                cmapsPUB('topo',1);
                h=colorbar;
                cx = caxis;
                caxis([-cx(2)/length(colormap) cx(2)]);
                ylabel(h,'Number of time steps pixel occupied by channel')
                axis equal
                drawnow
            subplot(2,1,2)
                imagesc( chan_occ_all )
                xl = xlim; yl = ylim;
                set(gca,'xtick',[0:500:xl(2)],'ytick',[0:500:yl(2)],...
                    'xticklabel',[0:100:(xl(2)/pix_per_chan)],...
                    'yticklabel',[0:100:(yl(2)/pix_per_chan)])
                title(sprintf('Channel occupation after %d time steps' ,i))
                xlabel('Down valley distance (channel widths)')
                ylabel('Cross valley distance (channel widths)')
                cmapsPUB('topo',1);
                h=colorbar;
                cx = caxis;
                caxis([-cx(2)/length(colormap) cx(2)]);
                ylabel(h,'Number of time steps pixel occupied by channel')
                axis equal
                drawnow
%       figure showing the total area visited so far (various regions)
        harea = figure(2); 
            set(harea,'color','w','position',[0 300 560 420]);
            plot(fp_area)
            title(sprintf('Time step %d ' ,i))
            xlabel('Simulation steps')
            ylabel('Number of pixels previously occupied by channel')
            leg = {'all';'analyzed area';'upstream';'downstream';...
                    'analyzed time';'upstream';'downstream'};
            legend(leg, 'location','northwest');
%       the age (birth step) of most recent point bar deposition       
        hage = figure(3); 
        clf
            set(hage,'color','w','position',[500 460 560 520]);
            imagesc( birth_step ) 
            xl = xlim; yl = ylim;
            set(gca,'xtick',[0:500:xl(2)],'ytick',[0:500:yl(2)],...
                'xticklabel',[0:100:(xl(2)/pix_per_chan)],...
                'yticklabel',[0:100:(yl(2)/pix_per_chan)])
            title(sprintf('Time step at deposition (age), time step %d',i))
            xlabel('Down valley distance (channel widths)')
            ylabel('Cross valley distance (channel widths)')
            cmapsPUB('age',0);
            h=colorbar;
                set(h,'Location','southoutside')
                ylabel(h,sprintf(strcat('Time step at deposition',...
                    ' (* %d for years)'), dt_save_years*dt) )
            axis equal
            drawnow
%       topography of floodplain surface with scale bar a and cross section
        htopo = figure(4); 
        clf
            set(htopo,'color','w','position',[500 700 560 520]);
            % new_elev is topo, bool_xsec adds scale bar and cross sec loc
            imagesc( new_elev -new_elev.*bool_xsec + bool_xsec*.97 )
            xl = xlim; yl = ylim;
            set(gca,'xtick',[0:500:xl(2)],'ytick',[0:500:yl(2)],...
                'xticklabel',[0:100:(xl(2)/pix_per_chan)],...
                'yticklabel',[0:100:(yl(2)/pix_per_chan)])
            xlabel('Down valley distance (channel widths)')
            ylabel('Cross valley distance (channel widths)')
            title(sprintf('Topography at time step %d ' ,i))
            cmapsPUB('topo',0);
            h= colorbar;
                caxis([ 0 1+.5/length(colormap)])
                set(h,'Location','southoutside')
                ylabel(h,'Normalized Floodplain Thickness')
            ht = text(980,500,'Stratigraphic Cross section',...
                    'HorizontalAlignment','center','Rotation',90);
            text(100,1720,'5 km scale bar')
            axis equal
            drawnow

%       stratigraphic cross section figure
        d_elev=(squeeze( stratigraphy( :,62,: )));  % choose one slice
%         d_elev = d_elev(99:351,:);  % not all the way accross the data
        % more dynamic method, fewer errors:
        de_hi = min(351,size(d_elev,1)); de_lo = max(de_hi-252 , 1 );
        d_elev = d_elev(de_lo:de_hi,:);  % not all the way accross the data
        
        tmp = sum(d_elev);
        beginning = find(tmp>0,1,'first');
        ending    = find(tmp>0,1,'last');
        if isempty(beginning) beginning = 1; end % catch empty values
        if isempty( ending) ending = size(tmp,2); end
            
        d_elev = (d_elev(:,beginning:ending));
        pt_bar_top = max(d_elev'); % collect the top of the point bar

        if size(d_elev,2)>200 % array is too big for stacked bar chart
            % 'stratigraphy' image in time and space, color is thickness
            hxsec = figure(51); 
                imagesc(flipud(log(d_elev')))
                colormap(jet)
                h=colorbar;
                set(h,'Direction','reverse')
                ylabel(h,strcat('log of deposit thickness from each',...
                    ' time step, pt bar is red'))
                title(strcat('Pseudo stratigraphy, large arrays will',...
                    ' not plot as stacked bar charts'))
                xlabel('Cross section width (pixels)')
                ylabel('Time steps since deposit created')
        else                % make a stacked bar chart
        % for unknown reasons, this bar chart plots with overlapping bars
        % and the figure for publication was made using a stacked bar chart
        % in Microsoft Excel
%       chronostratigraphy of selected cross section
        hxsec = figure(5);  % stratigraphy colored by age
        clf
            set(hxsec,'position',[400 300 660 420],'color','w')
            bar(([1:size(d_elev,1)]),(d_elev),'stacked','barwidth',1,...
                'edgecolor','none') %flipud switches left to right
            xlim([.5 size(d_elev,1)+.5])
            ylim([0 1])
            xlabel('Horizontal Distance (channel widths)')
            ylabel('Normalized Floodplain Thickness')
            title('Stratigraphic Cross Section',...
                'fontsize',14,'color',[197 90 17]/255)
            h=colorbar;
                cdata=cmapsPUB('xsec',1);
                cdata = cdata(3:end,:);
                if length(cdata)>size(d_elev,2)
                    colormap(cdata(1:size(d_elev,2),:))
                else
                    cdata=cmapsPUB('age',0);
                    rng = size(d_elev,2);
                    colormap(cdata(end-rng:end,:))
                end
                nos = get(h,'xtick');
                set(h,'xtick',sort(h.Limits(2)-nos),...
                    'xticklabel',fliplr(time_step_years*(nos+i-ending)))
            ylabel(h,sprintf('Deposit Age (years). Time step %d',i))
            set(hxsec.CurrentAxes,'xtick',[.5:50:255])
            set(hxsec.CurrentAxes,'XTickLabel',...
                num2cell((str2num(strjoin(get(hxsec.CurrentAxes,...
                'XTickLabel')))-.5)/5),'fontsize',14)
            hold on % add dots at the top of the point bar
                plot(pt_bar_top,'color','k','linestyle','none',...
                    'marker','.','markersize',1)
            hold off
            xlim( [ .5 250.5])
            drawnow
        end
%       Channel length (in pixels) mass & accommodation space (full pixels)
        hmass = figure(6); 
            set(hmass,'color','w','position',[40 100 1000 400]);
            subplot(2,1,1)
                plot(chan_area(1:i)/pix_per_chan)
                title('Channel length in pixels')
            subplot(2,1,2)
                px = 1:i;
                ha=plotyy(px,mass(1:i),px,max(mass)-mass(1:i));
                title('Mass in storage & accomodation space')
                legend('mass','accomodation')
                ylabel(ha(1),'Full pixels') %empty pixels
                ylabel(ha(2),'Full pixels')
            drawnow
%       Migration rate measured in pixels
        hmrat = figure(7);
        % recreate handle hmig in case the figure was deleted or
        % overwritten can happen if user clicks on or closes figures
        hmig = histogram( old_fp_distance( bwperim( ~fp_old_and_new ) ),...
                                                            mig_bins );
            BinWidth = hmig.BinWidth;
            BinEdges = hmig.BinEdges(1:end-1); 
            set(hmrat,'color','w','position',[40 300 1000 400]);
            % plotting by the smaller bin edge is appropriate for small
            % migration rates which are 0, 1, sqrt(2), etc so the first
            % frequency bar is for the migration rate of 0, the second
            % frequency bar has migration rates of 1 & sqrt(2), etc
            % hf = bar( hmig.BinEdges(1:end-1),hmig.Values ); % does not
            % work because hmig is a handle to a histogram figure which
            % can easily be overwritten or deleted
            hf = bar( BinEdges,mig_distance );
            title(sprintf('Migration distance. Bin width is %d',BinWidth))
            ylabel( 'Frequency')
            xlabel(sprintf(['Distance in pixels. %dm/pixel, channel',...
                            'width is %dm, Time step is %d years'],...
                            pixel_sz, 2*B, time_step_years))
%           This! is how axes should be adjusted, handles and 'children'                        
            hmrat.CurrentAxes.XLim = [min(mig_bins)-.5 max(mig_bins)+.5];
            drawnow


%% save output image files, the user is encouraged to save final figures as
%   .svg or .fig formats for easier editing
    if save_images
        ts = sprintf('Time step %d.png',i);
        tmp = round(disp_progress * dt_save_years * dt / 10 ) * 10 ;
        outfile=sprintf('%d year channel presence map. %s', tmp,ts);
        print(hocc,'-painters', '-dpng', '-r600',outfile)
        outfile=sprintf('Map of floodplain age. %s',ts);
        print(hage,'-painters', '-dpng', '-r600',outfile)
        outfile=sprintf('Map of floodplain topography. %s',ts);
        print(htopo,'-painters', '-dpng', '-r600',outfile)
        outfile=sprintf('Selected cross section. %s',ts);
        print(hxsec,'-painters', '-dpng', '-r600',outfile)
        outfile=sprintf('Mass, Accomodation, & Channel length. %s',ts);
        print(hmass,'-painters', '-dpng', '-r600',outfile)
        outfile=sprintf('Areas of various floodplain segments. %s',ts);
        print(harea,'-painters', '-dpng', '-r600',outfile)
        outfile=sprintf('Migration rate measured in pixels log. %s',ts);
        print(hmrat,'-painters', '-dpng', '-r600',outfile)
        % copy variable name, backward compatible for other code:
        chan_occupation = chan_occ_one; 
        save( sprintf( 'channel_occupation_%sk_%syr_%s',...
            num2str(sim_time_ky), num2str(time_step_years), ...
            num2str(i) ), 'chan_occupation' ,'-v7.3')
    end
%  reset chan_occ_one to count channel position until the next display stop  
    chan_occ_one = false(size(chan_occ_one));

% status update on program run time
    elapsed_algorithm_time = datenum(clock) - datenum(algorithm_start);
    est_complete = datenum( algorithm_start ) +...
        elapsed_algorithm_time / frac_complete ;
    disp(sprintf(strcat('Figures & saving complete \nCurrent Time:  %s',...
        '\nest completion %s'), datestr(clock) , datestr( est_complete)))


     end %dispaly stops
end %i
% save map of the entire floodplain ever visited (occupied) by the channel
outfile = sprintf('Previously Occupied by Channel. Time step %d.mat',i);
save(outfile,'prev_channel','-v7.3')
% save other data
outfile = sprintf('other_data_3chan_%sk_%syr', num2str( sim_time_ky ), ...
                                            num2str( time_step_years ) ) ;
save(outfile ,'fp_area', 'mass', 'chan_area', 'mig_distance',...
    'x_start', 'x_end', 'ymn', 'ymx', 'pix_per_chan','d_elev',...
    'new_elev','bool_xsec','birth_step','chan_occ_all','-v7.3')

algorithm_end =clock
agorithm_duration = algorithm_end - algorithm_start
%%
disp('You have reached the end of display_model.m')
disp('You have keyboard control of the variables.')
disp('You are encouraged to save variables and figures as desired')
disp('.svg and .fig are reccommended for editing.')
disp('Type dbcont or dbquit to end keyboard control')
keyboard % give user control dbquit or dbcontinue

end % function
