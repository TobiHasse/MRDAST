function deposit_storage_time(riv, riv2, B, show_figures, restart,...
    x_starts, x_ends, ymn, ymx, slice, ...
    pt_bar_dist, vert_dist, vd_age_dist, pb_age_dist )
% Purpose:  This function will take the meandering river simulation and
%           compute the storage time and age distributions for point bar
%           and vertically accreted overbank flood deposits.  
%           The function is designed to use a subset of the simulation
%           space defined by x and y limits and the end simulation time.
%           The model space is analyzed in small reaches called slices to
%           manage RAM requirements.
%           The modeler is encouraged to run the function display_model.m
%           before this function to help select the x, y, and time, limits
%           to analyze.
%           output data is saved to file during the simulation to enable
%           restarting the simulation if there is a crash or power outage
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     edited June 2021, written 2015

fprintf('Welcome to: "%s" \n',...
        fullfile(fileparts(mfilename('fullpath')),mfilename))

load params_storage.mat
load params_deposition.mat
end_sim = min(end_sim,numel(riv))
% the first slice to analyze, change this if restarting after crash
first_k = slice + 1; % start at 1 if not restart 
program_start = clock

% display model progress increment
disp_progress = round( display_progress_years / time_step_years )  
% disp_progress = 10            % DEBUG option display progress frequently

%% constants and user defined inputs
display_stops = [ disp_progress : disp_progress : end_sim, end_sim ];

%% find x y limits of simulation area
limits = domain_limits(riv,riv2,pix_per_chan,B);

%% if restarting the run, certain variables should not be initialized
if ~restart 

% make slices to analyze part of the simulation at a time (memory limits)
    if load_limits
        try
        % this file contains the extents of each slice for the dissertation
        % 211ky_xy_slice_box.mat is available at 
        % https://doi.org/10.5281/ZENODO.5651874. 
            infile = '211ky_xy_slice_box.mat';
            load(infile)  % load new x_steps and yrange
            disp( ['    x_starts','    x_ends','        ymn','       ymx'])
            [x_starts',x_ends',ymn',ymx'] % show the values
            
        catch %else
            sprintf(['Unable to load subreach extents from file: %s\n'...
                'calculating new subreach extents'],infile)
        end
    else
        show_slice_figures = false;
        % calculate the extents of the subreaches (slices)
        [x_starts, x_ends, ymn, ymx] = slices( array_size, end_sim, ...
            x_start, x_end, show_slice_figures);
    end
    
    if x_starts(1) > x_ends(end)
        fprintf(['Error inside function:\n%s\nThere is no region to', ...
            ' analyze \nX start: %d, X end: %d'],...
            fullfile(fileparts(mfilename('fullpath')),mfilename),...
            x_starts(1),x_ends(end))
        return % break % break for script, return for function (perhaps)
    end

    %% make some arrays for storing variables
    % if the program crashes (power outage, etc) it is possible to restart
    % skipping completed slices and load in appropriate variables
    pt_bar_dist    = zeros( end_sim ); % storage time distributions
    vert_dist      = zeros( end_sim );
    %     pbar_e_haz     = zeros( end_sim ); % erosion hazard distributions
    %     vert_e_haz     = zeros( end_sim );
    vd_age_dist    = zeros( end_sim ); % floodplain age distributions 
    pb_age_dist    = zeros( end_sim );
    
    algorithm_start = clock
    j = 0;   % counter for ticker times
    ticker_times = zeros(numel(x_starts)*numel(display_stops),13);
    t_new_isl = 0; t_rep = 0; % variables for tracking computation times
    t_perm = 0; t_clr = 0; t_dep = 0; t_reset = 0; t_clean = 0;

else % do this if restarting after slice# completed slices
    j_for_one_slice = round( end_sim / disp_progress );
    % increment for ticker times: steps per slice
    j = j_for_one_slice * (first_k - 1) - 1; 
    load tickers.mat
    temp = num2cell(ticker_times(j,:));
    [t_new_isl, t_rep, t_perm, t_clr, t_dep, t_reset, t_clean, ...
        ~,~,~,~,~,~] = deal(temp{:});
    % show user start time
    algorithm_start % should be loaded from tickers.mat on restart 
end

%% cycle through each slice, starting with the first_k^th slice (if 
% restarting the model analysis)

for k = first_k:length(x_starts)
    x_step = x_ends(k) - x_starts(k);
%% Make boolean array of channel location
% Initialize the array (needed outside of loop one time)
% doIntermp_TRH edited MATLAB function to return desired values
[ x, y ]     = doInterpm_TRH( riv( 1 ).Xcl, riv( 1 ).Ycl,...
                    pixel_sz/interp_nds_pix, 'lin' ); 
% bool_channel is a raster the size of the model domain, channel = true
bool_channel = node2grid( x, y, pixel_sz, limits); 

% Make some more arrays
vrt_accum      = zeros( ymx(k)-ymn(k)+1 , x_step, end_sim );
lat_accum      = zeros( ymx(k)-ymn(k)+1 , x_step, end_sim );
accum          = zeros( ymx(k)-ymn(k)+1 , x_step );
old_elev       = zeros( size( accum ) );
new_elev       = zeros( size( accum ) );
bool_pt_bar    = false( size( accum ) );
fp_old_and_new = false( size( accum ) );
bool_erosion   = false( size( vrt_accum ) );

%% Distance transform of channel offset so 0 is the edge of channel pixels
old_fp_distance = double( bwdist( bool_channel ) - ...
                    floor( pix_per_chan / 2) );  
old_fp_distance( old_fp_distance < 0 ) = 0; % make all channel pixels = 0
% cut down to extents of the reach (slice), remove overlap on x axis (-1)
old_fpd_slice = old_fp_distance( ymn(k):ymx(k) , x_starts(k):x_ends(k)-1 );  

% Start floodplain entirely full so change in elevation = 0 at beginning
old_elev( old_fpd_slice > 0 ) = 1;
% pb_age_dist(1) = 0; % not sure what these 2 lines are for
% vd_age_dist(1) = 0; % [ deposit fp_age_dist(i-1)-storage_dist(i) ];
% innitialize quantities of sediment in the age distributions
vd_age_curr = 0;    % starts empty
pb_age_curr = numel(accum);% starts full, 'terrace' assigned to point bar

lat_accum( :,:,1 ) = old_elev(:,:) ; % store the initial 'terrace' surface

tic

for i = 2 :1:  end_sim % size( riv, 2 ) % 6845
    % channel node spacing might be larger than pixels, interpolate more
    [ x, y ]        = doInterpm_TRH( riv( i ).Xcl, riv( i ).Ycl,...
                        pixel_sz/interp_nds_pix, 'lin' );
    % take the x y points and create pixel map with the channel = ture
    bool_channel    = node2grid( x, y, pixel_sz, limits);
    % Distance from channel edge offset so 0 is the edge of channel pixels
    new_fp_distance = double( bwdist( bool_channel ) - ...
                        floor( pix_per_chan / 2) );  
    new_fp_distance( new_fp_distance < 0 ) = 0; % make all chan pixels = 0
    % cut down to extents of the reach (slice)
    new_fpd_slice = new_fp_distance( ymn(k):ymx(k) , ...
                                            x_starts(k):x_ends(k) -1 );  

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
    % cut down to extents of the reach (slice)
    fp_oan_slice = fp_old_and_new( ymn(k):ymx(k),x_starts(k):x_ends(k) -1);

    t_new_isl = t_new_isl + toc; tic                % computation timers

    % compute storage times separately for point bar and vertically 
    % accreted deposits.  To do this we need to:
    % make a 3D array of recently eroded or channel locations
    bool_erosion = repmat( ~fp_oan_slice, [ 1,1,size( lat_accum, 3 ) ] );
    t_rep = t_rep + toc; tic                        % computation timers

    % Note: it is possible that the storage time distributions could be 
    % calculated after this function completed based on the age 
    % distributions (computed later about 100 lines below in this code).
    % The age time distributions (via the implicit method) are calculated
    % based on the storage time distributions calculated here.
    % This would reducing calculations here and free up memory.  
    % This code is left as is as it represents how the storage time 
    % distributions were computed for Tobias Hasse's dissertation Chap 2
    
    % to calculate storage distributions based on the age distributions:
    % storage = zeros(size(age));
    % storage(2:end,1:end-1) = age(1:end-1,2:end)-age(2:end,1:end-1);
    
    % The simple way to select the eroded material doesnt' work:
    % ARRAY(Bolean_array) becomes linearized, ruining the age structure.
    % Other methods: Indexing into the arrays is faster where it cuts out a 
    % lot of the array.  Possibly anonymous functions with structs, or cell
    % arrays and a 2D bool_erosion.  This could cut down on the amount of 
    % 3D stuff and free up memory.... these optimizations were not persued
    % and simple .* element multiplication is used

    if i < .7 * end_sim %3500   % speed optimization takes ? i * .004s 
        % in the sum function, think of a book of graph paper, the order is
        % the same as the output of size()
        % sum(sum(sum... vertical, horizontal, page... 
        pt_bar_curr = permute( sum( sum( lat_accum   (:,:,1:i ) .* ...
            bool_erosion(:,:,1:i ) ) ), [ 3 2 1 ] ); % from 1:i is all 
            % stored, from 2:i gets rid of the terrace (I do that later)
        vert_curr   = permute( sum( sum( vrt_accum   (:,:,1:i ) .* ...
            bool_erosion(:,:,1:i ) ) ), [ 3 2 1 ] );
    else          % else always takes ? 8s (based on 100ka sim, dt = 20yr)
        pt_bar_curr = permute( sum( sum( lat_accum .* bool_erosion ) ), ...
            [ 3 2 1 ] ); % this method can't get rid of terrace (quickly)
        vert_curr   = permute( sum( sum( vrt_accum .* bool_erosion ) ), ...
            [ 3 2 1 ] );
        pt_bar_curr = pt_bar_curr(1:i);
        vert_curr = vert_curr(1:i);
    end
    % combine current distributions to previously analyzed distributions +=
    pt_bar_dist( i , end-i+1:end ) = pt_bar_dist( i , end-i+1:end ) + ...
        pt_bar_curr';  %+= because analyzing in slices
    vert_dist  ( i , end-i+1:end ) = vert_dist  ( i , end-i+1:end ) + ...
        vert_curr'; 
    
%%%%%%%%%%%%%%%%%%% HOW TO STRIP OUT TERRACE ******************************
% for i = 1:size(pt_bar_dist,1) % (only do this once at the end)        %**
%     id=find(pt_bar_dist(i,:),1,'first');                              %**
%     pt_bar_dist(i,id)=0;                                              %**
% end                                                                   %**
%              OR STRIP THE TERRACE IN ONE LINE OF CODE                 %**
% % get rid of the 'terrace' initial floodplain height                  %**
%     pb_dist = fliplr(tril(fliplr(pt_bar_dist),-2));                   %**
%%%%%%%%%%____________________________________________________*************    

    t_perm = t_perm + toc; tic                      % computation timers
    
    if i < .55 * end_sim %2700% delete eroded data (3x faster than permute)
        vrt_accum( bool_erosion( :,:,1:i ) ) = 0; % delete the eroded data
        lat_accum( bool_erosion( :,:,1:i ) ) = 0; 
    else  % slow and steady( faster for large arrays )
        vrt_accum( bool_erosion ) = 0; % delete the eroded data
        lat_accum( bool_erosion ) = 0; 
    end
    old_elev( ~fp_oan_slice ) = 0; % delete the eroded data

    t_clr = t_clr + toc; tic                        % computation timers
  
    % determine point bar locations (this includes cutoff channels which
    % may aggrade quickly before cutoff)
    bool_pt_bar( ~ (fp_oan_slice | new_fpd_slice == 0 ) ) = 1;       
    %     figure(9)
    %         imagesc(bool_pt_bar)
    %         axis equal
    %         colorbar
    %         drawnow

    % Vertical accretion based on equation 1 from Chapter 2 of Tobias
    % Hasse's dissertation, adapted from Alan Howard 1992 Lowland 
    % Floodplain Rivers Modeling Channel Migration and Floodplain 
    % Sedimentation in Meandering Streams (eq 1.17)
    % I considered including the new point bar elevation with the old
    % elevation for computing the remaining accomomdation space which would
    % have reduced the topographic variability on the pseudo scroll bars.
    % old_elev = 0 for new point bars. I kept the 'big' scroll bars

    % calculate the thickness of newly accumulated sediment: pseudocode:
    % newly formed point bars +
    % remaining acommodation space with exponential decay (vertical)
    % pelagic + diffusive sedimentation * expon decay (horizontal distance)
    accum = bool_pt_bar * pt_bar_elev + ...      
        ( 1-old_elev ) ./ exp( old_elev ) .* ... 
        (Nu + Mu * exp( -new_fpd_slice / ( 2*B/pixel_sz * lambda ) ) );
    
    % compute new elevation and change in elevation
    new_elev = old_elev .*(1-bool_pt_bar) + accum;     % find new elevation
    new_elev( new_fpd_slice == 0 ) = 0 ;  % cut new channel b.c. exp(0) = 1

    % store changes in elevation
    vrt_accum( :,:,i ) = new_elev .* ~bool_pt_bar - old_elev;  
    lat_accum( :,:,i ) = new_elev .*  bool_pt_bar;  
    
    t_dep = t_dep + toc; tic                        % computation timers
    
    % compute age distributions and erosion hazard

    % Vertical accumulation, position in the array determins age:
    %     ii=1: nothing deposited - no accommodation space
    %     ii=2: only pt bar deposited
    %     ii=3: overbank deposits, but no erosion
    %     ii=4: overbank deposits can be eroded

    if implicit_age_dist 
        % compute the current age distribution based on storage time
        % distributions above
        % pseudocode: current = [ previous - removed, newly deposited]
        vd_age_curr = [ vd_age_curr - vert_curr( 1:end-1 )', ...
            sum( sum( vrt_accum( :,:,i ) ) ) ];
        vd_age_curr( vd_age_curr < 10^ -10 ) = 0; % erase tiny amounts
        % store the age distribution
        vd_age_dist( i , end-i+1:end ) = vd_age_dist( i , end-i+1:end )+...
                                                        vd_age_curr ;
        % repeat for the point bar
        pb_age_curr = [ pb_age_curr - pt_bar_curr( 1:end-1 )', ...
            sum( sum( lat_accum( :,:,i ) ) ) ];
        pb_age_curr( pb_age_curr < 10^ -10 ) = 0 ;
        pb_age_dist( i , end-i+1:end ) = pb_age_dist( i , end-i+1:end )+...
                                                        pb_age_curr ;
    else % find the volume of sediment remaining at each time step
        % pseudocode: age distribution from previous slices at this time
        % step + age distribution from current slice at this time step
        % TEST this method: it might be faster than computing the storage
        % times 
        vd_age_dist( i , end-i+1:end ) = vd_age_dist( i , end-i+1:end )+...
            permute( sum( sum( vrt_accum( :,:,1:i ) ) ), [ 3 2 1 ] )';
        pb_age_dist( i , end-i+1:end ) = pb_age_dist( i , end-i+1:end )+...
            permute( sum( sum( lat_accum( :,:,1:i ) ) ), [ 3 2 1 ] )';
    end
    
    % Compute erosion hazard functions (optional)
    %   vert_e_haz( i-1 , 1:end-1 ) = vert_dist  ( i , 1:end-1 ) ./ ...
    %       vd_age_dist( i-1 , 2:end ) ;
    %   pbar_e_haz( i-1 , 1:end-1 ) = pt_bar_dist( i , 1:end-1 ) ./ ...
    %       pb_age_dist( i-1 , 2:end ) ;
    
    % reset arrays for next iteration
    bool_pt_bar     = bool_pt_bar * 0;
    old_fp_distance = new_fp_distance;
    old_elev        = new_elev;

    t_reset = t_reset + toc; tic                    % computation timers

%% plot progress incrementally or at the end, save output
    if any( display_stops == i ) 
        pct = i / display_stops(end)*100;
        frac_complete = pct / 100 / numel( x_starts ) + ( k-1 ) / ...
                                                    numel( x_starts );
        fprintf('Starting figures, \nCurrent time:  %s', ...
            datestr(datenum(clock)))

%% Make a figure showing current floodplain topography
        figure(1)
            clf
            set(figure(1),'color','w');
            imagesc( new_elev )
            title(sprintf('Topography at time step %d slice #%d' ,i,k))
            xlabel(sprintf('Slice from %d to %d pixels',x_starts(k),x_ends(k)))
            ylabel(sprintf('Slice from %d to %d pixels',ymn(k),ymx(k)))
            cmapsPUB('topo',0);
            colorbar
            axis equal
            drawnow
%% (optional) store the simulation computation time variables
        j=j+1;
        ticker_times(j,:) = [t_new_isl, t_rep, t_perm,t_clr,t_dep, ...
            t_reset, t_clean, frac_complete * 100 , pct, k, ...
            numel( x_starts ), nnz( vrt_accum )/ numel( vrt_accum ),...
            nnz( lat_accum )/ numel( lat_accum )];
        ticker_variables = sprintf(['  t_new_isl,  t_rep,    t_perm, ',...
            '   t_clr,    t_dep, t_reset,    t_clean, %%fin total,',...
            '%%fin slice,slice, of slices, vsparsity, lsparsity'])
        disp(ticker_times(j,:))
        ticker_plot(ticker_times,ticker_variables)
% Saving ticker timers moved to end of slices since many slices calculated
%     save('tickers.mat','ticker_times','ticker_variables')

%% Update progress
        elapsed_algorithm_time = datenum(clock) - datenum(algorithm_start);
        fprintf(['Program is %s%% complete. (%s%% finished ',...
            'calculating loop %d of %d)\n'],...
            num2str( frac_complete*100 ),num2str(pct), k, numel(x_starts))
        fprintf(['Figures & saving complete \nCurrent Time:  %s \nest ',...
            'completion %s\n'], datestr(clock) ,...
            datestr( datenum( algorithm_start ) + ...
            elapsed_algorithm_time / frac_complete ) ) 

    end %dispaly stops
    t_clean = t_clean + toc; tic                    % computation timers
end %i
%% at end of slice, save and show figures
slice = k;
outfile = sprintf('distr_3chan_%sk_%syr_%s', num2str( sim_time_ky ), ...
                                num2str( time_step_years ), num2str( k ) );
save(outfile, 'pb_age_dist', 'pt_bar_dist', 'vd_age_dist', 'vert_dist', ...
    'x_starts', 'x_ends', 'ymn', 'ymx', 'slice', 'end_sim', ...
    'array_size', '-v7.3')
save('tickers.mat','ticker_times','ticker_variables','algorithm_start')

fprintf('Slice variables saved. Elapsed time %0.2f minutes\n',toc/60)
%% Make figures of the distributions through the most recent slice 
    if show_figures
        tic
        ts = sprintf('Slice: %d',k);

        figure(20) % point bar distribution
            imagesc(fliplr(pt_bar_dist))
            title(sprintf('Point bar storage time. With terrace. %s',ts))
            cmapsPUB('topo',1);
            colorbar
            drawnow
        figure(21)
            % get rid of the 'terrace' initial floodplain height
            pb_dist = fliplr(tril(fliplr(pt_bar_dist),-2));  
            imagesc(fliplr(pb_dist))
            title(sprintf(['Point bar storage time.', ...
                                            ' Without terrace. %s'],ts))
            cmapsPUB('topo',1);
            colorbar
            drawnow
        figure(22)
            imagesc(fliplr(vert_dist))
            title(sprintf(['Overbank (vertically acreted) storage ',...
                                                        'time. %s'],ts))
            cmapsPUB('topo',1);
            colorbar
            drawnow
        figure(23)
            pb_dist = fliplr(tril(fliplr(pb_age_dist),-2));  
            imagesc(fliplr(pb_dist))
            title(sprintf(['Age distribution of point bar material. %s',...
                                            ' Without terrace. %s'],ts))
            cmapsPUB('topo',1);
            colorbar
            drawnow
        figure(24)
            imagesc(fliplr(pb_age_dist))
            title(sprintf(['Age distribution of point bar material. ',...
                                            ' Without terrace. %s'],ts))
            cmapsPUB('topo',1);
            colorbar
            drawnow
        figure(25)
            imagesc(fliplr(vd_age_dist))
            title(sprintf(['Age distribution of overbank (vertically ',...
                                        'acreted) material. %s'],ts))
            cmapsPUB('topo',1);
            colorbar
            drawnow
        figure(26) % Erosion Hazards, see Bradley & Tucker 2013
            imagesc(fliplr((vert_dist(2:end,1:end-1)./...
                vd_age_dist(1:end-1,2:end))))
            title(sprintf(['Erosion hazard of overbank (vertically ',...
                                        'acreted) material. %s'],ts))
            cmapsPUB('topo',1);
            colorbar
            drawnow
        figure(27)
            imagesc(fliplr((pt_bar_dist(2:end,1:end-1)./...
                pb_age_dist(1:end-1,2:end))))
            title(sprintf(['Erosion hazard of point bar material. %s'],ts))
            cmapsPUB('topo',1);
            colorbar
            drawnow
        % to make the minimum at the threshold between the first 2 colors
        %             ca = caxis;
        %             caxis([-ca(2)/(size(colordata,1)-1) ca(2)])
        %%
        % for ii=20:27 % DEBUG set axis limits if inspecting short run
        %     figure(ii)
        %     xlim([0 300])
        %     ylim([0 300])
        % end
        
        fprintf(['Distribution figures complete. ',...
            'Elapsed time %0.2f minutes\n'], toc/60 )
        tic
    end % if show_figures


end %k finished all the slices
toc
algorithm_end =clock
agorithm_duration = algorithm_end - algorithm_start

[t_new_isl, t_rep,     t_perm,   t_clr,     t_dep,   t_reset,  t_clean]
disp('t_new_isl, t_rep,     t_perm,   t_clr,     t_dep, t_reset,  t_clean')

disp('You have reached the end of deposit_storage_time.m')
disp('You have keyboard control of the variables.')
disp('You are encouraged to save variables and figures as desired')
disp('.svg and .fig are reccommended for editing.')
disp('Type dbcont or dbquit to end keyboard control')
keyboard % give user control dbquit or dbcontinue

end % function