function [x_starts, x_ends, ymn, ymx] = slices( array_size, end_sim, ...
    x_start, x_end, show_figures)
% Purpose:  This function finds the minimum bounding box for each slice of
%           the model based on a sillouette of the alluvial valley.  This
%           is required to make the RAM requirements manageble
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     June 2021


x_starts(1) = x_start;
x_ends=[]; ymn=[]; ymx=[];
array_sz_2d = floor(array_size/end_sim);

no_infile = true;
try
    infile = 'Previously Occupied by Channel. Time step 6845.mat';
    load(infile)
    no_infile = false;
catch
    fprintf('Error loading: "%s" inside function:\n%s\n',...
        infile,fullfile(fileparts(mfilename('fullpath')),mfilename))
    fprintf(['Make sure a map of the floodplain visited by the'...
        ' channel exists\n\n Attempting to create the map \n\n'])
    try
        [~,prev_channel]= chan_occ('Chan_occ_map.gif',inf);
        disp('Success!!')
    catch
        disp('Attempt failed, function aborted before completion')
    return
    end
end

% save the infile for next time to save time
if no_infile
    save(infile, 'prev_channel')
end

% Find the first and last y value in each column of the matrix
% When multiple values are tied for the max, max returns the first index 
[~,yf]=max(prev_channel,[],1);  
[~,yl]=max(flipud(prev_channel),[],1);
yl = size(prev_channel,1)-yl; % flip the indexes

%% OPTIONAL plot the map of visited floodplain and the perimeter limits
if show_figures
    px=1:length(yf); % vector for plotting output
    figure
        h = subplot(2,1,1);
            plot(px,yf,px,yl);
            xlim([0 size(prev_channel,2) ])
            set(h,'yDir','reverse')
        subplot(2,1,2)
            imagesc(prev_channel)
        drawnow
end

i = 1;
while x_starts(i) < x_end
    slice_size = 0;
    j = 0;
    while slice_size < array_sz_2d
        y_max = max( yl ( x_starts(i):x_starts(i)+j ) );
        y_min = min( yf ( x_starts(i):x_starts(i)+j ) );
        slice_size = (j-1) * ( y_max - y_min );
        j=j+1;
    end
    i=i+1;
    % record the boundaries of the slice into vectors
    % since we don't know a priori how big these arrays are, keep making
    % them bigger each time, memory inefficient, but small arrays
    ymn(i-1)    = y_min;
    ymx(i-1)    = y_max;
    x_starts(i) = x_starts(i-1) + j-1 ;
    
end

x_starts(x_starts >= x_end) = [];   % remove any starts after the end
x_ends = [ x_starts(2:end) x_end ]; % create vector for ends of the slice
% note: x_ends could have a -1 correction here to avoid overlapping because
% the slice is inclusive of the ends, but the -1 correction is done
% within the deposit_storage_time.m code

%% display map of slices, array size, slice width 
if show_figures
    slice_map = zeros( size( prev_channel ) );
    for i = 1:2:numel(x_starts)
        try % try since it is likely that i+1 will exceed numel of x_starts
            slice_map(ymn(i):ymx(i),x_starts(i):x_ends(i))         = .3;
            slice_area( i ) =   (    ymx( i ) -      ymn( i ) ) * ...
                                ( x_ends( i ) - x_starts( i ) ); 
            slice_map(ymn(i+1):ymx(i+1),x_starts(i+1):x_ends(i+1)) = .7;
            slice_area(i+1) =   (    ymx(i+1) -      ymn(i+1) ) * ...
                                ( x_ends(i+1) - x_starts(i+1) );
        catch
        end
    end
    %
    figure(33)
        imagesc(slice_map + prev_channel)
        title('Map of the slices to be analyzed for storage time')
        axis equal
        drawnow
    figure(34)
        px = 1:numel(x_ends);
        hp = plotyy( px,slice_area*end_sim , px,x_ends-x_starts );
        title('Number of elements in each array of sediment')
        ylabel(hp(1),'Elements in the array')    
        ylabel(hp(2),'Width of slice (pixels)')
        xlabel(hp(1),'Slice number')
        drawnow
end % show figures
        
end % function
