function [ not_channel ] = rm_islands9( new_fp, mid_fp, old_fp, min_px_at_cutoff, channel_w_in_pix )
% function [ not_channel ] = rm_islands9( new_fp, mid_fp, old_fp, min_px_at_cutoff, region_index, channel_w_in_pix )
% this is a major overhaul of rm_islands February 2016 by Tobias Ackerman
% it should be much more efficient and uses new methods including analyzing
% the grayscale image
% this is modified from rm_islands873.m 
% and excludes the manual island identification and figure handling from rm_islands7.m
% and includes 3 channel inputs (one at the halfway point) for better
% isalnd identification

% including a partial dt channel position means nearly all small areas
% should be removed, and only the big ones with high change in distance are
% cutoffs and may need strings removed

% INPUTS
%     new_fp  - 2D array, distance transform from current channel position
%     old_fp  - 2D array, distance transform from previous channel position
%     mid_fp  - 2D array, distance transform from more recent than previous channel position
%     min_px. - Scalar, The smallest likely neck cutoff ( in pixels )

%     SCALARS:  ( not currently passed )
%     neck_dist_thresh     - min value to id neck cutoff
%     crescent_dist_thresh - max value to id crescents
%     crescent_expansion   - scalar on min_px. to id big crescents
%     channel_w_in_pix     - for use with fat_channel to avoid resolution issues

% OUTPUTS   
%     not_channel - 2D array, Boolean of areas that were not channel through both time steps

neck_dist_thresh = 8;  % 8.6 .. 13.3
% crescent_dist_thresh = 0.5;
% crescent_expansion = 3;  % regions x * bigger than min_pix_@_cutoff may be crescents
% channel_w_in_pix = 3;
% remove small islands
% not_channel = bwareaopen( old_fp & new_fp & mid_fp, min_px_at_cutoff/5 );  %DEBUG /5 shows smaller regions in case they should be kept
not_channel = bwareaopen( old_fp & new_fp & mid_fp, min_px_at_cutoff ); 

% NOTE: it does not work to find cutoffs before and after separately
% because the cutoff can happen between the before and after observations.
% When this happens, the cutoff is not observed without looking at both
% channel positions simultaneously.... back to the statistics

% This threshold should be the minimum area of floodplain that is cut off
% measured in pixels.  From Schwenk 2015 Area at cutoff/Lo^2 = 8 to 100
% Lo = ±4B, with three pixel wide channels Lo = ±6px
% ********** NOTE ********************July 2016********************** %
% since Schwenk meanders are 3x too big, this is likely wrong         %
% ******************************************************************* %
% Therefore the minimum pixels (area) at cutoff = 6^2 * 8 = 288
% min_pix_at_cutoff = 200;

% for plotting out figure for manual check
% one_region   = false(size(not_channel));

% change in distance from channel
del_distance = new_fp - old_fp + 100;        % + 100 b.c. must be > 0 for bwlabel
% del_distance = new_fp - old_fp;
% del_distance_shift = min(del_distance(:))+1
% del_distance = del_distance + del_distance_shift;
del_distance( ~not_channel ) = 0;            % burn channel into image
L = bwlabel(del_distance);                   % identify regions

% remove biggest region from regionprops calculation to save time
L(L==1)=0;
del_distance = del_distance - 100;           % fix del_distance

% find the properties of the grayscale regions
area_props = regionprops(L,del_distance,'Area','MeanIntensity','PixelIdxList');
% % DEBUG redo removal of small areas when min_px_at_cutoff is reduced above
% not_channel = bwareaopen( old_fp & new_fp & mid_fp, min_px_at_cutoff );  
% rm_islands9_L69 = numel(area_props)
% categorize regions as neck cutoffs or big crescents
% NOTE small regions have already been eroded, and large regions that are
% not neck cutoffs don't need anything done to them
neck_cutoffs = find( [ area_props.MeanIntensity ] > neck_dist_thresh );
%                                                  * 1 if dt = 10    * 3 if dt = 20
% big_crescents = find( [ area_props.Area         ] < min_px_at_cutoff * crescent_expansion & ...
%                       [ area_props.MeanIntensity] < crescent_dist_thresh );

% remove big crescents
% not_channel( vertcat( area_props( [ big_crescents ] ).PixelIdxList ) ) = 0;

% remove strings from extant neck cutoffs, skip all other regions
for i = 1:numel(neck_cutoffs) 
%     i_is_rm8_L107 = i
%     numel(neck_cutoffs)
    % for steps with multiple cutoffs, it would be faster to do fat new
    % channel outside of loop, but if there is no neck cutoff, fat new
    % channel is not created at all
    fat_new_channel = false( size( not_channel ) );
%     fat_new_channel( new_fp < (channel_w_in_pix * 2) ) = 1;  % new_fp < 6
    fat_new_channel( new_fp < (channel_w_in_pix + 2) ) = 1;
    balloon_new  = false( size( not_channel ) );
    balloon_orig = balloon_new; % = false(size(not_channel)); % see beginning of for loop
    balloon_orig( area_props( neck_cutoffs(i) ).PixelIdxList) = 1 ;
    balloon_edit = balloon_orig;
    balloon_new = balloon_orig;

    % combine algorithms shrink from rm_islands4, swell, fat_channel
    % the removed pixels in the fat channel should be eroded
    balloon_edit = bwdist(~balloon_edit);    % distance of balloon from channel
    balloon_edit( balloon_edit < channel_w_in_pix ) = 0;    % shrunken balloon
    balloon_edit = bwdist( balloon_edit );   % distance from shrunken balloon
%     balloon_new( balloon_edit > ( channel_w_in_pix * 2 ) & fat_new_channel ) = 0; % remove most of the string 
    balloon_new( balloon_edit > ( channel_w_in_pix + 2 ) & fat_new_channel ) = 0; % remove most of the string 
    balloon_new = bwareaopen( balloon_new, min_px_at_cutoff ); % remove remaining string blobs
    not_channel( balloon_orig & ~balloon_new ) = 0; %remove the string in from the channel

%     figure(31)
%     %     hold on
%     %     imagesc([ (old_fp & new_fp) ; ~not_channel + ~(old_fp & new_fp) ; not_channel ; balloons])
%         imagesc([ (old_fp & new_fp) , not_channel ; balloon_new , balloon_orig ])
%     %     hold off
%         axis equal
%         colorbar
%         pause(1)
% %         waitforbuttonpress

% from rm_islands7.m manual check of a region
% if i == region_index %this is the region we are looking for
%         one_region( area_ratio(i).PixelIdxList ) = true;
%         one_centroid = centroids( i-1,: )
%         idx = i;
% %         user_type = 
% end

end
% %% User classification section
% one_region( area_props( region_index ).PixelIdxList ) = true;
% % figure(9)
% % imagesc(one_region)
% one_reg_stat = regionprops(one_region);
% 
% % from rm_islands7.m manual check of a region
% figure(11)
%     set(gcf,'color','w')% 'position',[250 150 1000 800]) 
%     xmn = floor( max( one_reg_stat.Centroid(1) - 300, 1 ) ); % floor to make int
%     ymn = floor( max( one_reg_stat.Centroid(2) - 250, 1 ) );
%     xmx = floor( min( xmn + 600, size( not_channel, 2 ) ) );
%     ymx = floor( min( ymn + 500, size( not_channel, 1 ) ) );
% %     xmn = 1;
% %     ymn = 1;
% %     xmx = size(not_channel, 2)
% %     ymx = size(not_channel, 1)
%     d_dist = del_distance(ymn:ymx,xmn:xmx) .* one_region(ymn:ymx,xmn:xmx) ;
%     sc = max(max( d_dist ) );
%     channels = ~( old_fp & new_fp & mid_fp );
%     imagesc( [ channels( ymn:ymx,xmn:xmx ), not_channel( ymn:ymx,xmn:xmx );...
%              one_region( ymn:ymx,xmn:xmx ), d_dist / sc ] )
% 
%          title(sprintf(' area is %1.0f', one_reg_stat.Area) )
%     prompt = sprintf('Please enter a value:\n1 = crescent\n2 = neck cutoff\n3 = sm island\n4 = big island\n5 = compound island\n6 = other\n');
% %     prompt = 'What is the original value? ';
% % figure(gcf)
% % pause(0.01)
% user_type = input(prompt)
% 
%% this section is new July 2016, to remove small patches after cutting strings 
% removing strings can leave small patches, remove them now
% this was not included in rm_islands8.m but maybe should be
% not_channel = bwareaopen( not_channel, min_px_at_cutoff );% redundant with balloon_new = bwareaopen(... above

% the following version is more robust, but not needed bc of mid_fp method

% L = bwlabel( not_channel );                   % identify regions
% % remove biggest region from regionprops calculation to save time
% L(L==1)=0;
% area_props = regionprops(L,del_distance,'Area','MeanIntensity','PixelIdxList');
% big_crescents = find( [ area_props.Area         ] < min_px_at_cutoff * crescent_expansion & ...
%                       [ area_props.MeanIntensity] < crescent_dist_thresh );
% size(big_crescents)
% % remove big crescents
% not_channel( vertcat( area_props( [ big_crescents ] ).PixelIdxList ) ) = 0;


end