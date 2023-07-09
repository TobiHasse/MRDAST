function meander_belt_age( lam_vc, mb_age_lim, mb_start_step, riv, riv2, B)
% Purpose:  This file is set up to analyze the river planform output from 
%           the Schwenk / Hasse meander model, determine the meander belt
%           using two methods, one by Camporeale 2005, and the other 
%           tangent to the widest meander bends.
%           The age distribution of point bar deposits within the meander 
%           belt is detrmined up to an age set by mb_age_lim.

%           Figures will be saved in the .png format and the user is 
%           encouraged save the final figure in .svg format for further 
%           editing in a vector graphics editor such as inkscape.
% Author:   Tobias Hasse tobiack@udel.edu 
% Date:     2016 - revised 2021
% former file name: test_howard434s3box3chanAgeOnlyMB
% *************Figure 15 chapter 2 *****************

fprintf('Welcome to: "%s" \n',...
        fullfile(fileparts(mfilename('fullpath')),mfilename))
time_now = datestr(clock)

load params_storage.mat
load params_deposition.mat
end_sim = min(end_sim,numel(riv))

% lam_vc = 10.5 % characteristic meander wavelength (down valley axis)
% mb_age_lim = 4 % meander belt maximum age limit to check, should be 333 
                % for 10,000 year meander belt (since dt = 30 this give you
                % 9,990 years ? 10,000
% Adapted January 2019 by Tobias Hasse to compute the age distribution of
% the floodplain in the meander belt,

program_start = clock

disp_progress = round( display_progress_years / time_step_years )  
% disp_progress = 33   % 1 DEBUG show me every time

display_stops = [ disp_progress : disp_progress : size( riv, 2 ) ,...
    size( riv, 2 ) ];

%% find x y limits of simulation area
limits = domain_limits(riv,riv2,pix_per_chan,B);
%     
algorithm_start = clock
% pick start time for loop
i=mb_start_step; % 4000 for starting to measure the meander belt at 120kyr
i_start = i-mb_age_lim-6 % start early so floodplain has age up to 
                         % mb_age_limit at the 'start' of measuring mb age
% i_start = 4266 - 1-333;%ard code restart to fix select meander belt maps
i_start = max(2, i_start);%USUALLY start at 2 can't start at negative index
% i_start = 2; % for code checking 2021

%% Make boolean array of channel location
% Initialize the array (needed before loop for prev & occupation
[ x, y ] = doInterpm_TRH( riv( i_start-1 ).Xcl, riv( i_start-1 ).Ycl,...
    pixel_sz/interp_nds_pix, 'lin' ); % interpolate nodes along centerline 
                                      % ensure pixel channel is continuous
bool_channel = node2grid( x, y, pixel_sz, limits); % boolean raster

% Make some more arrays
deposit_step   = zeros( size( bool_channel ) ); % time step at depostion
fp_old_and_new = false( size( bool_channel ) ); % all floosplain

% initialize variables for calculating meander belt age distribution
mbsolid    = zeros(end_sim,333);   % Camporeale 2005 age distribution
closesolid = mbsolid;  % tangent age distribution 'close' to channel extent
mba        = zeros(end_sim,1);     % area of Camporeale meander belt 
closeMBa   = mba;                  % area of tangent meander belt

% Distance transform on channel centerline, make channel pixels negative
old_fp_distance = double( bwdist( bool_channel ) - ...
    floor( pix_per_chan / 2) );  
old_fp_distance( old_fp_distance < 0 ) = 0; % make channel pixels zero

tic
%% ****************** why reload? is it to touch up the run?
% load('TRH 205ky meander belt mb solidity and maps Jan 2020 Complete.mat')
%     mb(6845,333) =0; %since variables were not preallocated, expand array
%     ar(6845,333) =0;
%     solid(6845,333) =0;
%     mbsolid(6845,333) =0; % already preallocated above
%     closesolid(6845,333) =0;
% for i = i_start :1: 67%  % USUALLY 2:6845 OR 2:numel(riv) 
for i = i_start :1: end_sim % OR 4000-333:6845 FOR meander belt age
    % rasterize channel position
    [ x, y ]    = doInterpm_TRH( riv( i ).Xcl, riv( i ).Ycl, ...
        pixel_sz/interp_nds_pix, 'lin' );
    bool_channel = node2grid( x, y, pixel_sz, limits);
    %Distance transform on channel centerline, make channel pixels negative
    new_fp_distance = double( bwdist( bool_channel ) - ...
        floor( pix_per_chan / 2) );  
    new_fp_distance( new_fp_distance < 0 ) = 0; % make channel pixels zero

    if any( display_stops == i ) 
       disp(sprintf('Current Time:  %s \n i = %d', datestr(clock) , i ))
    end

    % make Camporeale meander belt based on new channel planform
    % for ppt discussion, used 16 as wavelength (lambda) using the fft
    % method
%   mbmap = Camporeale_mb(bool_channel,16*pix_per_chan); % chan width 5 pix
    % updated wavelength (10.5) using excess velocity method
    mbmap = Camporeale_mb(bool_channel,round(lam_vc*pix_per_chan));
%% tangent meander belt map, initially I used a structuring element
%     se=strel('disk',200);   % ********* STRUCTURING ELEMENT SHENANNIGANS 
    % the MATLAB algorithm used to make the disk structuring element run 
    % faster creates a non equilateral octagonal structuring element close 
    % in shape to an exribed irregular octagon.  To get a true circle I 
    % tried using strel( 'arbitrary' , nhood ) (see below)
    % 400 is the diameter not the radius which the strel command uses
%     circle_nhood = false(400);  
%     circle_nhood( floor( ( size( circle_nhood ) + 1 )/2 ) ) = true;
%     circle_nhood( bwdist ( circle_nhood, 'euclidean' ) < ...
%         length( circle_nhood )/2 ) = true;
%     se = strel( 'arbitrary', circle_nhood);
%   % The 'arbitrary' se is VERY slow (117 seconds).  It can be reproduced 
    % with BW dist (1.05 seconds): make the channel belt 2x wider (say, 
    % 400 pixels) then make it 1x narrower (say 200 pixels) only central 
    % portions remain.
%   % Algorithm for circular imclose on a river channel (using bwdist)
    temp=bwdist(bool_channel,'euclidean');  % distance from the channel
    temp(temp<200)=0;                       % erase the desired area
    temp=bwdist(temp,'euclidean');  % find the distance from the margin
    temp(temp<200)=0;               % erase the margin
    temp(temp>0)=1;                 % make the image binary
    
    closeMB = temp;                 % store in my variable
%     closeMB = imclose(bool_channel,se); % structuring element method

    % get the channel HALF WAY between old and new channel position 
    [ x, y ]    = doInterpm_TRH( riv2( i-1 ).Xcl, riv2( i-1 ).Ycl, ...
        pixel_sz/interp_nds_pix, 'lin' );
    bool_channel = node2grid( x, y, pixel_sz, limits);
    %Distance transform on channel centerline, make channel pixels negative
    mid_fp_distance = double( bwdist( bool_channel ) - ...
        floor( pix_per_chan / 2) );  
    mid_fp_distance( mid_fp_distance < 0 ) = 0; % make channel pixels zero

    % remove previous channel and spurrious 'islands' created by fast
    % migration AKA 'balloons' vs 'crescents'
    fp_old_and_new = rm_islands9( new_fp_distance, mid_fp_distance,...
        old_fp_distance, px_size_thresh, pix_per_chan );
    % ~fp_old_and_new is the channel and point bar and any channel which is
    % assigned a deposition age here will be overwritten when it is
    % abandoned in the future
    deposit_step( ~fp_old_and_new ) = i; % deposition 'birthday'

    % the index range 345:2240 is the upstream and downstream reach
    mba(i)=sum(sum(mbmap(:,345:2240)));  %Camporeale 2005 meander belt area
    closeMBa(i)=sum(sum(closeMB(:,345:2240))); % tangent belt area
    
    current_age = i-deposit_step(:,345:2240); % age now at time = i
    tr=false(size(current_age)); % make a boolean for 'young' floodplain

if i> i_start + mb_age_lim %33 % 333 % 4000
    for age = 2:min( mb_age_lim , max( max ( deposit_step ) ) ) %333
        tr(current_age<age)= true; % region of floodplain younger than age
        tr = bwareaopen(tr,65000);  % remove bends poking into the reach 
                                    % from outside (up or down)stream

        mbsolid(i,age)=sum(sum( mbmap(:,345:2240) & tr))/mba(i);
        closesolid(i,age)=sum(sum( closeMB(:,345:2240) & tr))/closeMBa(i);
    end % for age
end % if i_start

% reset arrays for next iteration
    old_fp_distance = new_fp_distance;
%%    
    if any( display_stops == i ) 
        pct = i / display_stops(end)*100;
        frac_complete = pct / 100 
        elapsed_algorithm_time = datenum(clock) - datenum(algorithm_start);
        disp(sprintf('Current Time:  %s \nest completion %s',...
            datestr(clock) , datestr( datenum( algorithm_start ) +...
            elapsed_algorithm_time / frac_complete ) ) )

        try
        figure(4)
%             imagesc(mbmap+bool_mb+bool_channel)
            imagesc(mbmap(:,345:2240)+tr)
            colorbar
            axis equal
            title('meander belt map')
            drawnow
        figure(44)
            imagesc(closeMB(:,345:2240)+tr)
            colorbar
            axis equal
            title('close meander belt map')
            drawnow
        end
    if i > mb_start_step % errors 
        %         try
            Ch2fig15_meander_belt(closesolid,mbsolid,mb_start_step,i )
%         catch
%             fprintf('Skipping figure inside: "%s" \n',...
%                 fullfile(fileparts(mfilename('fullpath')),mfilename))
%         end
    end % if
        figure(11)
            px = [1:length(mba)];
%             plot(mba)
            plot(px,mba,px,closeMBa)
            xlim([floor(i_start/10)*10 ceil(i/10)*10])
            legend('Camporeale meander belt area',...
                'tangent menader belt','location','south')
            title(sprintf(' meander belt area. iteration %d',i))
            drawnow
     end %dispaly stops
end %i
algorithm_end =clock
agorithm_duration = algorithm_end - algorithm_start
%%

descr = ['mb ar and solid are the area/perimeter, area, and solididty of the convex shell of the floodplain deposits less than a given age at a time step'...
' mbsolid & closesolid are the fractional age distribution of the nearby floodplain'...
' mbsolid is the solidity of the re-explored portion of the floodplain which is within the mender belt as defined by Comporeale 2005 as within 3.4 meander wavelengths from the valley axis.'...
' closesolid is computed using bwdist to mimic image processing closing an image with a circular structuring element of 200 pixel radius'... 
' (MATLAB uses a disk which is octagonal not circular)'... 
' mbmap, closeMB & bool_channel are the final Camporeale meander belt, final close meander belt and channel.'...
' RX a plot of solid" (Transposed) with only markers but also see plotSolid.m'];


save('TRH 205ky meander belt mb solidity and maps Jan 2020 Complete 2.mat','i','closesolid','mbsolid','mbmap','closeMB','bool_channel','descr','-v7.3')

beep on % tell me you're finished
beep;pause(1);beep;
disp('You have reached the end of meander_belt_age.m')
disp('You have keyboard control of the variables.')
disp('You are encouraged to save variables and figures as desired')
disp('.svg and .fig are reccommended for editing.')
disp('Type dbcont or dbquit to end keyboard control')
keyboard % give user control dbquit or dbcontinue


end % function 

