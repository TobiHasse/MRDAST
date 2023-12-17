%% Main CHAPTER ************ 2 ********************
% Purpose:  This script will call other functions to create the meandering
%           river simulation and place sediment deposits and track the
%           ages of the sediment eroded from those deposits.  It will step
%           through all the simulation and analysis of Chapter 2 of Tobias
%           Hasse's Dissertation.
% Attribution: The meandering river simulation code has been edited from
%           the original form which was published by Jon Shwenk:
%           http://onlinelibrary.wiley.com/doi/10.1002/2014JF003252/full
%           along with other functions as supplemetary info for the 
%           article: Schwenk, J., Lanzoni, S., & Foufoula?Georgiou, E. 
%           (2015). The life of a meander bend: Connecting shape and 
%           dynamics via analysis of a numerical model. 
%           Journal of Geophysical Research: Earth Surface.
%
% File repository: Some of the steps in the simulation take a lot of
%           computational time, particularly due to high RAM demands.  To
%           shortcut this, the user is encouraged to download select output
%           files from the file repository:
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     2015 - October 2021

%% generate parameter input files
% Rather than putting parameters directly in *.mat files, I have created
% commented code which will generate input parameter files.  
% This allows you to see commentary and variable definitions, or
% make adjustements for your simulations, But...
% if you play with parameters or node spacings, it is likely that you will
% encounter numerical instability. Instability will often cause problems
% in the 'enforce_spacing' subroutine, and Schwenk has added a warning and
% 'keyboard' command for debugging if that happens.

close all
clear
% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
% cd 'C:\Users\thasse\Documents\MATLAB\meander code\Ch2\bug fix' %
cd 'C:\Users\thasse\Documents\MATLAB\test' %

dt_save_years = 2;                  % save river planform every ## years

save_params_meander()               % parameters for meander model
% save_params_meander_Schwenk()     % as used by Schwenk (with commentary)

save_params_storage(dt_save_years)  % for computing storage times
save_params_deposition()            % sediment deposition 
save_initial_planform()             % initial random planform 
                                    % included for (reproducibility)

%% initialize variables for meander model

do_you_have_stats_toolbox = 0;      % enter 1 if you have the Statistics 
                                    % Toolbox, else enter 0
sim_time_ky = 15; %211;             % simulation time in thousands of years
                                    % note orig sim deleted after 205.55
sim_time = sim_time_ky * 1000;
dt_save = 12;                       % save every ## iterations
dt =     dt_save_years / dt_save;

%% Run the meander migration model


% Edit this code block ***********************************************

% outfile = sprintf('Ackerman_3chan_%ska',num2str(sim_time_ky))
% outfile = sprintf('Hasse_3chan_120-130 %ska',num2str(sim_time_ky))
% outfile = 'double check 2016 params'
% [river, B, mxub, freq_mig] = check2016run( do_you_have_stats_toolbox,...
%     outfile,sim_time,dt,dt_save); %

% i=1; CFO = 0.0036; Ain = 10;
% outfile = 'Schwenk low slope'
    
% i=1; CFO = 0.024; Ain = 3;
% outfile = 'Hasse 2016 low slope'
% [river, nodecount, B, mxub, freq_mig] = migration_model_TRH(CFO(i),...
%     Ain(i), do_you_have_stats_toolbox,outfile,sim_time,dt,dt_save); %

% i=1; CFO = 0.024; Ain = 3;     % Hasse
% i=1; CFO = 0.01; Ain = 16;     % Testing
% i=1; CFO = 0.0036; Ain = 10; % Schwenk
i=1; CFO = 0.005; Ain = 12; % new → lambda = 11.5 ± 7.9

outfile = '15ka_2023' %Ch2 2016 params'
[river, nodecount, B, freq_mig] = migration_model_TRH_Ch2(CFO(i), ...
    Ain(i), do_you_have_stats_toolbox, outfile, sim_time, dt, dt_save); 

% outfile = 'double check 2016 params'
% [river, B, mxub, freq_mig] = check2016run( do_you_have_stats_toolbox,...
%     outfile,sim_time,dt,dt_save); %

%% Save output
% original save code:
% save(sprintf('%s_2yr_A3_Cfo24_2Eo.mat',outfile),'river','B','-v7.3') 
% more flexible save code:
save(sprintf('%s_2yr_A%d_Cfo%0.3f_2Eo.mat',...
    outfile,Ain,CFO),'river','B','-v7.3') 

%% Schwenk's visualization of output (optional)
% this script plays animations of the entire simulation and each cutoff
% bend
visualize_Schwenk


%% save output smaller (dissertation settings)
% take 2 year river, convert to 30 year 'riv' and 'riv2' (offset)
% dt = 15 for Hasse dissertation.  With 'river' saved every 2 years, dt = 
% 15 yeilds a 30 year time step for deposition
% load from file in case user changes dt for storage time computation
load params_storage dt 

start_time = 100 + dt ;
end_time   = numel(river); %50001; %25001
riv = river ([start_time:dt:end_time]);
% original save code
% save(sprintf('%s_30yr_A3_Cfo24_2Eo.mat',outfile),'riv','B',...
%     'start_time','-v7.3') 
% more flexible save code
save(sprintf('%s_30yr_A%d_Cfo%0.3f_2Eo.mat',...
    outfile,Ain,CFO),'riv','B','start_time','-v7.3') 

% save output with 
start_time = 100 + dt+round(dt/2); 
end_time   = numel(river); %50001; %25001
riv2 = river ([start_time:dt:end_time]);
% original save code
% save(sprintf('Ackerman_211ka_30yr_A3_Cfo24_2Eo_offset.mat'),...
%     'riv2','B','start_time','-v7.3') 
% more flexible save code
save(sprintf('%s_30yr_A%d_Cfo%0.3f_2Eo_offset.mat',...
    outfile,Ain,CFO),'riv2','B','start_time','-v7.3') 



%% The following code is for analysing
% Meandering River Dynamics and Storage Time github repo: MRDAST
% https://github.com/TobiHasse/MRDAST
% MRDAST depends on the LOMR repo to generate the meandering river
% planform evolution and history





%% Read in the saved channel planforms and view

% ********** LOAD riv, B, riv2, TO RUN CODE BLOCKS BELOW ************

clear
cd 'C:\Users\thasse\Documents\MATLAB\test' %
% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
% dissertation files
% Hasse_211* files available at: https://doi.org/10.5281/ZENODO.5651841.
load('Hasse_211ka_30yr_A3_Cfo24_2Eo_offset.mat')
load('Hasse_211ka_30yr_A3_Cfo24_2Eo.mat')

% other files for checking
% load('code cleanup_30yr_A16_Cfo0.010_2Eo.mat')
% load('code cleanup_30yr_A16_Cfo0.010_2Eo_offset.mat')
% load('code cleanup_30yr_A10_Cfo0.004_2Eo.mat')
% load('code cleanup_30yr_A10_Cfo0.004_2Eo_offset.mat')

% load('5ka_2023_30yr_A3_Cfo0.024_2Eo.mat')
% load('5ka_2023_30yr_A3_Cfo0.024_2Eo_offset.mat')
% load('50ka_2023_30yr_A12_Cfo0.005_2Eo.mat')
% load('50ka_2023_30yr_A12_Cfo0.005_2Eo_offset.mat')
% load('15ka_2023_30yr_A12_Cfo0.005_2Eo.mat')
% load('15ka_2023_30yr_A12_Cfo0.005_2Eo_offset.mat')

beep; pause(1); beep



%% start and ending nodes of simulation

% estimated run time ::::::::::::::::::: 2 minutes

load params_storage.mat
[ x_start_rx, x_start_buffered, x_end_rx, end_sim_rx ] = ...
                            starts_and_ends( riv, riv2, pix_per_chan, B);

fprintf(['Farthest downstream starting pixel is %d but the suggested\n',...
    ' upstream starting pixel is %d to buffer some odd meander bends',...
    ' near the starting node. \nStarting pixel for dissertation was',...
    ' %d \n \n'], x_start_rx, x_start_buffered ,x_start)

fprintf(['Suggested end simulation model step is %d. \nEnd simulation ',...
    ' for dissertation was %d\nWithin that time the farthest upstream',...
    ' pixel for the end of the channel was %d. \nLast pixel for',...
    ' dissertation was %d \n \n'], end_sim_rx, end_sim, x_end_rx, x_end )

%% Create figures and files showing an overview of the model 

% estimated run time ::::::::::::::::::: 18 hours

% This will create figures showing the topography, age, selected 
% stratigraphy, floodplain width, migration rate, etc. throughout the
% simulation.  It will also save some files including the figures, data
% files, and a data file map of the floodplain area which was visited by
% the channel at least once.

display_model(riv, riv2, B)
% display_model_topo_video(riv, riv2, B)

%% animate channel occupation 

% estimated run time ::::::::::::::::::: 1 minute

% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other

outfile = 'junk.gif'
slice_interval = 6;


[valley_w_pxl,occ_bool,~] = chan_occ(outfile,slice_interval);

% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other

%% slicing the model to manage memory demands

% estimated run time ::::::::::::::::::: 10 seconds

% Dependancy:   Must run after display_model.m has successfully completed
%               preferably there is a floodplain data file titled: 
%               "Previously Occupied by Channel. Time step 6845.mat"
%               you may need to edit the slices.m function or file name

% for a model domain containing about 2,000 pixels across the valley and
% 3,000 pixels down valley and 7,000 time steps the array containing that
% is 42 billion pixels, at double precision this is 40 to 80 GB of RAM !!!
% (per array). To manage these memory demands the model is sliced into 
% short reaches and sediment is built and the storage times measured.

% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other

load params_storage.mat

% end_sim = end_sim_rx; % to use the recommended simulation time
show_slice_figures = true;

% show the slices
[x_starts, x_ends, ymn, ymx] = slices( array_size, end_sim, ...
    x_start, x_end, show_slice_figures);
%% Zig Zags Ch 2 Figure 14
% roll through the slices in the model and make diagrams similar to
% sequence stratigraphy Wheeler Diagrams, save thos figures to file

% estimated run time ::::::::::::::::::: 16 minutes

% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other

% zigzags will require you to select the folder with the outputs from
% display_model.m 

zigzags


%% Create deposits and calculate storage times

% estimated run time ::::::::::::::::::: 40 to 100     ****** DAYS *******

% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other

restart = false
slice = 0; % zero for first slice

show_figures = false;           % figures require significant RAM and time

deposit_storage_time(riv, riv2, B, show_figures,  restart,...
     [], [], [], [], slice, [],[],[],[])


%% Restart deposit_storage_time after some slices have been calculated

% estimated run time ::::::::::::::::::: 40 to 100     ****** DAYS *******
% run time is the remainder of the run started above

% this functionality is useful if there is a power outage or crash or if 
% the program must be stopped and restarted.  This is set up to restart 
% after the most recently completed floodplain slice and the input file
% name must be updated. This program will ignore an incomplete slice.


clear
restart = true
show_figures = false;           % figures require significant RAM and time
cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
% load planform files
load('Hasse_211ka_30yr_A3_Cfo24_2Eo_offset.mat')
load('Hasse_211ka_30yr_A3_Cfo24_2Eo.mat')
% load storage time files for the MOST RECENT slice
load distr_3chan_205k_30yr_3.mat    % "_##.mat"  UPDATE THE SLICE NUMBER ##

deposit_storage_time(riv, riv2, B, show_figures, restart,...
    x_starts, x_ends, ymn, ymx, slice, ...
    pt_bar_dist, vert_dist, vd_age_dist, pb_age_dist )


%% age distribution of area nearby to the channel

% estimated run time ::::::::::::::::::: 12-24 hours

lam_vc = 10.5 % characteristic meander wavelength (down valley axis)
mb_age_lim = 333; %333 for 9990 years ? 10,000 years
mb_start_step = 4000; % 4000 to start measuring the meander belt at 120kyr

meander_belt_age( lam_vc, mb_age_lim, mb_start_step, riv, riv2, B)


%% Partition storage and age distributions upstream, downstream & combined

% run the block of code above for compute_dists prior to this block of code
% so that the outputs match.

% estimated run time ::::::::::::::::::: 40 minutes

% after deposit_storage_time has completed the user can select the farthest
% downstream slice of both the upstream and downstream reach.
% extract_triangles will create and save 3 files including the point bar
% and floodplain storage time and age distributions for the upstream reach,
% the downstream reach, and the combination of both reaches.  
% extract_triangles also cuts off the final few thousand years of
% simulation because the final node of the simulated river drifted far
% upstream
clear
% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other

file_name = 'distr_3chan_211k_30yr_'; % input file
% last slice of each reach
upstream_slice = 16;
downstream_slice = 36;

load('params_storage.mat','x_start','x_end')
save('params_reaches.mat','upstream_slice','downstream_slice',...
    'x_start','x_end','-v7.3')

extract_triangles( file_name, upstream_slice, downstream_slice )
%% Extract distributions for chapter 1 figures and fits CDF & PDF

% estimated run time ::::::::::::::::::: 30 seconds

% after deposit_storage_time has completed and generated output files the 
% user can select the farthest downstream slice of both the upstream and 
% downstream reach and create storage time distributions accordingly
% compute_dists_triangles also cuts off the final few thousand years of
% simulation because the final node of the simulated river drifted far
% upstream
% compute dists will save the distributions in files titled like this:
% 'TRH 205k ALL storage time cdf pdf.mat' and also for up and downstream

% requires params_reaches.mat and output from extract_triangles

compute_dists_triangle( file_name, upstream_slice, downstream_slice )

%% Fit CCDF surviovor function distributions with various models to
% generate figures 8, 9 & 10

% estimated run time ::::::::::::::::::: 30 seconds

distribution_fits
%% For each of the survivor functions, adjust one fit parameter and track
% R squared.  Then find the ranges within a small tolerance of the maximum
% R squared

% estimated run time ::::::::::::::::::: 75 minutes

[RS, ranges5, coef] = distr_r_squared_ranges()

descr=strcat('These variables were created in distr_r_squared_ranges.m',...
    ' and are used to make figures 3, 11, & 12');
save('R_squared_ranges.mat','RS', 'ranges5', 'coef','descr','-v7.3')

%% Create figures showing model fits and ranges for the survivor functions
% using both the matlab fit and the R squared method 
% figure 3 shows the sensitivity of R squared as the scale break parameter
% (T in the dissertation, xo in the code) is forced to change over a wide
% range of values

% estimated run time ::::::::::::::::::: 10 seconds

% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other

load('R_squared_ranges.mat')
% Make a figure showing how R^2 varies for the transition paramter xo (T)
Ch2fig3_R2_ranges(RS)


% Make figures showing the ranges and the maximum R^2
Ch2fig11_fig12_ranges(ranges5, coef)


%% Make figure showing probability density functions

% estimated run time ::::::::::::::::::: 10 seconds

Ch2fig13_pdf

% analysis of storage time continues in Main_Ch4.m

