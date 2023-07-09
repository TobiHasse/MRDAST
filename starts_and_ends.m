function [ x_start, x_start_buffered, x_end, end_sim ] = ...
                            starts_and_ends( riv, riv2, pix_per_chan, B)
% Purpose:  Show the start and end nodes of the river planform throughout
%           the river simulation
%           Return suggested starting and ending pixels in the x (down
%           valley) direction
%           Return suggested end simulation time
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     June 2021

pixel_sz = 2 * B / pix_per_chan;  

% find x y limits of simulation area
limits = domain_limits(riv,riv2,pix_per_chan,B);
% initialize some variables
x_start_m  = zeros(1,numel(riv));
x_start_px = zeros(1,numel(riv));
x_end_m    = zeros(1,numel(riv));
x_end_px   = zeros(1,numel(riv));

% loop through each riv planform and extract the start and end coordinate
for i = 1:numel(riv)
    x_start_m(i) = riv(i).Xcl(1);
    x_end_m(i)   = riv(i).Xcl(end);
    % convert the coordinates in meters to pixel coordinates
    % bool_channel = node2grid( x, y, pixel_sz, limits);
    bool_channel = node2grid( [x_start_m(i) x_end_m(i)], [10 10],...
        pixel_sz, limits);
    x_start_px(i) = find(sum(bool_channel)>0,1,'first');
    x_end_px(i)   = find(sum(bool_channel)>0,1,'last');
end

% suggested limits
% suggested x_starting pixel is downstream of the farthest downstream
% starting node of the river channel
x_start = max(x_start_px); 
x_start_buffered = x_start + 20*pix_per_chan; % excluding odd initial bends

%% time_space_peak is the point where the product of x_end and simulation
% length are maximized
time_space_peak = ( 1:length(x_end_px) ) .* cummin( x_end_px );
[v, idx] = max(time_space_peak);    % returns the value and index
x_end = v/idx;                      % division because v is the product 
end_sim = idx;                      % the index is the suggested end time

%% Make figures to show results
hf = figure(10);
clf
set(hf,'color','w','position',[100 100 800 500])
px = 1:numel(riv);
[ax, h1, h2 ] = plotyy(px,x_start_m,px,x_start_px);
set(ax(1),'ylim',[0 2000])
set(ax,{'ycolor'},{'k'})
ylabel(ax(1),'X coordinate (meters) of starting node')
ylabel(ax(2),'X coordinate (pixels) of starting node')
title('Starting coordinate drifts down river valley')
legend('meters','pixels','location','northwest')


hf = figure(11);
clf
set(hf,'color','w','position',[100 100 800 500])
px = 1:numel(riv);
[ax, h1, h2 ] = plotyy(px,x_end_m,px,x_end_px);
yl = get(ax(1),'ylim');
set(ax(1),'ylim',[yl(1) yl(2)*1.5])
set(ax,{'ycolor'},{'k'})
ylabel(ax(1),'X coordinate (meters) of last river node')
ylabel(ax(2),'X coordinate (pixels) of last river node')
title('End of river is not always at the end of the valley')
text(6845,yl(2)*1.07,'End storage time simulation at 6845','Rotation',90);
legend('meters','pixels','location','northwest')

end