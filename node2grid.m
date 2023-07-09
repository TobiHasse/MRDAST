function [b_A] = node2grid( x, y, pxsz, lim )
% This function makes a binary raster out of x,y point data 

% Written by Tobias Hasse tobiack@udel.edu May, 2015

% INPUTS
%     x,y       - the x and y point coordinates, vectors of matching length
%     lim       - [ Xmin Xmax Ymin Ymax ]
%     pxsz      - the pixel size dimension
% OUTPUTS
%     b_A       - a boolean array true where (x,y) falls inside the pixel

%% make some test data
% x = [1:1:1000];
% y = 100 * sin(x/100);
% lim(3)= min(y);
% y = y - lim(3);  % translates y into positive territory
% 
%% translate coordiantes into the positive domain
x = x - lim(1);
y = y - lim(3);

% create 2d indecies for the pixels that contain nodes
i = floor (( pxsz + x )/pxsz );
j = floor (( pxsz + y )/pxsz );


%% find the number of rows & cols & create a blank boolean array
% nrows = max(j);  
% ncols = max(i);
nrows = ceil( ( lim(4)-lim(3) )/pxsz );
ncols = ceil( ( lim(2)-lim(1) )/pxsz );
b_A = false(nrows,ncols);
% b_A = zeros(nrows,ncols);

%% linearize the indecies and create the bool array
id =  j + nrows*(i-1);
b_A( id ) = 1;


end