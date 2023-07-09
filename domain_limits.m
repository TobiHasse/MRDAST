function limits = domain_limits(riv,riv2,pix_per_chan,B)
% Purpose:  find the absolute minimum and maximum X and Y coordinates of
%           river centerline for converting centerlines to rasterized 
%           floodplain map
% Author:   Tobias Hasse    tobiack@udel.edu
% Date:     June 2021

% find x y limits of simulation area
% regular, one riv variable limits
Xmin = min(vertcat(riv(1:end).Xcl)); Ymin = min(vertcat(riv(1:end).Ycl));
Xmax = max(vertcat(riv(1:end).Xcl)); Ymax = max(vertcat(riv(1:end).Ycl));
limits1 = [ Xmin Xmax Ymin Ymax ];
% two riv variable limits
Xmin = min([ Xmin;vertcat(riv2(1:end).Xcl) ]); 
Ymin = min([ Ymin;vertcat(riv2(1:end).Ycl) ]);
Xmax = max([ Xmax;vertcat(riv2(1:end).Xcl) ]); 
Ymax = max([ Ymax;vertcat(riv2(1:end).Ycl) ]);
limits = [ Xmin Xmax Ymin Ymax ]; 
% additional limit allignment
% if new x_min is to the left, shift an even number of pixels to the right
shift_px = ceil(max(limits1(1)-limits(1),0)/14); 
% shift the xmin an even number of pixels to the left 
% (4 pixels in a 5 pixel channel)
limits(1) = min( limits(1), limits1(1) - shift_px/pix_per_chan * 2 * B ); 

end