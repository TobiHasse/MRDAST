function [mb] = Camporeale_mb(chan_cl,lambda_px)
% this function is written by Tobi Hasse January 16 2020 to take a channel
% centerline and return a meander belt swath that is 3.4 x the down valley
% meander wavelength as per Camporeale 2005
% inputs      chan_cl      boolean raster of the channel
%             lambda_px    meander wavelength in pixels (down valley aka cartesian)
% outputs     mb           boolean raster of the meander belt size(chan_cl)

mb = false(size(chan_cl));  % place holder for output
bwd = bwdist(chan_cl);
% figure(2)
%     imagesc(bwd)

sausage = false(size(mb));
sausage(bwd< 2.5*lambda_px )=true;  % 200 for lambda_px = 80
bb = regionprops(sausage,'boundingbox');
bb = bb.BoundingBox;
bb(1:2) = ceil(bb(1:2));   % order is c1 r1 r_width c_height
bb(3:4) = bb(1:2)+bb(3:4)-1;  % -1 required or there is a single zero pixel at the right & bottom edges
% clip the image to the bounding box, and pad the top and bottom with 0
sausage = padarray( sausage(bb(2):bb(4),bb(1):bb(3)) ,[1,0], 'both');  
% Note: padarray probably renders some other code below redundant (Jan 26 2020)
% Also: will give a slightly different result when the meander belt
% 'sausage' touches the top or bottom edges...

% figure(3)
%     imagesc(sausage)

skel = ~watershed(sausage); % only works if sausage touches at least 2 image borders
% clearvars sausage
% figure(5)
%     imagesc(skel)    
%     title('skelaton')
% 
% no_spur = bwmorph(skel,'spur',1000);
% figure(6)
% imagesc(no_spur)
% title('no spurs')
%
% p = true(size(skel));             % make perimeter boolean 
% p(2:end-1,2:end-1) = false;    % make only the perimeter true
% p_skel = unique(skel(p));
% ht = size(skel,1);
% ends = skel;                       % copy 
% ends(:,2:end-1) = false;    % make only the left & right end pixels true
% [row col] = find(ends);     % find the perimeter pixels of the skeleton
% cols = unique(col);
% lhs = row(col==cols(1))-ht/2;
% max(abs(lhs));
% rhs = row(col==cols(2))-ht/2;
% max(abs(rhs))
% Take watershed and find path along skeleton from central contact with 
% left and right edges of binary image
% find coordinates for end nodes
tmp=skel(:,1);
ids = find(tmp);
idctr = abs(ids - length(tmp)/2);
r1=ids(idctr==min(idctr));
c1 = 1;
tmp=skel(:,end);
ids = find(tmp);
idctr = abs(ids - length(tmp)/2);
r2=ids(idctr==min(idctr));
c2 = size(skel,2);
% skel = skel;
% if the skel does not touch either the left or right side
if isempty(r1*c1*r2*c2)
    
    bbs = regionprops(skel,'boundingbox');
    bbs = bbs.BoundingBox;
    bbs(1:2) = ceil(bbs(1:2));   % order is c1 r1 r_width c_height
    bbs(3:4) = bbs(1:2)+bbs(3:4)-1;  % -1 required or there is a single zero pixel at the right & bottom edges
    skel_sm = skel(bbs(2):bbs(4),bbs(1):bbs(3));  % clips the image to the bounding box

tmp=skel_sm(:,1);
ids = find(tmp);
idctr = abs(ids - length(tmp)/2);
r1=ids(idctr==min(idctr));
c1 = 1;
tmp=skel_sm(:,end);
ids = find(tmp);
idctr = abs(ids - length(tmp)/2);
r2=ids(idctr==min(idctr));
c2 = size(skel_sm,2);
% before_adjustment = [c1 r1 c2 r2]
c1= c1 + bbs(1) -1;
c2= c2 + bbs(1) -1;
r1= r1 + bbs(2) -1;
r2= r2 + bbs(2) -1;


end


% https://blogs.mathworks.com/steve/2011/12/13/exploring-shortest-paths-part-5/
D1 = bwdistgeodesic(skel, c1, r1, 'quasi-euclidean');
D2 = bwdistgeodesic(skel, c2, r2, 'quasi-euclidean');

D = D1 + D2;
D = round(D * 8) / 8;

D(isnan(D)) = inf;
skeleton_path = imregionalmin(floor(D/10));
% P = imoverlay(skel, imdilate(skeleton_path, ones(3,3)), [1 0 0]);
% P = skel + skeleton_path;
% figure(6)
%     imagesc(P)%, 'InitialMagnification', 200)
%     hold on
%     plot(c1, r1, 'g*', 'MarkerSize', 15)
%     plot(c2, r2, 'g*', 'MarkerSize', 15)
%     hold off
%  CAMPOREALE 2005 MEANDER BELT
skel_dist = bwdist(skeleton_path);
% figure(7)
%     imagesc(skel_dist)
skel_bin=false(size(skel_dist));
%
% down valley wavelength = 16 w (32B) for the expected model input
% Camporeale 2005 suggests a meander belt width of 3.4 lambda =
% 3.4*16w=3.4*16*5pixels=272 pixels for Hasse model
% mb_w_pix = 272;
mb_w_pix = 3.4 * lambda_px;
mb_hw_pix = mb_w_pix/2; % since using anodi to delineate only half width of meander belt needed
skel_dist(skel_dist<mb_hw_pix)=1;
skel_dist(skel_dist>  1)=0;

mb(bb(2):bb(4),bb(1):bb(3)) = skel_dist(2:end-1,:);
end
