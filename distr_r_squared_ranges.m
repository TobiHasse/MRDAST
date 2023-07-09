function [RS, ranges5, coef] = distr_r_squared_ranges()
% Purpose:  This function will find the ranges for each of the survivor
%           function (CCDF) parameters figures 3, 11 & 12.  It is reordered
%           and cleaned up code from the file curv_fits_array.m
%           originally written in December 2018
%
%           Dependancies:  input files created by compute_dists.m
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     edited October 2021, written December 2018
clear

% Load files and add CDF variables to an array of CCDF's aka survior
% functions or Complementary Cumulative Distribution Functions
load('TRH 205k upstream storage time cdf pdf.mat')
ccd(2,:)=1-CDFp;
ccd(3,:)=1-CDFv;
load('TRH 205k downstream storage time cdf pdf.mat')
ccd(4,:)=1-CDFp;
ccd(5,:)=1-CDFv;
% load this one last to use other variables in file for plots
load('TRH 205k ALL storage time cdf pdf.mat')
ccd(1,:)=1-CDFall;

fx=px(2:end);   % grab the x coordinates for plotting

% loop through each coefficient, incrementing one coefficient slowly
% accross a range and storing the resulting R^2 value of the fit
tic
for i = 1:5
    [RS{i},coef] = make_fits_r_squared(fx,ccd,i);
    i
    datestr(clock)
end
toc
% DEBUG
% load('temp RS.mat')

% obtain fit objects for all the survivor functions
[fits,coef] = survivor_fit_eqn_8(fx,ccd,1);


% find the minimum and maximum of the range on the abscissa where R2 is
% within my_tolerance of the maximum R^2
my_tolerance = 0.001;
range_xo = find_range(RS{1},my_tolerance);
range_b  = find_range(RS{2},my_tolerance);
range_c  = find_range(RS{3},my_tolerance);
range_m  = find_range(RS{4},my_tolerance);
range_n  = find_range(RS{5},my_tolerance);

% store the ranges in an array
ranges5=[];
ranges5(:,:,1)=range_b;
ranges5(:,:,2)=range_c;
ranges5(:,:,3)=range_m;
ranges5(:,:,4)=range_n;
ranges5(:,:,5)=range_xo;  % NOTE ranges are reordered *****!!

end % end function distr_r_squared_ranges


function range_out = find_range(rs, my_tolerance)
% find the minimum and maximum of the range on the abscissa where R2 is
% within my_tolerance of the maximum R^2


search_value=max(rs(1:5,:)') - my_tolerance;

for i = 1:5%size(rs,1)-1
    % find 
    [bool,idm(i)]=find(rs(i,:)==max(rs(i,:)));
    [bool,idf]=find(rs(i,:)>search_value(i),1,'first');
    [bool,idl]=find(fliplr(rs(i,:))>search_value(i),1,'first');

    idl = length(rs)-idl+1;
    x = [rs(i,idf),rs(i,idf-1)];
    v = [rs(6,idf),rs(6,idf-1)];
    vf(i)=interp1(x,v,search_value(i));
    try
        x = [rs(i,idl),rs(i,idl+1)];
        v = [rs(6,idl),rs(6,idl+1)];
        vl(i)=interp1(x,v,search_value(i));
    catch
        vl(i)=rs(6,idl);
        disp('top of range at limit')
    end
end
vm=rs(6,idm);
range_out=[vl',vm',vf'];

end % find_range


function [rs,coef] = make_fits_r_squared(fx,ccd,my_case)
% fit each survivor function in ccd using equation 8 while adjusting one
% parameter and tracking the R^2 of the fit

jj = 0;

switch my_case
    case 1              % increment xo
        jmn = 200; jstp = 1; jmx = 2000;
%         jmn = 200; jstp = 800; jmx = 2000; %fast version
    case 2              % increment b
        jmn = 50000; jstp = 5000; jmx=240000;
%         jmn = 50000; jstp = 100000; jmx=240000;
    case 3              % increment c
        jmn = 0.4; jstp = 0.03; jmx = 1.2;
%         jmn = 0.4; jstp = 0.4; jmx = 1.2;
    case 4              % increment m
        jmn = 0.8; jstp = 0.03; jmx = 1.1;
%         jmn = 0.8; jstp = 0.1; jmx = 1.1;
    case 5              % increment n
        jmn = -0.002; jstp = 0.0001; jmx = -0.0001;
%         jmn = -0.002; jstp = 0.0009; jmx = -0.0001;
    otherwise
        disp(strcat('Please choose a switch between 1 and 5.',...
            ' ~dist_r_squared_ranges.m'))
        rs = nan;
        return
        
end
%  for j=200:1:2000 %increment xo
% for j=50000:5000:240000 %increment b
% for j=.4:.03:1.2 %increment c
% for j=.8:.03:1.1 %increment m
% for j=-.002:.0001:-.0001 %increment n
for j = jmn:jstp:jmx
    jj=jj+1;
    
    for i = 1:size(ccd,1)
        % equation 8
        ft = fittype(strcat(' (m*exp(n*x))./( 1+exp( 1 *(x-xo) ) ) + ',...
            'c*((b^c*m*exp(n*xo)/(c*(xo^-c-b^-c)))/(b^c+m*exp(n*xo) / ',...
            '(c*(xo^-c-b^-c))))^(1/c)^c*(x^-c-b^-c)./',...
            '( (1-(((b^c*m*exp(n*xo)/(c*(xo^-c-b^-c)))/(b^c+m*exp(n*xo)/',...
            '(c*(xo^-c-b^-c))))^(1/c)/b)^c) *( 1+exp( 1 *(xo-x) ) ) ) '));
        % variable titles, xo below is scale break T in dissertation
        %        b     c       m     n      xo
        sp = [156000 .845      1 -0.00138 700];  % fit bounds
        up = [205000  2        2    0    5000];
        lo = [ 50000 .2      .05    -1    200];
        switch my_case
            case 1              % bounds for forcing xo
                %        b     c       m     n      xo
                sp = [156000 .845      1 -0.00138 j];  % fit bounds
                up = [205000  2        2    0    j+1];
                lo = [ 50000 .2      .05    -1    j-1];
                
            case 2              % bounds for forcing b
                %        b     c       m     n      xo
                sp = [j      .845      1 -0.00138 700];  % fit bounds
                up = [j*1.001 2        2    0    5000];
                lo = [j*.999  .2      .05    -1    200];
            case 3              % bounds for forcing c
                %        b     c       m     n      xo
                sp = [156000   j          1 -0.00138 700];  % fit bounds
                up = [205000  j*1.001      2    0    5000];
                lo = [ 50000  j*.999      .05    -1    200];
            case 4              % bounds for forcing m
                %        b     c       m     n      xo
                sp = [156000 .845      j    -0.00138 700];  % fit bounds
                up = [205000  2      j*1.001    0    5000];
                lo = [ 50000 .2      j*.999    -1    200];
            case 5              % bounds for forcing n
                %        b     c       m     n      xo
                sp = [156000 .845      1     j     700];  % fit bounds
                up = [205000  2        2  j*.999    5000];
                lo = [ 50000 .2      .05  j*1.001   200];
            otherwise
        end
        
        % do the fitting
        [fobj,gof(i),out(i)] = fit(fx',ccd(i,:)',ft,'Startpoint',sp,...
            'lower',lo,'upper',up,'tolfun',1e-10);%,'maxfunevals',2000);
        %     fobj % optional, print the fit to screen
        rs(i,jj)=gof(i).rsquare;    % store r squared
        rs(size(ccd,1)+1,jj)=j;     % j is empty
        %     fits{i}=fobj;               % store the fit object in a cell array
        % use the struct coef to store useful info
        coef.names = coeffnames(fobj);
        coef.vals(i,:)=coeffvalues(fobj);
        temp = confint(fobj);       % confidence interval
        coef.intervalLo(i,:) = temp(1,:);
        coef.intervalHi(i,:) = temp(2,:);
        
        % show the results in a figure
        %     figure(99)
        %     plot(fobj,fx,ccd(i,:))
        %     ylim([10^-7 10^2])
        %     xlim([10 10^6])
        %     title(sprintf('survivor function %d of %d', i, size(ccd,1) ) )
        %     set(gca,'yscale','log',...
        %         'xscale','log')
        %     xlabel(sprintf('%s',out(i).message))
    end %i
    % j
end %j

end % function make_fits_r_squared

