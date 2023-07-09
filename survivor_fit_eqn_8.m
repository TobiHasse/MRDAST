function [fits, coef] = survivor_fit_eqn_8(fx,ccd,jj)
% Purpose:  This function will fit survivor functions (CCDF) using equation
%           8 from Tobias Hasse's dissertation
%           It is reordered and cleaned up code from the files curve_fits.m
%           and curv_fits_array.m both originally written in December 2018.
%           
%           Dependancies:   input files created by compute_dists.m
%
%           Called by:      distribution_fits.m & distr_r_squared_ranges.m
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     edited October 2021, written December 2018
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
    % do the fitting
    [fobj,gof(i),out(i)] = fit(fx',ccd(i,:)',ft,'Startpoint',sp,...
        'lower',lo,'upper',up,'tolfun',1e-10);%,'maxfunevals',2000);

    fobj % optional, print the fit to screen
    rs(i,jj)=gof(i).rsquare;    % store r squared
    fits{i}=fobj;               % store the fit object in a cell array
    % use the struct coef to store useful info
    coef.names = coeffnames(fobj);
    coef.vals(i,:)=coeffvalues(fobj);
    temp = confint(fobj);       % confidence interval
    coef.intervalLo(i,:) = temp(1,:);
    coef.intervalHi(i,:) = temp(2,:);
    
    % show the results in a figure
    figure(99)
    plot(fobj,fx,ccd(i,:))
    ylim([10^-7 10^2])
    xlim([10 10^6])
    title(sprintf('survivor function %d of %d', i, size(ccd,1) ) )
    set(gca,'yscale','log',...
        'xscale','log')
    xlabel(sprintf('%s',out(i).message))
end %i

end % function survivor_fit_eqn_8

