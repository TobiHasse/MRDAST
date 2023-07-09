function [PDF CDF MED] = pdfcdfpctile(sd,dt,pctile)
    % Purpose    Return the Probability density functino, Cumulative
    %            distribution function and a percentile, such as the median
    %            for an input array assumed to be a lower triangular matrix
    %   INPUT    sd   a lower triangular matrix of storage quantities with
    %                 model age on rows and sediment age on columns
    %            dt   model step
    %            pctile the percentile value sought 
    %   OUTPUT   PDF  a probability distribution map
    %            CDF  a cumulative distribution map
    %            MED  a vector of median age
    % Author:    Tobias Hasse tobiack@udel.edu
    % Date:      edited October 2021, written November 2018

    if pctile>1
        disp('percentile is greater than 1, dividing by 100')
        pctile = pctile / 100;
    end
    % make a probability density function (PDF) of the distribution
    PDF=bsxfun(@rdivide,sd,sum(sd,2));
    
    % make a cumulative distribution function
    CDF=cumsum(PDF,2);
    % if pctile is 50% then this will find the median
    [jnk sd_median]=min(abs(CDF-pctile),[],2);    % index of the percentile
    MED = sd_median * dt;             % bc each index = dt years
end

