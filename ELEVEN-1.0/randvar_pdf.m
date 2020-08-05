function [bins,freq,area] = randvar_pdf(data,numbins,opt)

    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 3
        error('Too many inputs.')
    elseif nargin == 2
        opt = 0;
    elseif  nargin == 3
        xlow = prctile(data,0.5*(100 - 95));
        xupp = prctile(data,0.5*(100 + 95));
    end

    Ns = length(data);

    data_max = max(data);
    data_min = min(data);
    
    if opt == 0
        binwidth = (data_max-data_min)/(numbins-1);
        bins     = (data_min:binwidth:data_max);
    elseif opt == 1
        binwidth = (xupp-xlow)/(numbins-1);
        bins = (xlow:binwidth:xupp);
    elseif opt == 2
        binwidth = (xupp-xlow)/(numbins-1);
        bins = [data_min xlow:binwidth:xupp];
    else
        binwidth = (xupp-xlow)/(numbins-1);
        bins = [xlow:binwidth:xupp data_max];
    end

freq     = histc(data,bins);
freq     = freq/(Ns*binwidth);
area     = binwidth*sum(freq);

return