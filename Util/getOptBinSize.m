function [dxOpt] = getOptBinSize(data,type)
% Function `getOptBinSize' returns the optimal number of bins in a histogram
% based on a set of given criterium
%
% Example usage:
% dxOpt = getOptBinSize(data,'Sturges'); hist(data,0:dxOpt:max(data));
%
% Inputs
% data: data vector
% type: type of data
%   options: Shimizaki: only assumes independence (great for poisson)
%            Sturges: assumes ~normal
%            Doane: Sturges with skewedness (better for non-normal)
%            Scott: Best for normal data
%            Freedman-Diaconis: Based on IQR = less sensitive to outliers
%
% Outputs
% dxOpt: optimal bin size

% 2021 Liangyu Tao

if ~isempty(data)
    n = numel(data);
    switch type
        case 'Shimizaki'
            n2cons = ceil(range(data)./2)-1;
            for k = 2:n2cons+1
                for it = 1:51
                    nBin = k-1;
                    shift = range(data)./nBin.*(it-1)./50;
                    bins = linspace(min(data)+shift-range(data)./nBin,max(data)+shift-range(data)./nBin,k);
                    tmp = histcounts(data,bins);
                    mu = mean(tmp);
                    var = sum((tmp-mu).^2)./(nBin);
                    bSize = diff(bins(1:2));
                    c(it,k) = (2*mu-var)./(bSize.^2);
                end
                bSizeAll(k) = bSize;
                
            end
            dxOpt = bSizeAll(mean(c)==min(mean(c)));
            
        case 'Sturges'
            k = log2(n)+1;
            dxOpt = range(data)./k;
        case 'Doane'
            g1 = skewness(data);
            sig = sqrt(6.*(n-2)./(n+1)./(n+3));
            k = 1+log2(n)+log2(1+abs(g1)./sig);
            dxOpt = range(data)./k;
        case 'Scott'
            dxOpt = 3.49*std(data)./(nthroot(n,3));
        case 'Freedman–Diaconis'
            dxOpt = 2*iqr(data)./(nthroot(n,3));
            
    end
else
    dxOpt = [];
end

end