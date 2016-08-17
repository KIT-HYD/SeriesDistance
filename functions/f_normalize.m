function [k] = f_normalize(k)
% function normalizes the range of the input vector k to the range [0 1]
    k = (k-min(k(:)))/(max(k(:))-min(k(:)));
    % return zeros if all are NaN (e.g. relevant for e_sd_t if time series are shifted horizontally
    if(all(isnan(k))) 
        k = zeros(length(k),1);
    end
end
