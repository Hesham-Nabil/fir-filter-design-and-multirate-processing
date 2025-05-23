function q_subband = quantize_subband(subband, levels)
    q_subband = subband;
    thresholds = (levels(1:end-1) + levels(2:end)) / 2;
    for i = 1:length(thresholds)
        q_subband(subband >= thresholds(i)) = levels(i+1);
    end
end