%% Find LF and HF
function [LF, HF, FreqmaxP, maxHFD, LFHFratio, inter] = findLFHF(psd, w)
VLFpsd = [];
VLFw = [];
LFpsd = [];
LFw = [];
HFpsd = [];
HFw = [];
[m,~] = size(w);

for i = 1:m
    if w(i) <= 0.05
        VLFpsd(end+1) = psd(i);
        VLFw(end+1) = w(i);
    end
    if w(i) > 0.05 && w(i) <= 0.15
        LFpsd(end+1) = psd(i);
        LFw(end+1) = w(i);
    end
    if w(i) > 0.15 && w(i) <= 0.4
        HFpsd(end+1) = psd(i);
        HFw(end+1) = w(i);
    end
end

LF = trapz(LFw, LFpsd) ./ (trapz(w, psd) - trapz(VLFw, VLFpsd));
HF = trapz(HFw, HFpsd) ./ (trapz(w, psd) - trapz(VLFw, VLFpsd));
LFHFratio = LF ./ HF;
inter = LF ./ (LF + HF);
[maxHFD, maxIndex] = max(HFpsd);
FreqmaxP = HFw(maxIndex);



