function [hrv,locs] = filterOutlier(hrv,locs)
rr.hrv = hrv;
rr.mean = mean(hrv);
rr.std = std(hrv);
[m,~] = size(hrv);
numstd = 1.96;  % 95% confidence interval
hrvtemp = [];
locstemp = [];

for i = 1:m
    if (hrv(i) >= rr.mean - numstd * rr.std) && (hrv(i) <= rr.mean + numstd * rr.std)
        hrvtemp(end+1) = hrv(i);
        locstemp(end+1) = locs(i);
    end
end

hrv = hrvtemp;
locs = locstemp;
figure
plot(locs, hrv)