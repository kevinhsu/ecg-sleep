%% Configuration
db = '/slpdb/slp01a';
display(db);
[~,config]=wfdbloadlib;
echo off
warning off

%% Load ECG Signal
display('Reading samples ECG signal from MIT-BIH Arrhythmia Database');
[siginfo,Fs] = wfdbdesc(db);
Fs = Fs(1); % Sampling Frequency
LengthSamples = siginfo.LengthSamples;

if LengthSamples > 1800000  % Avoid overflow. The maximum sampling lenght is 2 hours
    LengthSamples = 1800000;
end

[tm, signal] = rdsamp(db,[1],LengthSamples);
ecg = signal(:,1);
ecg = cmddenoise(ecg,'db1',4); % Wavelet Shrinkage Denoise
ecg = ecg';
N = length(ecg);
figure
plot(tm,ecg);hold on

%% Peak Detection
display('Reading and plotting annotations (human labels) of QRS complexes performend on the signals');
sortecg = sort(ecg,'descend');
thresholdratio = 0.5;
threshold = thresholdratio * mean(sortecg(1:round(N/Fs)));  %Only one R-wave in one sampling point
[peaks, locs, widths, prominence] = findpeaks(ecg,tm,'MinPeakHeight',threshold);
plot(locs, peaks, 'ro');

%% Convert ECG to HRV
display('Converting ECG to HRV');
m = length(locs);
hrv = locs(2:m) - locs(1:(m-1));
figure
plot(locs(2:m), hrv)

%% Filter Outlier
display('Filtering Outlier');
[hrv,locs] = filterOutlier(hrv,locs);
figure
plot(locs, hrv);hold on

% Smoothing Data (No Improvements)
% hrv = smooth(hrv);
% hrv = hrv';
% plot(locs, hrv);

%% Global Feature Extraction
display('Global Feature Extraction');
% Time Series
rr.locs = locs';
rr.hrv = hrv';
rr.dethrv = detrend(rr.hrv);
rr.nolhrv = zscore(rr.hrv);

% Statistical Feature Extraction
rr.mean = mean(rr.hrv);
rr.std = std(rr.hrv);
rr.CV = rr.std / rr.mean;

% Power Spectrum Density
[rr.psd, rr.w] = periodogram(rr.dethrv, hamming(length(rr.dethrv)));

figure
plot(rr.w, rr.psd);hold on
line([0.04 0.04], [0 .1], 'Color',[.8 .8 .8]);
line([0.15 0.15], [0 .1], 'Color',[.8 .8 .8]);
line([0.4 0.4], [0 .1], 'Color',[.8 .8 .8]);

[rr.LF, rr.HF, rr.FreqmaxP, rr.maxHFD, rr.LFHFratio, rr.inter] = findLFHF(rr.psd, rr.w);

% Detrended Fluctuation Analysis
[H,pval95,d,p] = dfaecg(rr.dethrv);
rr.H = H;
rr.pval951 = pval95(1);
rr.pval952 = pval95(2);
if ~isempty(d) && ~isempty(p)
    fitting = polyfit(log10(d),log10(p),1);
    rr.dfaslope = fitting(1);
    rr.dfaintercept = fitting(2);
end

% Theoretical Detrended Fluctuation Analysis
if ~isempty(d) && ~isempty(p)
    fitting = polyfit(log10(d),log10((d.^(0.5))/(d(1)^(0.5)/p(1))),1);
    rr.dfaslopeT = fitting(1);
    rr.dfainterceptT = fitting(2);
end

% Piecewise Detrended Fluctuation Analysis
if length(d) > 1
    [rr.alpha1, rr.alpha2, rr.alpha3, rr.alpha1flag, rr.alpha2flag, rr.alpha3flag] = dfapiece(d,p);
else
    rr.alpha1 = 0;
    rr.alpha2 = 0;
    rr.alpha3 = 0;
    rr.alpha1flag = 0;
    rr.alpha2flag = 0;
    rr.alpha3flag = 0;
end

% Fast Fourier Transformation
Yfft = fft(rr.hrv);
L = length(rr.hrv);
P2 = abs(Yfft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
[sortedP1,sortIndex] = sort(P1,'descend');

if numel(sortedP1) >= 1
    rr.P11 = sortedP1(1);  % Pick the Top 3 Elements
    rr.freq1 = f(sortIndex(1));
else
    rr.P11 = 0;
    rr.freq1 = 0;
end
if numel(sortedP1) >= 2
    rr.P12 = sortedP1(2);
    rr.freq2 = f(sortIndex(2));
else
    rr.P12 = 0;
    rr.freq2 = 0;
end
if numel(sortedP1) >= 3
    rr.P13 = sortedP1(3);
    rr.freq3 = f(sortIndex(3));
else
    rr.P13 = 0;
    rr.freq3 = 0;
end

%% Rolling Time Window Feature Extraction
display('Rolling Time Window Feature Extraction');
[ann,type,subtype,chan,num,comments]=rdann(db,'st',[],N);
numTimeWindow = length(comments);
timeWindow = floor(N ./ numTimeWindow);

for i = 1:numTimeWindow
    % Time Series and Time Window
    rrtw.hrv = [];
    rrtw.locs = [];
    for j = 1:length(rr.locs)
        if (rr.locs(j) >= tm(timeWindow * (i-1) + 1)) && (rr.locs(j) < tm(timeWindow * i)) % 30s, 60s, 90s, 120s
            rrtw.hrv(end + 1) = rr.hrv(j);
            rrtw.locs(end + 1) = rr.locs(j);
        end
    end
    
    rrtw.hrv = rrtw.hrv';
    rrtw.locs = rrtw.locs';
    
    % No Available Time Window
    if isempty(rrtw.hrv) || length(rrtw.hrv) == 1
        if i ~= 1
            rrtw.mean(i) = rrtw.mean(i-1);
            rrtw.std(i) = rrtw.std(i-1);
            rrtw.CV(i) = rrtw.CV(i-1);
            rrtw.LF(i) = rrtw.LF(i-1);
            rrtw.HF(i) = rrtw.HF(i-1);
            rrtw.FreqmaxP(i) = rrtw.FreqmaxP(i-1);
            rrtw.maxHFD(i) = rrtw.maxHFD(i-1);
            rrtw.LFHFratio(i) = rrtw.LFHFratio(i-1);
            rrtw.inter(i) = rrtw.inter(i-1);
            rrtw.H(i) = rrtw.H(i-1);
            rrtw.pval951(i) = rrtw.pval951(i-1);
            rrtw.pval952(i) = rrtw.pval952(i-1);
            rrtw.dfaslope(i) = rrtw.dfaslope(i-1);
            rrtw.dfaintercept(i) = rrtw.dfaintercept(i-1);
            rrtw.dfaslopeT(i) = rrtw.dfaslopeT(i-1);
            rrtw.dfainterceptT(i) = rrtw.dfainterceptT(i-1);
            rrtw.alpha1(i) = rrtw.alpha1(i-1);
            rrtw.alpha2(i) = rrtw.alpha2(i-1);
            rrtw.alpha3(i) = rrtw.alpha3(i-1);
            rrtw.alpha1flag(i) = rrtw.alpha1flag(i-1);
            rrtw.alpha2flag(i) = rrtw.alpha2flag(i-1);
            rrtw.alpha3flag(i) = rrtw.alpha3flag(i-1);
            rrtw.P11(i) = rrtw.P11(i-1);
            rrtw.P12(i) = rrtw.P12(i-1);
            rrtw.P13(i) = rrtw.P13(i-1);
            rrtw.freq1(i) = rrtw.freq1(i-1);
            rrtw.freq2(i) = rrtw.freq2(i-1);
            rrtw.freq3(i) = rrtw.freq3(i-1);
        else
            rrtw.mean(i) = 0;
            rrtw.std(i) = 0;
            rrtw.CV(i) = 0;
            rrtw.LF(i) = 0;
            rrtw.HF(i) = 0;
            rrtw.FreqmaxP(i) = 0;
            rrtw.maxHFD(i) = 0;
            rrtw.LFHFratio(i) = 0;
            rrtw.inter(i) = 0;
            rrtw.H(i) = 0;
            rrtw.pval951(i) = 0;
            rrtw.pval952(i) = 0;
            rrtw.dfaslope(i) = 0;
            rrtw.dfaintercept(i) = 0;
            rrtw.dfaslopeT(i) = 0;
            rrtw.dfainterceptT(i) = 0;
            rrtw.alpha1flag(i) = rrtw.alpha1flag(i-1);
            rrtw.alpha1(i) = 0;
            rrtw.alpha2(i) = 0;
            rrtw.alpha3(i) = 0;
            rrtw.alpha1flag(i) = 0;
            rrtw.alpha2flag(i) = 0;
            rrtw.alpha3flag(i) = 0;
            rrtw.P11(i) = 0;
            rrtw.P12(i) = 0;
            rrtw.P13(i) = 0;
            rrtw.freq1(i) = 0;
            rrtw.freq2(i) = 0;
            rrtw.freq3(i) = 0;
        end
        
        continue;
        
    end
    
    % Normalize Time Series
    rrtw.dethrv = detrend(rrtw.hrv);
    rrtw.nolhrv = zscore(rrtw.hrv);
    
    % Statistical Feature Extraction
    rrtw.mean(i) = mean(rrtw.hrv);
    rrtw.std(i) = std(rrtw.hrv);
    rrtw.CV(i) = rrtw.std(i) ./ rrtw.mean(i);
    
    % Power Spectrum Density
    [rrtw.psd, rrtw.w] = periodogram(rrtw.dethrv, hamming(length(rrtw.dethrv)));
    [rrtw.LF(i), rrtw.HF(i), rrtw.FreqmaxP(i), rrtw.maxHFD(i), rrtw.LFHFratio(i), rrtw.inter(i)] = findLFHF(rrtw.psd, rrtw.w);
    
    % Detrended Fluctuation Analysis
    [H,pval95,d,p] = dfaecg(rrtw.dethrv);
    rrtw.H(i) = H;
    rrtw.pval951(i) = pval95(1);
    rrtw.pval952(i) = pval95(2);
    if ~isempty(d) && ~isempty(p)
        fitting = polyfit(log10(d),log10(p),1);
        rrtw.dfaslope(i) = fitting(1);
        rrtw.dfaintercept(i) = fitting(2);
    else
        if i ~= 1
            rrtw.dfaslope(i) = rrtw.dfaslope(i-1);
            rrtw.dfaintercept(i) = rrtw.dfaintercept(i-1);
        else
            rrtw.dfaslope(i) = 0;
            rrtw.dfaintercept(i) = 0;
        end
    end
    
    % Theoretical Detrended Fluctuation Analysis
    if ~isempty(d) && ~isempty(p)
        fitting = polyfit(log10(d),log10((d.^(0.5))/(d(1)^(0.5)/p(1))),1);
        rrtw.dfaslopeT(i) = fitting(1);
        rrtw.dfainterceptT(i) = fitting(2);
    else
        if i ~= 1
            rrtw.dfaslopeT(i) = rrtw.dfaslopeT(i-1);
            rrtw.dfainterceptT(i) = rrtw.dfainterceptT(i-1);
        else
            rrtw.dfaslopeT(i) = 0;
            rrtw.dfainterceptT(i) = 0;
        end
    end
    
    % Piecewise Detrended Fluctuation Analysis
    if length(d) > 1
        [rrtw.alpha1(i), rrtw.alpha2(i), rrtw.alpha3(i), rrtw.alpha1flag(i), rrtw.alpha2flag(i), rrtw.alpha3flag(i)] = dfapiece(d,p);
    else
        rrtw.alpha1(i) = 0;
        rrtw.alpha2(i) = 0;
        rrtw.alpha3(i) = 0;
        rrtw.alpha1flag(i) = 0;
        rrtw.alpha2flag(i) = 0;
        rrtw.alpha3flag(i) = 0;
    end
    
    % Fast Fourier Transformation
    Yfft = fft(rrtw.hrv);
    L = length(rrtw.hrv);
    P2 = abs(Yfft/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    [sortedP1,sortIndex] = sort(P1,'descend');
    
    if numel(sortedP1) >= 1
        rrtw.P11(i) = sortedP1(1);  % Pick the Top 3 Elements
        rrtw.freq1(i) = f(sortIndex(1));
    else
        rrtw.P11(i) = 0;
        rrtw.freq1(i) = 0;
    end
    if numel(sortedP1) >= 2
        rrtw.P12(i) = sortedP1(2);
        rrtw.freq2(i) = f(sortIndex(2));
    else
        rrtw.P12(i) = 0;
        rrtw.freq2(i) = 0;
    end
    if numel(sortedP1) >= 3
        rrtw.P13(i) = sortedP1(3);
        rrtw.freq3(i) = f(sortIndex(3));
    else
        rrtw.P13(i) = 0;
        rrtw.freq3(i) = 0;
    end
    
end

% Clear rrtw Signal
rrtw.hrv = [];
rrtw.locs = [];
rrtw.dethrv = [];
rrtw.nolhrv = [];
rrtw.psd = [];
rrtw.w = [];

% Clear rr Signal
rr.hrv = [];
rr.locs = [];
rr.dethrv = [];
rr.nolhrv = [];
rr.psd = [];
rr.w = [];


%% Label Extraction
display('Label Extraction');
% Load Sleep Stage Annotations
for i = 1:length(comments)
    C = char(comments{i});
    switch C(1)
        case 'W'
            S = -1;  % Subject is awake
        case 'R'
            S = 0;   % REM sleep
        case '1'
            S = 1;   % Sleep stage 1
        case '2'
            S = 2;   % Sleep stage 2
        case '3'
            S = 3;   % Sleep stage 3
        case '4'
            S = 4;   % Sleep stage 4
        otherwise
            S = -1;  % Other status
    end
    
    rrtw.labels(i) = S;
    
    switch C(1)
        case 'W'
            S = -1;  % Subject is awake
        case 'R'
            S = 1;   % REM sleep
        case '1'
            S = 1;   % Sleep stage 1
        case '2'
            S = 1;   % Sleep stage 2
        case '3'
            S = 1;   % Sleep stage 3
        case '4'
            S = 1;   % Sleep stage 4
        otherwise
            S = -1;  % Other status
    end
    
    rrtw.binarylabels(i) = S;
end

rrtw.labels =  rrtw.labels';
rrtw.binarylabels = rrtw.binarylabels';

%% K-folds Cross Validation
display('k-folds Cross Validation');
globalfeatures = [];
features = [];
labels = [];
binarylabels = [];

% Global Features: Used to match the most similar patient in database
globalfeatures(1:numTimeWindow,end+1) = rr.mean';
globalfeatures(1:numTimeWindow,end+1) = rr.std';
globalfeatures(1:numTimeWindow,end+1) = rr.CV';
globalfeatures(1:numTimeWindow,end+1) = rr.LF';
globalfeatures(1:numTimeWindow,end+1) = rr.HF';
globalfeatures(1:numTimeWindow,end+1) = rr.FreqmaxP';
globalfeatures(1:numTimeWindow,end+1) = rr.maxHFD';
globalfeatures(1:numTimeWindow,end+1) = rr.LFHFratio';
globalfeatures(1:numTimeWindow,end+1) = rr.inter';
globalfeatures(1:numTimeWindow,end+1) = rr.H';
globalfeatures(1:numTimeWindow,end+1) = rr.pval951';
globalfeatures(1:numTimeWindow,end+1) = rr.pval952';
globalfeatures(1:numTimeWindow,end+1) = rr.dfaslope';
globalfeatures(1:numTimeWindow,end+1) = rr.dfaintercept';
globalfeatures(1:numTimeWindow,end+1) = rr.dfaslopeT';
globalfeatures(1:numTimeWindow,end+1) = rr.dfainterceptT';
globalfeatures(1:numTimeWindow,end+1) = rr.alpha1';
globalfeatures(1:numTimeWindow,end+1) = rr.alpha2';
globalfeatures(1:numTimeWindow,end+1) = rr.alpha3';
% globalfeatures(1:numTimeWindow,end+1) = rr.alpha1flag';
% globalfeatures(1:numTimeWindow,end+1) = rr.alpha2flag';
% globalfeatures(1:numTimeWindow,end+1) = rr.alpha3flag';
globalfeatures(1:numTimeWindow,end+1) = rr.P11';
globalfeatures(1:numTimeWindow,end+1) = rr.P12';
globalfeatures(1:numTimeWindow,end+1) = rr.P13';
globalfeatures(1:numTimeWindow,end+1) = rr.freq1';
globalfeatures(1:numTimeWindow,end+1) = rr.freq2';
globalfeatures(1:numTimeWindow,end+1) = rr.freq3';
% Please add more features here ...

% Windowed Feautures: Used to classify sleep stages of selected patient
features(:,end+1) = rrtw.mean';
features(:,end+1) = rrtw.std';
features(:,end+1) = rrtw.CV';
features(:,end+1) = rrtw.LF';
features(:,end+1) = rrtw.HF';
features(:,end+1) = rrtw.FreqmaxP';
features(:,end+1) = rrtw.maxHFD';
features(:,end+1) = rrtw.LFHFratio';
features(:,end+1) = rrtw.inter';
features(:,end+1) = rrtw.H';
features(:,end+1) = rrtw.pval951';
features(:,end+1) = rrtw.pval952';
features(:,end+1) = rrtw.dfaslope';
features(:,end+1) = rrtw.dfaintercept';
features(:,end+1) = rrtw.dfaslopeT';
features(:,end+1) = rrtw.dfainterceptT';
features(:,end+1) = rrtw.alpha1';
features(:,end+1) = rrtw.alpha2';
features(:,end+1) = rrtw.alpha3';
% features(:,end+1) = rrtw.alpha1flag';
% features(:,end+1) = rrtw.alpha2flag';
% features(:,end+1) = rrtw.alpha3flag';
features(:,end+1) = rrtw.P11';
features(:,end+1) = rrtw.P12';
features(:,end+1) = rrtw.P13';
features(:,end+1) = rrtw.freq1';
features(:,end+1) = rrtw.freq2';
features(:,end+1) = rrtw.freq3';
% Please add more features here ...

% Load binary class labels and multi-class labels
labels = rrtw.labels';
binarylabels = rrtw.binarylabels';

% Write Features and Labels in .txt files
save 'globlefeatures.txt' globalfeatures -ascii
save 'features.txt' features -ascii
save 'labels.txt' labels -ascii
save 'binarylabels.txt' binarylabels -ascii;

% Normailzation
features = mapminmax(features', 0, 1);
features = features';

% Superivised Learning
svmbinarycl = fitcsvm(features,binarylabels,'KernelFunction','rbf','BoxConstraint',Inf,'ClassNames',[-1,1]);
svmbinarycl2 = crossval(svmbinarycl);
binarymisclass = kfoldLoss(svmbinarycl2,'mode','individual');   % 10-folds cross validation
svmbinarymisclass = mean(binarymisclass);
display(svmbinarymisclass);

svmmulticl = fitcsvm(features,labels,'KernelFunction','rbf','BoxConstraint',Inf,'ClassNames',[-1,1]);
svmmulticl2 = crossval(svmmulticl);
multimisclass = kfoldLoss(svmmulticl2,'mode','individual');   % 10-folds cross validation
svmmultimisclass = mean(multimisclass);
display(svmmultimisclass);

% Unsupervised Learning for 2 Clusteriing
[u centroid] = kmeans(features, 2);
u = u';
temp = unique(u);
count = 0;

for i = 1:length(u)
    if u(i) == temp(1)
        u(i) = 1;
    else
        u(i) = -1;
    end
    if u(i) == binarylabels(i)
        count = count + 1;
    end
end

kmeansbinarymisclass1 = (length(u) - count) ./ length(u);
count = 0;

for i = 1:length(u)
    if u(i) == temp(1)
        u(i) = -1;
    else
        u(i) = 1;
    end
    if u(i) == binarylabels(i)
        count = count + 1;
    end
end

kmeansbinarymisclass2 = (length(u) - count) ./ length(u);
kmeansbinarymisclass = min(kmeansbinarymisclass1, kmeansbinarymisclass2);
display(kmeansbinarymisclass);






