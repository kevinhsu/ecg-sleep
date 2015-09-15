%% Database
db1 = '/slpdb/slp01a';
db2 = '/slpdb/slp02b';
db3 = '/slpdb/slp01b';
db4 = '/slpdb/slp04';

%% Configuration
db = db1;
[~,config]=wfdbloadlib;
echo on

%% Load ECG Signal
display('Reading samples ECG signal from MIT-BIH Arrhythmia Database')
[siginfo,Fs] = wfdbdesc(db);
Fs = Fs(1);
[tm, signal]=rdsamp(db,1);
ecg = signal(:,1);
ecg = cmddenoise(ecg,'db1',4); % Wavelet Shrinkage Denoise
ecg = ecg';
[N,~] = size(ecg);
figure
plot(tm,ecg);hold on

%% Peak Detection
display('Reading and plotting annotations (human labels) of QRS complexes performend on the signals')
sortecg = sort(ecg,'descend');
thresholdratio = 0.5;
threshold = thresholdratio * mean(sortecg(1:round(N/Fs)));  %Only one R-wave in one sampling point
[peaks, locs, widths, prominence] = findpeaks(ecg,tm,'MinPeakHeight',threshold);
plot(locs, peaks, 'ro');

%% Convert ECG to HRV
[m,n] = size(locs);
hrv = locs(2:m) - locs(1:(m-1));
figure
plot(locs(2:m), hrv)

%% Filter Outlier
[hrv,locs] = filterOutlier(hrv,locs);

%% Global Feature Extraction
rr.locs = locs;
rr.hrv = hrv;
rr.dethrv = detrend(rr.hrv);
rr.nolhrv = zscore(rr.hrv);
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
fitting = polyfit(log10(d),log10(p),1);
rr.dfaslope = fitting(1);
rr.dfaintercept = fitting(2);

% Piecewise Detrended Fluctuation Analysis
[rr.alpha1, rr.alpha2, rr.alpha3] = dfapiece(d,p);

% Theoretical Detrended Fluctuation Analysis
fitting = polyfit(log10(d),log10((d.^(0.5))/(d(1)^(0.5)/p(1))),1);
rr.dfaslopeT = fitting(1);
rr.dfainterceptT = fitting(2);

%% Rolling Time Window Feature Extraction
[ann,type,subtype,chan,num,comments]=rdann(db,'st',[],N);
[numTimeWindow,~] = size(comments);
timeWindow = floor(N ./ numTimeWindow);

for i = 1:numTimeWindow
    % Time Series
    for 
    % rrtw.hrv = rr.hrv((timeWindow * (i-1) + 1):(timeWindow * i));
    % rrtw.locs = rr.locs((timeWindow * (i-1) + 1):(timeWindow * i));
    rrtw.dethrv = detrend(rrtw.hrv);
    rrtw.nolhrv = zscore(rrtw.hrv);
    figure
    plot(rrtw.locs, rrtw.hrv);
    break;
    
    % Feature Extraction
    rrtw.mean(i) = mean(rrtw.hrv);
    rrtw.std(i) = std(rrtw.hrv);
    rrtw.CV(i) = rrtw.std / rrtw.mean;
    
    % Power Spectrum Density
    [rrtw.psd, rrtw.w] = periodogram(rrtw.dethrv, hamming(length(rrtw.dethrv)));
    [rrtw.LF(i), rrtw.HF(i), rrtw.FreqmaxP(i), rrtw.maxHFD(i), rrtw.LFHFratio(i), rrtw.inter(i)] = findLFHF(rrtw.psd, rrtw.w);
    
    % Detrended Fluctuation Analysis
    [H,pval95,d,p] = dfaecg(rrtw.dethrv);
    rrtw.H(i) = H;
    rrtw.pval951(i) = pval95(1);
    rrtw.pval952(i) = pval95(2);
    fitting = polyfit(log10(d),log10(p),1);
    rrtw.dfaslope(i) = fitting(1);
    rrtw.dfaintercept(i) = fitting(2);
    
    % Piecewise Detrended Fluctuation Analysis
    [rr.alpha1(i), rr.alpha2(i), rr.alpha3(i)] = dfapiece(d,p);
    
    % Theoretical Detrended Fluctuation Analysis
    fitting = polyfit(log10(d),log10((d.^(0.5))/(d(1)^(0.5)/p(1))),1);
    rrtw.dfaslopeT(i) = fitting(1);
    rrtw.dfainterceptT(i) = fitting(2);
end


%% Cross Validation
% Load Sleep Stage ann







