%% Database names
db1 = '/slpdb/slp01a';
db2 = '/slpdb/slp02b';
db3 = '/slpdb/slp01b';
db4 = '/slpdb/slp04';

%% Configuration
db = db1;
[~,config]=wfdbloadlib;
echo on

%% Load ECG
display('Reading samples ECG signal from MIT-BIH Arrhythmia Database')
[siginfo,Fs] = wfdbdesc(db);
Fs = Fs(1);
[tm, signal]=rdsamp(db,1);
ecg = signal(:,1);
[N,~] = size(ecg);
figure
plot(tm,ecg);hold on;grid on;

%% Peak Detection
display('Reading and plotting annotations (human labels) of QRS complexes performend on the signals')
sortecg = sort(ecg,'descend');
thresholdratio = 0.5;
threshold = thresholdratio * mean(sortecg(1:round(N/Fs)));  %Only one R-wave in one sampling point
[peaks, locs, widths, prominence] = findpeaks(ecg,tm,'MinPeakHeight',threshold);
plot(locs, peaks, 'ro');

%% ECG to HRV
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

% Spectral features - Power Spectrum Density
[rr.psd, rr.w] = periodogram(rr.dethrv, hamming(length(rr.dethrv)));

figure
plot(rr.w, rr.psd);hold on;
line([0.04 0.04], [0 .1], 'Color',[.8 .8 .8]);
line([0.15 0.15], [0 .1], 'Color',[.8 .8 .8]);
line([0.4 0.4], [0 .1], 'Color',[.8 .8 .8]);

[rr.LF, rr.HF, rr.FreqmaxP, rr.maxHFD, rr.LFHFratio, rr.inter] = findLFHF(rr.psd, rr.w);



%% Rolling time window Feature Extraction
rrtw;





% x = rr.nolhrv;
% N = length(x);
% xdft = fft(x);
% xdft = xdft(1:N/2+1);
% psdx = (1/(Fs*N)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:Fs/length(x):Fs/2;
% figure
% plot(freq,10*log10(psdx))
% grid on
% title('Periodogram Using FFT')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% Power spectrum density
% LF
% PF


%%Sleep stage ann
[ann,type,subtype,chan,num,comments]=rdann(db,'st',[],N);






