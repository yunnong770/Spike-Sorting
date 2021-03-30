close all
clear
clc
raw_data = csvread('EMG_example_1_90s_fs_2k.csv');
% raw_data = csvread('EMG_example_2_fs_2k.csv');
% raw_data = csvread('EMG_example_20s_2000Hz-2020.csv');

%% Read data
emg1 = raw_data(:,1);
emg2 = raw_data(:,2);
figure(1)
plot(emg1);
figure(2)
plot(emg2);

Fs = 2000;
L = length(emg2);
T = 1/Fs;
t = (0:L-1)*T;
f = Fs*((-L/2):(L/2)-1)/L;
duration = (L-1)/Fs;

%% Filter
chan1.sig = emg2 - mean(emg2);
chan1.filt_sig = bandpass(chan1.sig, [20 200], Fs);
chan1.fft = fft(chan1.filt_sig);
d1 = designfilt('lowpassiir','FilterOrder',12, ...
'HalfPowerFrequency',0.15,'DesignMethod','butter');
chan1.filt_sig = filtfilt(d1, chan1.filt_sig);
figure(1);
plot(chan1.sig);
hold on
plot(chan1.filt_sig);

%% Baseline Removal
[c,l] = wavedec(chan1.filt_sig,4,'db8');
c(1:l(1)) = 0;
chan1.filt_br_sig = waverec(c, l , 'db8');
figure(1)
plot(chan1.filt_br_sig)
hold on
plot(chan1.filt_sig)
% chan1.filt_br_sig = chan1.filt_sig;
%% Thresholding
coef = 3;
sigma = median(abs(chan1.filt_br_sig)/0.6745);
thr = coef * sigma;
[peaks, loc] = findpeaks(abs(chan1.filt_br_sig), t, 'MinPeakHeight', thr, 'MinPeakDistance', 0.01);
loc_ind = int64(loc/duration*(L-1))+1;
peaks = chan1.filt_br_sig(loc_ind);
figure(3)
plot(t, chan1.filt_br_sig);
hold on
plot(loc, peaks, 'x');

%% Segmentation and Alignment
L_spike = 40;
N_spikes = length(loc);
spikes = zeros(N_spikes, L_spike);
for i = 1:N_spikes
    spikes(i,:) = chan1.filt_br_sig(loc_ind(i)-(L_spike/2-1):loc_ind(i)+L_spike/2);
end
% plot(spikes.');
%%
% Interpolation for better alignment
x = 0:T:T*(L_spike-1);
xq = 0:T/8:T*(L_spike-1);% Query point for interpolation
spikes_interp = interp1(x, spikes.', xq, 'spline');

% Alignment Adjustment
N = 120:1:180;
[M loc] = max(abs(spikes_interp(N,:)));
L_spike_interp = (L_spike-1)*8+1;
del_ind = []; % del possible incorrect thresholded spikes
for i = 1:N_spikes
   if (loc(i)<=10)
       del_ind = [del_ind i];       
   end
end
loc(del_ind) = [];
spikes_interp(:,del_ind) = [];
dim = size(spikes_interp);
N_spikes = dim(2);
spikes_aligned = zeros(250, N_spikes);
for i = 1:N_spikes
    spikes_aligned(:,i) = spikes_interp(loc(i)+120-1-124:loc(i)+120-1+125, i); 
end

% 2nd order derivative
% xq = xq(1:250);
% dy = diff(spikes_aligned);
% dt = diff(xq);
% dydt = (dy.')./dt;
% d_dydt = diff(dydt.');
% dt(end) = [];
% d2ydt2 = (d_dydt.')./dt;

%% Feature Extraction
feature = [];

% % 5 largest coefficient in wavelet decomposition
% wl_coef = [];
% ind = zeros(1,3);
% for i = 1:length(loc)
%     [cA cD] = dwt(spikes(i,:).','db4');
%     wl_coef = [wl_coef; [cA; cD].'];
% end
% V = var(wl_coef); % Obtain the variance of each wavelet coefficient
% for i = 1:3    
%     [M ind(i)] = max(V);
%     V(ind(i)) = 0;
%     feature = [feature spikes(:,ind(i))];
% end

% Width of the main spike

% % Find the width
% spikes_width = [];
% for i = 1:length(spikes)
%     n = 125;
%     r_ind = 0;
%     l_ind = 0;
%     while (spikes_interp(i,n)*spikes_interp(i,n+1)>0)
%         n = n+1;
%     end
%     if abs(spikes_interp(i,n)) > abs(spikes_interp(i,n+1))
%         r_ind = n+1;
%     else
%         r_ind = n;
%     end
% 
%     n = 125;
%     while (spikes_interp(i,n)*spikes_interp(i,n-1)>0)
%         n = n-1;
%     end
%     if abs(spikes_interp(i,n)) > abs(spikes_interp(i,n+1))
%         l_ind = n-1;
%     else
%         l_ind = n;
%     end
%     spikes_width = [spikes_width ; (r_ind - l_ind)*T];
% end

feature = [feature spikes_width];


%% Feature Normalization
colmin = min(feature);
colmax = max(feature);
feature = rescale(feature,'InputMin',colmin,'InputMax',colmax);

%% t-SNE
covariance_matrix = cov(feature);
[V, D] = eig(covariance_matrix);
dr_spikes = feature*V(:,39:end);
scatter3(dr_spikes(:,1).', dr_spikes(:,2).', dr_spikes(:,3).');