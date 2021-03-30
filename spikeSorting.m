close all
clear
clc
%% Input file
M= csvread('EMG_example_2_fs_2k.csv'); %read in csv file
time= M(:,1); % first column is the time series
fs= (time(2)-time(1))^-1; % calculate the sample frequency
channel_number= size(M,2)-1; % num of channels in the database
for i=1:channel_number,
figure('Color',[1 1 1]);plot(time,M(:,i+1)); %plot each channel
str= sprintf('Channel %d',i);
xlabel('seconds');title(str);xlim([time(1) time(size(time,1))]); % label and title each plots
end
channel_select= 1; % select channel for testing. channel_select<= channel_number
test_input= M(:,channel_select+1); % test_input will go through all the individual sections

%% Filter Signal
chan1.sig = test_input - mean(test_input); % DC removal
[c, l] = wavedec(chan1.sig, 6, 'db4'); % Remove low frequency noise by wavelet decomposition
c(1:l(1)) = 0; 
chan1.filt_br = waverec(c, l, 'db4'); % Reconstruct the signal with low frequency noise removed
chan1.filt_sig = wden(chan1.filt_br, 'minimaxi', 'h', 'one', 6, 'db4'); % Wavelet denoising
plot(chan1.filt_sig)

%% Detect Spikes
vNEO = chan1.filt_sig(2:end-1) .* chan1.filt_sig(2:end-1) - ...
    chan1.filt_sig(3:end).*chan1.filt_sig(1:end-2);
vNEO = [0; vNEO; 0];
[peaks_energy, locs_energy] = findpeaks(vNEO, 'MinPeakHeight', 0.03, 'MinPeakDistance', 18);
locs_spikes = [];
for i = 1:length(locs_energy)
    [M,loc] = max(abs(chan1.filt_sig((locs_energy(i)-5):(locs_energy(i)+5))));
    locs_spikes = [locs_spikes (loc + locs_energy(i)-6)];
end
plot(chan1.filt_sig);
hold on
plot(locs_spikes, chan1.filt_sig(locs_spikes), 'o')

% chan1.sig_analytical = abs(sqrt(chan1.filt_sig.^2 + hilbert(chan1.filt_sig).^2)); % Analytical Signal
% chan1.noise_var = var(chan1.sig_analytical);
% Z0 = linspace(0, max(chan1.sig_analytical), 1000);
% fz = 1/(2*chan1.noise_var)*exp(-Z0/(2*chan1.noise_var));
% plot(Z0/chan1.noise_var, fz)
% [p,x] = hist(chan1.sig_analytical, 150); plot(x,p/sum(p)); %PDF

%% Align Spikes
aligned_spikes = [];
for i = 1:length(locs_spikes)
    if locs_spikes(i)+20 <= length(chan1.filt_sig)
        append_spike = chan1.filt_sig(locs_spikes(i)-20:locs_spikes(i)+20);
        aligned_spikes = [aligned_spikes; append_spike.'];
    else
        append_spike = chan1.filt_sig(locs_spikes(i)-20:end);
        append_spike = [append_spike; zeros(41-length(append_spike), 1)];
        aligned_spikes = [aligned_spikes; append_spike.'];
    end
end
plot(aligned_spikes.')

%% Extract Features
colmin = min(aligned_spikes);
colmax = max(aligned_spikes);
feature = rescale(aligned_spikes,'InputMin',colmin,'InputMax',colmax);
plot(feature.')

%% Cluster spikes
Y = tsne(feature, 'Exaggeration',3, 'NumDimensions', 2, 'Perplexity', 20, 'LearnRate', 100);
% scatter(Y(:,1), Y(:,2));
s = [];
idx_mat = [];
for i = 3:15
    [idx,C,sumd,D] = kmeans(Y,i,'MaxIter',10000,...
        'Replicates',10);
    s = [s silhouette(Y,idx)];
    idx_mat = [idx_mat idx];
end
s = mean(s);
[M, num_group] = max(s);
idx = idx_mat(:,num_group);
num_group = num_group + 1;
figure(4)
for i = 1:num_group
    spike = mean(aligned_spikes(idx == i, :));
    plot(spike)
    hold on;
end
figure(5)
for i = 1:num_group
    scatter(Y(idx == i, 1), Y(idx == i, 2))
    hold on;
end

