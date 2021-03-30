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
ca6 = appcoef(c, l, 'db4');
[cd1, cd2, cd3, cd4, cd5, cd6] = detcoef(c, l, [1 2 3 4 5 6]);
ca6(1:end) = 0;
cd1(cd1<max(cd1)) = 0;
cd2(cd2<max(cd2)) = 0;
cd3(cd3<max(cd3)) = 0;
cd4(cd4<max(cd4)) = 0;
cd5(cd5<=max(cd5)) = 0;
cd6(cd6<=max(cd6)) = 0;
c = [ca6; cd6; cd5; cd4; cd3; cd2; cd1];
chan1.filt_sig = waverec(c, l, 'db4'); % Reconstruct the signal with low frequency noise removed
% chan1.filt_sig = wden(chan1.filt_sig, 'minimaxi', 'h', 'one', 6, 'db4'); % Wavelet denoising
plot(chan1.filt_sig)

%% Detect Spikes
vNEO = chan1.filt_sig(2:end-1) .* chan1.filt_sig(2:end-1) - ...
    chan1.filt_sig(3:end).*chan1.filt_sig(1:end-2);
vNEO = [0; vNEO; 0];
[peaks_energy, locs_energy] = findpeaks(vNEO, 'MinPeakHeight', 0.03, 'MinPeakDistance', 20);
locs_spikes = [];
for i = 1:length(locs_energy)
    [M,loc] = max(abs(chan1.filt_sig((locs_energy(i)-5):(locs_energy(i)+5))));
    locs_spikes = [locs_spikes (loc + locs_energy(i)-6)];
end
plot(chan1.filt_sig);
hold on
plot(locs_spikes, chan1.filt_sig(locs_spikes), 'o')