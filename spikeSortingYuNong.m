close all
clear
clc
%% Input file
filename = 'EMG_example_1_90s_fs_2k.csv'; % Modify the filename to the desired file
channel_select= 1; % select channel for testing. channel_select<= channel_number
if filename == "EMG_example_1_90s_fs_2k.csv" || filename == "EMG_example_20s_2000Hz-2020.csv"
    M= csvread(filename);   
    fs = 2000;
    time= (1:length(M)).'/2000;
    channel_number= size(M,2); % num of channels in the database
    test_input= M(:,channel_select); % test_input will go through all the individual sections
     for i=1:channel_number
        figure('Color',[1 1 1]);plot(time,M(:,i)); %plot each channel
        str= sprintf('Channel %d',i);
        xlabel('seconds');title(str);xlim([time(1) time(size(time,1))]); % label and title each plots
    end
else
    M= csvread(filename); %read in csv file
    time= M(:,1); % first column is the time series
    fs= (time(2)-time(1))^-1; % calculate the sample frequency
    channel_number= size(M,2)-1; % num of channels in the database
    for i=1:channel_number
        figure('Color',[1 1 1]);plot(time,M(:,i+1)); %plot each channel
        str= sprintf('Channel %d',i);
        xlabel('seconds');title(str);xlim([time(1) time(size(time,1))]); % label and title each plots
    end
    test_input= M(:,channel_select+1); % test_input will go through all the individual sections
end

%% Filter Signal
close all
chan.sig = test_input - mean(test_input); % DC removal 
[c, l] = wavedec(chan.sig, 6, 'db4'); % Remove low frequency noise by wavelet decomposition
c(1:l(1)) = 0; % Setting approximate coefficients to zero
chan.filt_br = waverec(c, l, 'db4'); % Reconstruct the signal with low frequency noise removed
chan.filt_sig = wdenoise(chan.filt_br, 5,'DenoisingMethod', 'Minimax', 'Wavelet',...
    'db4', 'ThresholdRule', 'Hard', 'NoiseEstimate', 'LevelDependent'); % Wavelet denoising in detailed levels

% Plot raw vs filtered signal
figure
subplot(2,1,1)
plot(time, chan.filt_sig, 'r')
title("Filtered Signal with Wavelet Denoising")
xlabel("Time (s)")
ylabel("Voltage (mV)")
xlim([0, max(time)])
subplot(2,1,2)
plot(time, chan.sig)
title('Raw Signal')
xlabel("Time (s)")
ylabel("Voltage (mV)")
xlim([0, max(time)])
%% Detect Spikes
vNEO = chan.filt_sig(2:end-1) .* chan.filt_sig(2:end-1) - ...% Calculate NEO of the filtered signal
    chan.filt_sig(3:end).*chan.filt_sig(1:end-2);
vNEO = [0; vNEO; 0]; % zero padding to make the NEO same length as the filtered signal
rms_vNEO = sqrt(sum(vNEO.^2)/length(vNEO)); % Calculate RMS of NEO for thresholding
[peaks_energy, locs_energy] = findpeaks(vNEO,...
    'MinPeakHeight', 0.5*rms_vNEO, 'MinPeakDistance', 15); % Find peaks of NEO that is greater than the threshold

% Adjust the center of the spikes in case the peak of NEO is different than
% the peak of the spikes
locs_spikes = [];
for i = 1:length(locs_energy)
    if (locs_energy(i)-5) >= 1 && (locs_energy(i)+5) <= length(chan.filt_sig) 
        [M,loc] = max(abs(chan.filt_sig((locs_energy(i)-5):(locs_energy(i)+5))));
    elseif (locs_energy(i)-5) < 1
        [M,loc] = max(abs(chan.filt_sig(1:(locs_energy(i)+5))));% Handle the case where a spike is detected at the very beginning of the recording
    else
        [M,loc] = max(abs(chan.filt_sig((locs_energy(i)-5):end)));% Handle the case where a spike is detected at the end of the recording
    end
    locs_spikes = [locs_spikes (loc + locs_energy(i)-6)]; % Append the location of newly found spike to the end of the array that stores the location of every spike
end
figure
subplot(2,1,1)
plot(time, vNEO);
hold on
plot(time, ones(1, length(vNEO))*rms_vNEO, '--')
title('NEO vs Threshold')
xlabel("Time (s)")
ylabel("Energy (mW)")
xlim([0, max(time)])
ylim([0, rms_vNEO*10])
subplot(2,1,2)
plot(time, chan.filt_sig);
hold on
plot(locs_spikes/fs, chan.filt_sig(locs_spikes), 'o')
title('Identified Spikes')
xlabel("Time (s)")
ylabel("Voltage (mV)")
xlim([0, max(time)])

%% Align Spikes
close all

aligned_spikes = []; % Matrix that stores every spikes
% Loop through the location of the center of every spikes and store them in
% the matrix above
for i = 1:length(locs_spikes) 
    if locs_spikes(i)-20 < 1 % Handles the case where a spike appears at the very beginning of the signal
        append_spike = zeros(21-locs_spikes(i), 1); % Zero padding to make the length identical
        append_spike = [append_spike; chan.filt_sig(1:locs_spikes(i)+20)];
        aligned_spikes = [aligned_spikes; append_spike.']; 
    elseif locs_spikes(i)+20 <= length(chan.filt_sig) % Handles the normal case where zero padding is not needed
        append_spike = chan.filt_sig(locs_spikes(i)-20:locs_spikes(i)+20);
        aligned_spikes = [aligned_spikes; append_spike.'];
    else
        append_spike = chan.filt_sig(locs_spikes(i)-20:end);% Handles the case where a spike appears at the very end of the signal
        append_spike = [append_spike; zeros(41-length(append_spike), 1)];% Zero padding to make the length identical
        aligned_spikes = [aligned_spikes; append_spike.'];
    end
end
plot(aligned_spikes.')
title('Aligned Spikes')
xlabel("Samples")
ylabel("Voltage (mV)")
%% Extract Features
close all

feature = aligned_spikes(:,15:27); % Use only the 13 samples near the center as features since they are most relevant
% Perform T-SNE with small perplexity and large exaggeration to increase
% separability
Y = tsne(feature, 'Exaggeration', 10, 'NumDimensions', 3, 'Perplexity', 4, 'LearnRate', 100);

% Plot the feature space reduced by T-SNE
figure
scatter3(Y(:, 1), Y(:, 2), Y(:, 3))
title("Feature Space under T-SNE")
%% Cluster spikes
close all

Z = linkage(Y,'single'); % Create the hierarchical cluster tree for single-linkage clustering
idx = cluster(Z,'maxclust',round(length(aligned_spikes)/60)); % Perform clustering with adaptive maximum number of clusters

% Plot the clustering result in the feature space
for i = 1:max(idx)
    scatter3(Y(idx == i, 1), Y(idx == i, 2), Y(idx == i, 3))
    hold on;
end
title("Clustered Result in T-SNE space")
figure
hist(idx, round(length(aligned_spikes)/60))
title("Histogram of the Clustering Result")
figure
dendrogram(Z, 'ColorThreshold', 3)
title("Dendrogram of the Hierarchical Cluster Structure")
%% Classify Spikes
close all

% Apply adaptive thresholding to filter clusters with small number of
% samples
[idx_count,]=hist(idx,unique(idx)); % Get the number of samples in each cluster
spike_group = (idx_count >= length(aligned_spikes)/30 ).*(1:length(idx_count)); % Apply thresholding
spike_group(spike_group == 0) = []; % Remove group numbers that are filtered out

% Plot MUAPs and the corresponding MUAPTs
MUAPs = zeros(length(idx_count), 26); % Matrix that holds MUAPs
MUAPTs = zeros(length(idx_count),length(chan.sig));% Matrix that holds MUAPTs
for i = spike_group
    ind_spikes = find(idx == i);
        MUAPs(i,:) = mean(aligned_spikes(ind_spikes,8:33), 1);
    for j = 1:length(ind_spikes) 
       center = locs_spikes(ind_spikes(j));
       if (center - 13 < 1) % Handle the case where a spike appears at the very beginning of the signal
           MUAPTs(i, 14-center:center + 12) = MUAPs(i,27-center:end);
       elseif (center + 12 > length(time))% Handle the case where a spike appears at the very end of the signal
           MUAPTs(i, center - 13:end) = MUAPs(i,1:14+(length(time)-center));
       else
           MUAPTs(i, center - 13:center + 12) = MUAPs(i,:);
       end
    end
    
    % Plot each MUAPT with respect to the filtered signal
    figure
    plot(time,chan.filt_sig)
    hold on
    plot(time, MUAPTs(i,:)', 'r')
    title("MUAPT" + i + " vs Filtered Signal");
    xlabel("Time (s)")
    ylabel("Voltage (mV)")
    xlim([0, max(time)])
end

% Plot all resolved spikes
figure
plot(MUAPs(spike_group,:).')
title("All resolved spikes")
xlabel("Samples")
ylabel("Voltage (mV)")
 
%% Analysis
close all
num_neurons = length(spike_group);

% Firing rate of each neuron
bin_width_second = 0.1;% 0.1 second for each bin (for both Gaussian and linear window)
bin_width_sample = bin_width_second * fs; % For the optional binning method to calculate the firing rate over time
num_bins = (length(time) - 1) / fs / bin_width_second;

% Obtain the location of each spike for different neurons
spikes_locs_binary = zeros(num_neurons, length(chan.sig)); % matrix stores firing or not at a time stamp by binary coding
firing_rate = zeros(num_neurons,length(spikes_locs_binary)); % matrix stores firing rate
for i = 1:num_neurons
    ind_spikes = find(idx == spike_group(i)); % Get the index of spikes belong to a certain group
    for j = 1:length(ind_spikes) 
        center = locs_spikes(ind_spikes(j));
        spikes_locs_binary(i, center) = 1;
    end 
end
window = gausswin(bin_width_sample); % Create Gaussian window
for i = 1:num_neurons
    firing_rate(i,:) = conv(spikes_locs_binary(i,:).', window, 'same')/bin_width_second; % Perform convolution to get the frequency
%     Optional binning method to calculate the firing rate over time      
%     for j = 1:num_bins
%         firing_rate(i,j) = sum(window.'.*spikes_locs_binary(i, 1+(j-1)*bin_width_sample:j*bin_width_sample)...
%             /bin_width_second);
%     end

    % Plot the firing rate over time for each neuron
    figure
    subplot(2,1,1)
    plot(time, MUAPTs(spike_group(i),:), 'r')
    title("Firing Rate Vs Time for Neuron " + i)
    xlabel("Time (s)")
    ylabel("Voltage (mV)")
    subplot(2,1,2)
    plot(time, firing_rate(i,:))
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
end

% Interspike Interval Distribution
interspike_interval = zeros(num_neurons, length(aligned_spikes)); % Create redundancy in the length to acommondate different number of spikes
for i = 1:num_neurons
    spikes_locs = find(spikes_locs_binary(i,:) == 1);
    num_spikes = sum(spikes_locs_binary(i,:));
    interspike_interval(i, 1:num_spikes-1) = spikes_locs(2:end)-spikes_locs(1:end-1); % Subtract the location of a preceeding spike from the next spike
end
interspike_interval = interspike_interval/fs/1000; % Convert the unit from samples to milliseconds

% Plot the interspike interval distribution for each spike
for i = 1:num_neurons
    spikes_locs = find(spikes_locs_binary(i,:) == 1);
    num_spikes = sum(spikes_locs_binary(i,:));
    interval = interspike_interval(i, 1:num_spikes-1);
    figure
    hist(interval, 3000)
    title("Distribution of Inter-spike Interval of Neuron " + i)
    xlabel("Time (s)")
    ylabel("Count")
%     xlim([0 0.001]) % Uncomment for a zoom-in view
end

% Phase locking (Pick the first channel of the recording below as it 
% appears to be the result of a periodic stimulation)

% Only perform phase locking analysis to a certain signal since it is only 
% meaningful when periodic pattern exists
if (filename == "EMG_example_2_fs_2k.csv") && channel_select == 1
    figure
    most_spikes_group = mode(idx); % Perform phase locking analysis on the most frequently appearing spike
    truncated_signal = zeros(22,1301);
    % Manual Truncation since the periodicity is not perfect for automated
    % trucation
    truncated_signal(1,:) = MUAPTs(most_spikes_group, 1400:2700);
    truncated_signal(2,:) = MUAPTs(most_spikes_group, 3200:4500);
    truncated_signal(3,:) = MUAPTs(most_spikes_group, 4900:6200);
    truncated_signal(4,:) = MUAPTs(most_spikes_group, 6700:8000);
    truncated_signal(6,:) = MUAPTs(most_spikes_group, 8500:9800);
    truncated_signal(7,:) = MUAPTs(most_spikes_group, 10200:11500);
    truncated_signal(8,:) = MUAPTs(most_spikes_group, 12000:13300);
    truncated_signal(9,:) = MUAPTs(most_spikes_group, 13800:15100);
    truncated_signal(10,:) = MUAPTs(most_spikes_group, 15600:16900);
    truncated_signal(11,:) = MUAPTs(most_spikes_group, 17400:18700);
    truncated_signal(12,:) = MUAPTs(most_spikes_group, 19300:20600);
    truncated_signal(13,:) = MUAPTs(most_spikes_group, 21000:22300);
    truncated_signal(14,:) = MUAPTs(most_spikes_group, 22700:24000);
    truncated_signal(15,:) = MUAPTs(most_spikes_group, 24500:25800);
    truncated_signal(16,:) = MUAPTs(most_spikes_group, 26300:27600);
    truncated_signal(17,:) = MUAPTs(most_spikes_group, 28100:29400);
    truncated_signal(18,:) = MUAPTs(most_spikes_group, 30000:31300);
    truncated_signal(19,:) = MUAPTs(most_spikes_group, 31700:33000);
    truncated_signal(20,:) = MUAPTs(most_spikes_group, 33600:34900);
    truncated_signal(21,:) = MUAPTs(most_spikes_group, 35400:36700);
    truncated_signal(22,:) = MUAPTs(most_spikes_group, 37500:38800);
    
    % Plot the average fft
    fft_truncated = abs(fft(truncated_signal));
    freqHz = (0:1:length(fft_truncated)-1)*fs/length(truncated_signal);
    freqHz = freqHz - max(freqHz)/2;
    plot(freqHz, fftshift(mean(fft_truncated)))
    xlabel("Frequency (Hz)")
    ylabel("Amplitude")
    title("Average FFT")
end

