%
%
%simulate threshold, convolution and NEO
%detection of spikes
%
%define time
Fs = 10000;
dt = 1/Fs;
t = 0:dt:1;

%define APs
AP1amp = 1;
AP1dur1 = .001;
AP1dur2 = .002;
AP1prob = .001;
AP2amp = AP1amp*.5;

%make single AP
%by concantenating two 1/2 sin cycles
oneT = 0:dt:(AP1dur1+AP1dur2);
oneAP = AP1amp*sin(2*pi/(2*AP1dur1)*oneT).*(oneT<=AP1dur1) + ...
    AP1amp*0.5*sin(2*pi/(2*AP1dur2)*(oneT+AP1dur1)).*(oneT>AP1dur1) ;

%make a bunch of APs
locationAP1 = rand(size(t))<AP1prob; %poisson
%locationAP1 = mod([1:length(t)],457)==0; %uniform
vAP = conv(locationAP1,oneAP);
vAP = vAP(1:length(t)); %truncate
locationAP2 = rand(size(t))<AP1prob; %poisson
vAP2 = conv(locationAP2,oneAP*AP2amp);
vAP = vAP + vAP2(1:length(t)); %truncate

%make gaussian noise plus 60 Hz interference
noiseAmp = .3;
interAmp = .2;
vnoise = noiseAmp*randn(size(t)) + interAmp*sin(60*2*pi*t);

%now add niose to APs and 
%filter it the way most physiologists do
vtotal = vnoise + vAP ;
vin=vtotal;

tic

[b,a] = butter(4, 1000/(Fs/2)); %4 KHz, 4-pole lowpass
vtotal = filter(b,a,vtotal);
[b,a] = butter(2, 200/(Fs/2), 'high'); %200 Hz, 2-pole highpass
vtotal = filter(b,a,vtotal);

NumPlots = 3;
clf;

%The final waveforms
subplot(NumPlots,1,1)
title('raw data (butterworth filtered) and threshold')
hold on
plot(t,vtotal)  %the signal
plot(t,vAP*max(vtotal)*.5+max(vtotal)*1.5,'r') %the computed APs
plot(t,vin*max(vtotal)*1-max(vtotal)*1.5,'k')
%tLow=.5*max(vtotal);
%plot(t, ((vtotal>tLow))*max(vtotal)*.5-max(vtotal)*1.5,'k' ); %threshold of raw signal

%Matched template convolution 
subplot(NumPlots,1,2)
%oneHat = 0:dt:(AP1dur1+AP1dur2);
%oneHat = AP1amp*(oneT<=AP1dur1) - AP1amp*0.5*(oneT>AP1dur1) ;
%vMex = conv(vtotal,oneHat)/length(oneHat);
vMex = conv(vtotal,oneAP)/length(oneAP);
%recenter the convolution and truncate excess
vMex = vMex(length(oneAP)/2:length(vMex)-length(oneAP)/2);
plot(t,vMex); %the convolved signalx
hold on
plot(t,vAP*max(vMex)*.5+max(vMex)*1.5,'r') %the computed APs
tLow = 0.5*max(abs(vMex));
plot(t, (abs(vMex)>tLow)*max(vMex)*.5-max(vMex)*1.5,'k' ); %threshold of convolved
title('matched template')

subplot(NumPlots,1,3)
title('matched -> NEO')
%redo NEO using the output from the convolution filter
vNEO = vMex(2:end-1) .* vMex(2:end-1) - ...
    vMex(3:end).*vMex(1:end-2);
vNEO = [0 vNEO 0];
w = bartlett(length(oneAP)/2);
vNEO = conv(vNEO,w);
vNEO = vNEO(length(w)/2:length(vNEO)-length(w)/2);
hold on
plot(t,vNEO)
plot(t,vAP*max(vNEO)*.5+max(vNEO)*1.5,'r')
tLow = 0.33*max(vNEO);
%plot(t, ((vNEO>tLow)&(vNEO<tHigh))*max(vNEO)*.5-max(vNEO)*1.3,'g' );
plot(t, (vNEO>tLow)*max(vNEO)*.5-max(vNEO)*1.5,'k' );

toc





%