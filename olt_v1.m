clear; close all;
fSamp = 5e6;                            % Hz
blockSize=1024;                        % Samples
blockTime = 1/fSamp*blockSize;          % Seconds
simBlocks = 30;                         % Number of Simulated blocks
t=[0:1:(simBlocks*blockSize)-1]/fSamp;  % Seconds

%% Simulate carrier
carrier_freq = rand*10000-5000;          % Hz
carrier = exp(1i*2*pi*carrier_freq*t);  

%% Simulate channel
rxSignal=awgn(carrier,100,'measured');   % simulate AWGN Channel
%rxSignal = rxSignal.*exp(1i*pi*rand);   % add random phase


%% Estimate carrier using FFT
tmp = abs(fftshift(fft(rxSignal(1:blockSize))));
fScale = linspace(-fSamp/2,fSamp/2,blockSize);
[~,loc]=max(tmp);
fGuess=fScale(loc);

fprintf('Initial frequency error: %d\n',carrier_freq-fGuess);
%fGuess=carrier_freq;

%% Open Loop Carrier Tracking
fGuesses=zeros(simBlocks,1);
fErrors=zeros(simBlocks,1);
for idx=1:simBlocks

   localReplica = exp(1i*2*pi*fGuess*t(((idx-1)*simBlocks)+1:idx*simBlocks));
   rxBlock = rxSignal(((idx-1)*simBlocks)+1:idx*simBlocks);
   fGuesses(idx)=fGuess;
   xcorrVals(idx)=dot(localReplica,rxBlock);

   if idx>1
     dotVal = real(xcorrVals(idx-1))*real(xcorrVals(idx))+...
              imag(xcorrVals(idx-1))*imag(xcorrVals(idx));
     crossVal = real(xcorrVals(idx-1))*imag(xcorrVals(idx))-...
              imag(xcorrVals(idx-1))*real(xcorrVals(idx));
     fError = atan2(dotVal,crossVal)/(blockSize/fSamp);
     fError = (atan2(imag(xcorrVals(idx)),real(xcorrVals(idx))) - atan2(imag(xcorrVals(idx-1)),real(xcorrVals(idx-1))))/(blockSize/fSamp);
     fError = (atan2(imag(xcorrVals(idx)),real(xcorrVals(idx))))/(blockSize/fSamp);
%      mError=1e-15;
%      if fError<mError
%          fError=0;
%      end
     fErrors(idx)=fError;
     fGuess=fGuess+fError;
   end

end

figure(1)
plot(carrier_freq-fGuesses)