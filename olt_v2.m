clear;

%close all;

FSamp = 5e6; % Sampling rate Hz
blockSize = floor(FSamp*.001); % Samples per block
tBlockSize = blockSize/FSamp; % Seconds
simTime = 2; % Seconds
numBlocks = floor(FSamp*simTime/blockSize); % Number of simulated blocks
%% Generate carrier
t=[0:1:blockSize*numBlocks-1]/FSamp; % Seconds
%fCarrier = rand*10000-5000; % Hz
fCarrier=2000;
pDeltaCarr=rem(2*pi*fCarrier*t,2*pi); %exponent
carrier = exp(1i*pDeltaCarr); %no noise
carrier = awgn(carrier,-40,'measured'); %add noise
%load test;
%carrier=test;
%% Initial guess using FFT
% [~,loc]=max(fftshift(fft(carrier)));
% fGuess=-FSamp/2+loc/length(carrier)*FSamp;

% [~,loc]=max(fftshift(fft(carrier,2^23))); %Unrealistic fft resolution
% fGuess=-FSamp/2+loc/(2^23)*FSamp;
fGuess=fCarrier-200;

%fCarrier=fGuess;

%fGuess=fCarrier+50;

fprintf('Initial guess: %d(Hz)\n',fGuess);

%% Open loop carrier tracking

counter = 1;

while numBlocks>2
    angles = zeros(numBlocks,1);
    I=0;
    Q=0;
    pDeltaGuess=rem(2*pi*fGuess*t,2*pi);
    for idx=1:numBlocks
        rxSignal = carrier((idx-1)*blockSize+1:idx*blockSize);
        localReplica = exp(1i*pDeltaGuess((idx-1)*blockSize+1:idx*blockSize));
        
        I_d1 = I;
        Q_d1 = Q;
        
        tmp = localReplica*rxSignal';
        
        I=real(tmp);
        Q=imag(tmp);
        
        fError=0;
        if idx>1
            fError=fnFreqDiscrim(I_d1,Q_d1,I,Q,tBlockSize,1); %atan2
            %fError=fnFreqDiscrim(I_d1,Q_d1,I,Q,tBlockSize,2); %atan
        end
        
        %angles(idx)=atan2(Q,I);
        
        angles(idx)=fError; %Estimated for each block 
        
    end
    
    %errorfit=polyfit((1:1:numBlocks).'*tBlockSize,unwrap(angles),1); % use least squares to find rate of change (rad/s)
    
    meanError=mean(angles(2:numBlocks)); %Mean of angles
    
    fError = meanError(1)/(2*pi); % convert from rad/s
    fErrors(counter)=fError;
    
    fGuess=fGuess-fError;
    fGuesses(counter)=fGuess;
    
    % plot results
    figure (1);
    subplot(4,1,1);
    plot((1:1:counter),fGuesses);
    hold('on');
    xlabel('Simulations');
    ylabel('Frequency (Hz)');
    plot((1:1:counter),ones(counter,1)*fCarrier);
    legend('Estimated','Actual');
    hold('off');
    
    subplot(4,1,2);
    plot((1:1:counter),fGuesses-fCarrier);
    xlabel('Simulations');
    ylabel('FError (Hz)');
    
    subplot(4,1,3);
    plot((1:1:numBlocks).'*tBlockSize,unwrap(angles));
    hold('on');
    yFit=polyval(meanError,(1:1:numBlocks).'*tBlockSize);
    plot((1:1:numBlocks).'*tBlockSize,yFit);
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    legend('Measured','Fit');
    hold('off');
    
    subplot(4,1,4)
    [pxx,f]=periodogram(carrier,[],blockSize,FSamp,'centered');
    plot(f,10*log10(pxx))
    xlabel('dB/Hz');
    ylabel('Hz');
    axis([-10e3,10e3,-80,10]);
    drawnow;
    % end plot results
    counter=counter+1;
    blockSize=blockSize*2;
    numBlocks = floor(FSamp*simTime/blockSize); % Number of simulated blocks
    tBlockSize = blockSize/FSamp; % Seconds
    pause(.5)
    
end

fprintf('Final guess: %d(Hz)\n',fGuess);