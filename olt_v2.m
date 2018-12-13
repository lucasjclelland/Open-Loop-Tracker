clear; 
%close all;
FSamp = 5e6;                                    % Sampling rate Hz
blockSize = floor(FSamp*.001);                  % Samples per block 
tBlockSize = blockSize/FSamp;                   % Seconds
simTime = .2;                                   % Seconds
numBlocks = floor(FSamp*simTime/blockSize);     % Number of simulated blocks

%% Generate carrier
t=[0:1:blockSize*numBlocks-1]/FSamp;            % Seconds
fCarrier = rand*10000-5000;                     % Hz
pDeltaCarr=rem(2*pi*fCarrier*t,2*pi);
carrier = exp(1i*pDeltaCarr);      
carrier = awgn(carrier,-40,'measured');

%% Initial guess using FFT
[~,loc]=max(fftshift(fft(carrier)));
fGuess=-FSamp/2+loc/length(carrier)*FSamp;
fprintf('Initial error: %d(Hz)\n',fGuess-fCarrier);

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
       I_d1 = I; %Previous Block I 
       Q_d1 = Q; %Previous Block Q
       tmp = localReplica*rxSignal'; %correlate local with rx
       I=real(tmp);
       Q=imag(tmp);
       if idx>1
            angles(idx)=fnFreqDiscrim(I_d1,Q_d1,I,Q,tBlockSize,1)/(2*pi); %How does this change anything?
       end
      % angles(idx)=atan2(Q,I); %Are these the phase angles?
    end

    tmp=polyfit((1:1:numBlocks).'*tBlockSize,unwrap(angles),1); % use least squares to find rate of change (rad/s)
    fError = tmp(1)/(2*pi);
    fErrors(counter)=fError;
    fGuess=fGuess-fError;
    fGuesses(counter)=fGuess;
    
    %   plot results
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
    yFit=polyval(tmp,(1:1:numBlocks).'*tBlockSize);
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
    numBlocks = floor(FSamp*simTime/blockSize);          % Number of simulated blocks
    tBlockSize = blockSize/FSamp;                       % Seconds
    pause(.5)
end
fprintf('Final error: %d(Hz)\n',fGuess-fCarrier);