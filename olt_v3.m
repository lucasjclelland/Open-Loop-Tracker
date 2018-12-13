clear;

%close all;

FSamp = 5e6; % Sampling rate Hz
blockSize = floor(FSamp*.001); % Samples per block
tBlockSize = blockSize/FSamp; % Seconds
simTime = 2; % Seconds
numBlocks = floor(FSamp*simTime/blockSize); % Number of simulated blocks

%% Generate carrier
t=[0:1:blockSize*numBlocks-1]/FSamp; % Seconds
fCarrier = rand*10000-5000; % Hz
% fCarrier=333.333;
pDeltaCarr=rem(2*pi*fCarrier*t,2*pi);
carrier = exp(1i*pDeltaCarr);
carrier = awgn(carrier,100,'measured'); % Add noise

%% Simulate initial acquisition
fGuess=fCarrier+1000*rand()-500; % Initial guess within +-500 Hz
fGuess=fCarrier-5;
% [~,loc]=max(fftshift(fft(carrier)));
% fGuess=-FSamp/2+loc/length(carrier)*FSamp;
% fScale = linspace(-FSamp/2,FSamp/2,length(carrier));
fprintf('Initial error: %d(Hz)\n',fGuess-fCarrier);

%% Open loop carrier tracking

counter = 1;
fErrors = zeros(numBlocks,1);
fGuesses = ones(numBlocks,1)*fGuess;

while numBlocks>2
    
    counter=counter+1;
    
    I = 0;
    Q = 0;
    
    for idx=1:numBlocks
        
        % fGuess=fGuesses(idx);
        pDeltaGuess=rem(2*pi*fGuesses(idx)*t((idx-1)*blockSize+1:idx*blockSize),2*pi); %remainder of exponent divided by 2pi
        rxSignal = carrier((idx-1)*blockSize+1:idx*blockSize);
        localReplica = exp(1i*pDeltaGuess);
        I_d1 = I;
        Q_d1 = Q;
        tmp = localReplica*rxSignal';
        I=real(tmp);
        Q=imag(tmp);
        
        if idx>1
            fErrors(idx)=fnFreqDiscrim(I_d1,Q_d1,I,Q,tBlockSize,1)/(2*pi); % Freq error estimate Hz
        end
        
    end
    
    error_coeffs=polyfit((2:1:numBlocks).'*tBlockSize,fErrors(2:end),2);  % use least squares to find rate of change (rad/s)
    guess_coeffs=polyfit((1:1:numBlocks).'*tBlockSize,fGuesses(1:end),3); % fGuess coefficients
    % plot results
    figure (1);
%     plot(fGuesses-fCarrier); 
    plot(fErrors(2:end),'x')
    hold on
    errorfit=polyval(error_coeffs,(2:1:numBlocks).'*tBlockSize);
    plot(errorfit)
    %xlim([0 200])
    LineWidth=2;
    ylabel('Frequency Error (Hz)')
    xlabel('Blocks')
    title('Frequency Errors Fit to Curve')
    hold off
%     t2(counter)=mean(fGuesses-fCarrier);
    drawnow;
    % end plot results
    blockSize=blockSize*2; %double the blocksize each run
    numBlocks = floor(FSamp*simTime/blockSize); % Number of simulated blocks
    tBlockSize = blockSize/FSamp; % Seconds
    pause(.5)
    fErrors = polyval(error_coeffs,(1:1:numBlocks).'*tBlockSize); %Generate half as many guesses based on a best fit curve
    fGuesses = polyval(guess_coeffs,(1:1:numBlocks).'*tBlockSize);
    fGuesses= fGuesses-fErrors;
    fErrors = zeros(numBlocks,1);
    counter=counter+1;
    
end