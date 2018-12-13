% function df = fnFreqDiscrim(I_d1, Q_d1, I, Q, T, type)
% Frequency discriminator as defined in Phil Ward Chapter of
% Kaplan
%
% Arguments:
% I_d1 - Previous in-phase correlator output
% Q_d1 - Previous quadraphase correlator output
% I    - Current in-phase correlator output
% Q    - Current quadraphase correlator output
% T    - Pre-detection integration time [sec]
% type - atan2, atan
%
% Outputs
% df - frequency estimate in radians/sec

% EE607 - Navigation Receiver Design
% Summer 2009
%
%
% Copyright 2009 - Sanjeev Gunawardena
%
% August 20, 2009 - Initial Version
% May 27, 2018 - added type as an input argument


function DiscrFreq = fnFreqDiscrim(I_d1, Q_d1, I, Q, T, type)
switch type
  case 0
    DiscrFreq=0;
  case 1
    DiscrFreq=(atan2((I_d1*Q)-(Q_d1*I),(I_d1*I)+(Q_d1*Q)))/T;
  case 2
    DiscrFreq=(atan(((I_d1*Q)-(Q_d1*I))/((I_d1*I)+(Q_d1*Q))))/T;
    if isnan(DiscrFreq)
      DiscrFreq=0;
    end
  case 3
    DiscrFreq=(atan(((I_d1*I)+(Q_d1*Q)))/((I_d1*Q)-(Q_d1*I)))/T;
    if isnan(DiscrFreq)
      DiscrFreq=0;
    end
  otherwise
    error('unrecognized type: %d',type)
end
end
