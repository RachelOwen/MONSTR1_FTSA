function [Spec_Freq, FreqMin, FreqMax] = resampleWL2Freq(Spec_WL, WLMin, WLMax, NewLength)
% resampleWL2Freq resamples spectrum taken with spectrometer with wavelength axis 
% to make it evenly spaced in frequency domain, along the row
%    Spec_WL: Spectral data vector as a function of wavelength
%    WLMin, WLMax: the first and last wavelength
%    NewLength: the number of data points after resampling
%    Spec_Freq: Spectral data vector evenly spaced in frequency domain with
%               ascending frequency axis
%    FreqMin, FreqMax: the first and last frequency
                
% Wavelength: nm
% Frequency: THz
% If center WL set to 804nm, WLMin = 791.7008, WLMax = 814.5770

% Constant: Light speed in air (Unit: nm*Thz) 
%LightSpeedC = 2.99702547E+5;
LightSpeedC = 2.99709E+5;
RowSize = size(Spec_WL, 1);
ColSize = size(Spec_WL, 2);
WLIndex = linspace(WLMin, WLMax, RowSize);
FreqIndex = LightSpeedC ./ WLIndex; 
FreqMax = FreqIndex(1);
FreqMin = FreqIndex(end);
FreqIndexNew = linspace(FreqMax, FreqMin, NewLength);

% the resampling happens along the row for each column
Spec_Freq_temp = interp1( FreqIndex', Spec_WL, FreqIndexNew', 'linear');
% Reverse spectrum to have ascending frequency axis
Spec_Freq = Spec_Freq_temp(end:(-1):1, :);
