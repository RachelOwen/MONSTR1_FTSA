% 2018-02-27 CLS, version 3:
%   Merged FTSA and FTSA_Multi2DNew so that there will be less ambiguity as programs continue to be updated.
%   In shorter version and 1D scans, draws data from beginning of data set, not separately collected FWM signals.
%   On 2017-02-27, Eric added option in FTSA_v2f (for 2D data sets, not yet multi2D) to not resample t axis in evenly spaced points
%   Plots use frequency axis for f_t of unevenly sampled points directly
%   converted from lambda. This is now incorporated.
%   Note from CLS: When not resampling, there might be an issue with
%   correcting t0 and correcting GDD, since these parameters are parameters
%   are corrected by multiplying by a frequency-dependent phase. Something
%   to fix in the future?
% 2016-10-14 CLS, version 2e:
%   Previously didn't work with SiV sample... Fixed so that AbsFreqAxis is monotonic (see near line 574).
%   Also fixed bugs related to wrap-arounds.
% 2016-05-30 CLS, version 2d:
%   Fixed inconsistency with absorption axis limits between S1 and S2. (Could still use some cleanup.)
% 2016-03-10 CLS, version 2c:
%   Fixed bug with making tausteps finite.
%   Added option to generate only the first T-step.
% 2016-03-02 CLS, version 2b:
%   Eliminates the use of redundant "chopped" matrices.
% 2016-03-02 CLS, version 2a:
%   Modified for 3D from version 2D standard FTSA.
%   Works both for T-dependent 3D data and Correlation spectra.

format long g;   clear all;   clc;   tic;
%profile on;
HowProcessBtn = questdlg('Retrieve one spectrum or process whole data set?',...
    'Retrieve one or Process all', 'One', 'One 2D spectrum', 'Whole', 'Whole');
%if strcmp(HowProcessBtn, 'Cancel')|isempty(HowProcessBtn) return; end;
if isempty(HowProcessBtn) return; end;
%-------------------------------------------------------------------------
% constants and control parameters:
SpeedC = 2.99709E+5;       % nm / ps, speed of light in air (from Wolfram alpha)
HeNeFreq = 473.61338;          % HeNe frequency THz
Planck =   4.13566766;           % Planck's constant "h", meV/THz.
hc = Planck*SpeedC;  % h*c in eV nm

NEmiDim = 1024;  NAbsDim = 1024; %NTDim = 1024; %Length of the emission dimension (CCD length) and the absorption dimension (zero padding included). By construction, NAbsDim cannot get smaller than NAbsDim0 (see line 451).
isPrepulseOrDelayT = 'T'; %Changing prepulse power (P) or delay (T)
is2Dor3D = '3D';
ScanIndex = '16'; %Index number for 2D scans
UdrSmplRatio = 32; % The ratio of under-sampling, no under-sampling = 1
isRephasing = 1; % Determine whether this scan is rephasing (1) or non-rephasing (0)
isSRDTrephasing = 0; % Determine to do 2D-rephasing with SRDT spectrum (1) or not (0). If other value, then forces phase to this value.

ReferenceSubtractionAmplitude = 0.0;
ScatteringSubtractionAmplitude = 0.0;
ScatteringSubtractiont0 = 0.2; %ps;
ScatteringSubtractionPhase = 0.333; % radians/pi
PPfname = 'DT_data_A.dat';

% Select type of plots
Abs = 1;
Real = 0;

WLMin = 825.797;  WLMax = 852.914;    % nm. % Wavelengths of SPEX-750M spectrometer
EneLowLmt = ceil(hc/WLMax);  EneHiLmt = floor(hc/WLMin); %meV %Settings used when plotting 2D spectra: Plot range and contour levels. %Defaults unless you identify a smaller range in the next line.
%EneLowLmt = 1527;  EneHiLmt = 1545; 
FreqLowLmt = EneLowLmt/Planck; FreqHiLmt = EneHiLmt/Planck; %THz
NContourLevels = 20;
PhiSI = 0; PhiCam = 0; PhiTR = 0;  %radians. %Phaser factors for SI, Camera and tracer
TauSteps = 0; % Use all tau-steps spectra (0); provide number of steps otherwise.

Second = 1; % Second time cut (1) or not (0)?
TimeCut = 5.7;    % ps %Cut off residual part near time zero after IFT FWM
TimeCut2 = 12;  % ps.

Delay_t0 = 6.30;   %ps   % delay between third pulse and reference
isGDDcorrection = 1; % Determine whether to correct for GDD in the substrate (1) or not (0).
GDD = (9000)*(2)/10^6; % fs^2/mm * mm, with an extra factor to convert to ps^2 (ultimately compare to THz). From http://refractiveindex.info/?shelf=main&book=GaAs, GVD of GaAs is 9000.3 at 970 nm.
GDD0 = 1475.7/Planck; % Center Frequency (meV)/Planck's constant will give THz.

% Follow echo (1) or not (0)? ("Second" must also be selected above)
Echo = 1;
c = 60;     % Start following offset after 'c' tau steps
m = 0.2;    % Slope for following echo
d = 100;     % Number of points for each tau step

% Offset correction parameter
isOffsetCor = 0;    % 0 for no offsect correction, and yes windowing.
TukeyAlpha = 0.25;     % Select a decimal between 0 (no window) and 1 (Hanning window).
TauIndex = 0;       % 0 for last point; step number otherwise

% Save the images as pngs/pdfs? Execution will take longer.
isSaveImages = 0;

% Deconvolution parameters
isDeconvolve = 0;   % 0 for no deconvolution
SpexRes = 40;       % HWHM in ueV
SpexG = 2*pi*SpexRes/(1000*Planck);   % Angular width (x 10^12 rad/s)

isPhaseCycling = 1; % Determine whether scan has used 2-axis phase cycling (1) or not (0).  If 1-axis phase cycling is used, set to (0) also.
isResample = 1; % Resample WL to frequency?

% Nyquist frequency, would normally include a factor of 2 in the
% denominator, but NO FACTOR OF 2 here: a full interferometer fringe is half fringe at OC
NyquistFreq =  HeNeFreq / UdrSmplRatio;
% 2D matrix elments in the window defined by IndexLowLmt&IndexHiLmt numbers or
%FreqLowLmt&FreqHiLmt area plotted while those confined by NFirstCol&NLastCol numbers saved
IndexLowLmt = round(((FreqLowLmt/SpeedC-1/WLMax)*NEmiDim)/(1/WLMin-1/WLMax));
IndexHiLmt = round(((FreqHiLmt/SpeedC-1/WLMax)*NEmiDim)/(1/WLMin-1/WLMax));
Multiplier = 1.0 - 1.0;
NFirstCol = round(IndexLowLmt - Multiplier*(IndexHiLmt-IndexLowLmt)/2.0);
NLastCol = round(IndexHiLmt + Multiplier*(IndexHiLmt-IndexLowLmt)/2.0);
if (NFirstCol < 1) NFirstCol = 1; end
if (NLastCol > NEmiDim) NLastCol = NEmiDim; end

%% Data path and Filenames:

DataPath = strcat('./',is2Dor3D,num2str(str2num(ScanIndex)),'/');
FWMfname = 'FWM00.dat';
% SIfname  = 'SI00.dat'; %Obsolete, but can be used if desired not to pull phase-cycled data See PHSCL flags below.
Reffname = 'REF00.dat';
TIFWMfname = 'TI-FWM_tau.dat';
if isRephasing
    FName2Ddata = '2dmatrixS1';
    %FName2DFFT = '2dfftmatrixS1';
    FName2DFFT = 'OutputS1';
else
    FName2Ddata = '2dmatrixS2';
    %FName2DFFT = '2dfftmatrixS2';
    FName2DFFT = 'OutputS2';
end
if is2Dor3D == '2D'
    SIfname = 'SI0total1.dat';
elseif is2Dor3D == '3D'
    SIfname = strcat(FName2Ddata,'/SItotal',isPrepulseOrDelayT,'0.dat');
end

%NOTE: if using extension tube, comment out flipud(filename) operation!
% Also change this operation for matrix of all spectra (starting at line
% 214).

FWMSpec1 = dlmread(strcat(DataPath, FWMfname), '\t');
FWMSpec1 = flipud(FWMSpec1); % flipud(array) flips array up-down
RefSpec = dlmread(strcat(DataPath, Reffname), '\t');
RefSpec = flipud(RefSpec);
% Read full 2D matrix.
MatIntf1 = dlmread(strcat(DataPath, SIfname), '\t');
MatIntf1 = flipud(MatIntf1); % See near line 128(?) for when to comment this out!
SISpec = MatIntf1(:,1);

%% 
%------------------------------------------------------------------------
OutDataPath = strcat(DataPath, FName2DFFT, '/');
if ~isdir(OutDataPath) mkdir(DataPath, FName2DFFT); end
% change reference strength in analysis (RefRatio = 1.00 no change)
RefRatio = 1.0;
FWMSpec = FWMSpec1 - sum([FWMSpec1(1:10); FWMSpec1(1015:end)])/20; % subtract dark counts ave.
% FWMSpec = FWMSpec1; % uncomment if you are not going to use line above
%InterfSpec = SISpec - FWMSpec - RefRatio * RefSpec; % Use this instead of below if not pulling phase-cycled data. PHSCLPHSCL
InterfSpec = SISpec;

if isResample
    [RefSpec_THz, FreqMin, FreqMax] = resampleWL2Freq(RefSpec, WLMin, WLMax, NEmiDim);
    [FWMSpec_THz, FreqMin, FreqMax] = resampleWL2Freq(FWMSpec, WLMin, WLMax, NEmiDim);
    [SISpec_THz, FreqMin, FreqMax] = resampleWL2Freq(SISpec, WLMin, WLMax, NEmiDim);
    [InterfSpec_THz, FreqMin, FreqMax] = resampleWL2Freq(InterfSpec, WLMin, WLMax, NEmiDim);
    % create frequency axis
    FreqAxis = transpose(linspace(FreqMin, FreqMax, NEmiDim));
    dlmwrite( strcat(OutDataPath, 'FWM_Ref_InterfSpec.dat'),...
        [FreqAxis FWMSpec_THz RefSpec_THz InterfSpec_THz], '\t');    
else
    FreqMin = SpeedC / WLMax;
    FreqMax = SpeedC / WLMin;
    RefSpec_THz = RefSpec(end:(-1):1, :);
    FWMSpec_THz = FWMSpec(end:(-1):1, :);
    SISpec_THz = SISpec(end:(-1):1, :);
    InterfSpec_THz = InterfSpec(end:(-1):1, :);
    % create frequency axis
    RowSize = size(RefSpec, 1);
    WLIndex = linspace(WLMax, WLMin, RowSize);
    FreqAxis = SpeedC./WLIndex';
end

% Create time axis. If isResample is false, this axis is not quite right, but it doesn't matter: 
RTimeAxis = linspace( 0, (NEmiDim-1)/(FreqMax-FreqMin), NEmiDim );
EneAxis = FreqAxis*Planck;

% Clean up reference spectrum
Ref_tDom = fft(RefSpec_THz);
Ref_tDom(60:end-58) = 0;
RefSpec_THz = ifft(Ref_tDom);
[MaxRef,MaxRefIndex] = max(RefSpec_THz);
FreqCenter = FreqAxis(MaxRefIndex);
EneCenter = Planck*FreqCenter;

figure(100);
set(gcf, 'Units', 'inch');
% set(gcf, 'position', [ 0 5.75 4 4 ]);
set(gcf, 'position', [ 0 5.75 3 2.75 ]);
subplot(211);
semilogy(RTimeAxis,abs(Ref_tDom));
subplot(212);
plot(EneAxis,RefSpec_THz);
dlmwrite( strcat(OutDataPath, 'Ref_Filtered.dat'),[EneAxis RefSpec_THz/max(RefSpec_THz)], '\t');

if isSRDTrephasing ~= 0
    % Pump-probe data path and filename:
    PPDataPath = DataPath;
    PPSpec  = dlmread(strcat(PPDataPath, PPfname), '\t');
    PPSpec = flipud(PPSpec); %Original 3D version wasn't flipped (?), so make sure we get this right before trusting. TKTK
    if (isResample)
        [DTSpec, FreqMin, FreqMax] = resampleWL2Freq(PPSpec(:, 1), WLMin, WLMax, NEmiDim);
        [TrSpec, FreqMin, FreqMax] = resampleWL2Freq(PPSpec(:, 2), WLMin, WLMax, NEmiDim);
    else
        DTSpect = PPSpec(:,1);
        DTSpec = DTSpec(end:(-1):1, :);
        TrSpec = PPSpec(:,2);
        TrSpec = TrSpec(end:(-1):1, :);
    end
    dlmwrite( strcat(OutDataPath, 'PumpProbe.dat'),[EneAxis DTSpec], '\t');
        
    DTSpec = DTSpec - ReferenceSubtractionAmplitude*max(abs(DTSpec))*RefSpec_THz/max(RefSpec_THz) ...
        - ScatteringSubtractionAmplitude*max(abs(DTSpec))*RefSpec_THz/max(RefSpec_THz) ...
        .* cos(2*pi*ScatteringSubtractiont0*FreqAxis - pi*ScatteringSubtractionPhase);
        
    %Subtract off the average of first 100 points, then scale to max.
    DTSpec = DTSpec - mean(DTSpec(1:100));
    DTSpec = DTSpec ./ max(DTSpec);
    
    %     TrSpec = TrSpec - sum([TrSpec(1:10); TrSpec(1015:end)])/20;
    TrSpec = TrSpec ./(max(TrSpec));
end

%% Plot the spectra of a single FWM, Ref and Interferogram

fig1 = figure(1);
set(gcf, 'Units', 'inch');
% set(gcf, 'position', [ 4 5.75 4 4 ]);
set(gcf, 'position', [ 3.25 5.75 3 2.75 ]);
subplot(211);   % input spectra
plot(EneAxis, RefSpec_THz, 'b-', EneAxis, FWMSpec_THz*10.0, 'k-', EneAxis, SISpec_THz, 'c-');
AxisScale = axis;    axis([ EneLowLmt-1 EneHiLmt+1 AxisScale(3:4) ]);
xlabel('Emission Frequency (meV)', 'FontSize', 10);
ylabel('Spectral Intensity (a.u.)', 'FontSize', 10);
title('Spectra of FWM(x10), Reference and their interferogram', 'FontSize', 10);
subplot(212);   % determined interference term
plot(EneAxis, InterfSpec_THz, 'c-');
AxisScale = axis;    axis([ EneLowLmt-1 EneHiLmt+1 AxisScale(3:4) ]);
xlabel('Emission Frequency (meV)', 'FontSize', 10);
ylabel('Spectral Intensity (a.u.)', 'FontSize', 10);
title('Interferometric Term', 'FontSize', 10);
%set(1,'PaperPositionMode','auto'); % h is figure number
if isSaveImages
    saveas(gcf, strcat(OutDataPath, 'Fig1'), 'pdf');
end

% Inverse Fourier transform to time domain (use FFT)
IFT_InterfSpec = fft(fftshift(InterfSpec_THz, 1), NEmiDim, 1);
IFT_InterfSpec1 = IFT_InterfSpec;
% Cut off residual part near time zero
NEmiInitCut =round(TimeCut*(FreqMax-FreqMin));
IFT_InterfSpec1(1:NEmiInitCut) = 0;
if Second
    if Echo
        NEmiFinCut = NEmiInitCut + d+1;
    else
        NEmiFinCut = round(TimeCut2*(FreqMax-FreqMin));
    end
    IFT_InterfSpec1(NEmiFinCut+1:NEmiDim) = 0;
else
    % Remove negative time delays for proper KK transformation
    IFT_InterfSpec1(NEmiDim/2+1: NEmiDim) = 0;
end
% Plot the time-domain information
fig2 = figure(2);
set(gcf, 'Units', 'inch');
% set(gcf, 'position', [ 8 5.75 4 4 ]);
set(gcf, 'position', [ 6.5 5.75 3 2.75 ]);
plot(RTimeAxis, abs(IFT_InterfSpec), 'b-', RTimeAxis, abs(IFT_InterfSpec1), 'r-');
axis([ 0  RTimeAxis(end)/2  0  1.1*max(abs(IFT_InterfSpec1)) ]);
xlabel('Real Emission Time (ps)', 'FontSize', 10);
ylabel('Intensity (a.u.)', 'FontSize', 10);
title('Fourier transformed interferogram', 'FontSize', 10);
if isSaveImages
    saveas(gcf, strcat(OutDataPath, 'Fig2'), 'pdf');
end
% Fourier transform back to spectral domain (use IFFT)
FT1 = ifft(IFT_InterfSpec1, NEmiDim, 1);
FT1Shift = ifftshift(FT1, 1);
% Calculate retrieved signal
RetrFWM = FT1Shift./ sqrt(RefSpec_THz);
% Find phase
if isGDDcorrection
    PhaseRetrFWM = angle( RetrFWM .* exp(-i * 2*pi* FreqAxis .* Delay_t0) .* exp(-i * 4*pi^2*GDD/2 * (FreqAxis - GDD0).^2) );
else
    PhaseRetrFWM = angle( RetrFWM .* exp(-i * 2*pi* FreqAxis .* Delay_t0) );
end
PhaseRetrFWM1 = unwrap(PhaseRetrFWM(NFirstCol:NLastCol));
RetrFWM1 = RetrFWM(NFirstCol:NLastCol);
FWMSpec_THz1 = FWMSpec_THz(NFirstCol:NLastCol);
FreqAxis1 = FreqAxis(NFirstCol:NLastCol);
RefSpec_THz1 = RefSpec_THz(NFirstCol:NLastCol);
EneAxis1 = FreqAxis1 * Planck;
dlmwrite( strcat(OutDataPath, 'RetrEfield.dat'), [FreqAxis abs(RetrFWM) PhaseRetrFWM], '\t');

% Find the overall phase by matching Pump-probe measurement
if isSRDTrephasing == 0   % no SRDTrephasing;
    OverallPhiFactor = (PhiSI-PhiCam-PhiTR)/pi;      % varies in between 0 & 2.0, in unit of pi
else
    DTSpec1 = DTSpec(NFirstCol:NLastCol);       
    DTSpec1 = DTSpec1 - mean(DTSpec1(1:100)); %Maybe take it from the original DT spec instead?
    DTSpec1 = DTSpec1 ./ max(DTSpec1);  
    SRDT0 = sqrt(RefSpec_THz1) .* abs(RetrFWM1) .* exp(i *PhaseRetrFWM1); %TrSpec
    for IdxM = 0 : 1 : 200
        SRDT1 = real( SRDT0 .* exp(-i* 2*pi * IdxM/200) );
        %SRDT1 = SRDT1 - mean(SRDT1(1:100));
        SRDT1 = SRDT1 ./ max(SRDT1);  
        SRDTDevi(IdxM+1) = sqrt( sum(abs(SRDT1-DTSpec1).^2, 1) /(length(SRDT1)-1) );
    end
    [MinDevi, MinDeviInd] = min(SRDTDevi);
    if isSRDTrephasing == 1
        OverallPhiFactor = 2 * (MinDeviInd-1)/200  % in unit of pi
    else
        OverallPhiFactor = isSRDTrephasing % in unit of pi
    end
    SRDTfn = real(SRDT0 * exp(-i* pi* OverallPhiFactor));    
    %SRDTfn = SRDTfn - mean(SRDTfn(1:100));
    SRDTfn = SRDTfn ./ max(SRDTfn);
    
    fig10 = figure(10);
    set(gcf, 'Units', 'inch');
%     set(gcf, 'position', [ 12 5.75 4 4 ]);
    set(gcf, 'position', [ 9.75 5.75 3 2.75 ]);
    %plot( FreqAxis1, SRDTfn, 'r-', FreqAxis1, DTSpec1, 'k.' );
    %axis([ FreqLowLmt-1 FreqHiLmt+1 min(SRDTfn)-0.2 max(SRDTfn)+0.2 ]);
    plot( EneAxis1, SRDTfn, 'r-', EneAxis1, DTSpec1, 'k.');
        %,EneAxis1,abs(RetrFWM1)./max(abs(RetrFWM1)),'b-');
    %axis([ EneLowLmt-1 EneHiLmt+1 min(SRDTfn)-0.2 max(SRDTfn)+0.2 ]);
    xlabel('Emission Freq (THz)', 'FontSize', 12);
    ylabel('Matched SRDT (normalized)', 'FontSize', 12);
    % legend('Pump-Probe with adjusted phase', 'Experimental SRDT', 'FontSize', 12);
    title('Rephasing with SRDT', 'FontSize', 12);
    text(EneLowLmt-1.5, max(SRDTfn)-0.1, ...
        strcat('Overall Phase = ', num2str(OverallPhiFactor,'%5.2f'), '*pi'), 'FontSize', 10);
    text(EneLowLmt-1.5, max(SRDTfn)-0.15, ...
        strcat('Std Dev. = ', num2str(MinDevi,'%6.4f')), 'FontSize', 10);
    if isSaveImages
        saveas(gcf, strcat(OutDataPath, 'Fig10_mod'), 'pdf');
    end
        dlmwrite(strcat(OutDataPath, 'DTMatched.dat'), [FreqAxis1 DTSpec1 SRDTfn], '\t');
    fig11 = figure(11);    
    set(gcf, 'Units', 'inch');
%     set(gcf, 'position', [ 12 0.75 4 4 ]);  
    set(gcf, 'position', [ 9.75 2 3 2.75 ]);
    plot( (0: 1 : 200)/200 *2, SRDTDevi, '.');
    xlabel('Overall phase factor (pi)', 'FontSize', 12);
    ylabel('Weighted Standard Deviation', 'FontSize', 12);
end

RetrFWMR = abs(RetrFWM) .* real( exp(i*(PhaseRetrFWM-OverallPhiFactor*pi)) );
RetrFWMI = abs(RetrFWM) .* imag( exp(i*(PhaseRetrFWM-OverallPhiFactor*pi)) );
% Save retrieved FWM efield, FreqAxis AMplitude, phase, Real and Imag
dlmwrite( strcat(OutDataPath, 'RetrFWMefieldAll.dat'), ...
    [FreqAxis abs(RetrFWM) PhaseRetrFWM RetrFWMR RetrFWMI], '\t');
% Save comparison, FreqAxis, Retred FWM intensity, phase and Expr FWM intensity
dlmwrite( strcat(OutDataPath, 'CompRetrFWMefield_Expr.dat'),...
    [ FreqAxis1 abs(RetrFWM1.^2)./max(abs(RetrFWM1.^2)) PhaseRetrFWM1 ...
    FWMSpec_THz1./max(FWMSpec_THz1) ],  '\t');

fig3 = figure(3);
set(gcf, 'Units', 'inch');
% set(gcf, 'position', [ 0 0.75 4 4 ]);
set(gcf, 'position', [ 0 2 3 2.75 ]);
subplot(211);
% This is the factor used to normalize the 2D data
MaxSRFWMAmp = max(abs(RetrFWM1));
plot(EneAxis1, FWMSpec_THz1./max(FWMSpec_THz1), 'k-', ...
    EneAxis1, abs(RetrFWM1.^2)./max(abs(RetrFWM1.^2)), 'r-' );
AxisScale = axis;    axis([ EneLowLmt EneHiLmt 0 1.1 ]);
xlabel('Emission Frequency (meV)', 'FontSize', 10);
ylabel('Intensity (a.u.)', 'FontSize', 10);
title('Retrieved(Red) & Measured(Blk) FWM Intensity', 'FontSize', 10);
subplot(212);
IndexLowLmt1 = IndexLowLmt - NFirstCol+1;
IndexHiLmt1 = IndexHiLmt - NFirstCol+1;
plot(EneAxis1, PhaseRetrFWM1, 'r-');
axis([ EneLowLmt EneHiLmt...
    min(PhaseRetrFWM1(IndexLowLmt1:IndexHiLmt1))-pi/8 ...
    max(PhaseRetrFWM1(IndexLowLmt1:IndexHiLmt1))+pi/8]);
xlabel('Emission Frequency (meV)', 'FontSize', 10);
ylabel('Ripped Phase (rad.)', 'FontSize', 10);
title('Retrieved FWM Phase, unwrapped', 'FontSize', 10);
if isSaveImages
    saveas(gcf, strcat(OutDataPath, 'Fig3'), 'pdf');
end
    
%% Process the whole data if confirmed at beginning

if strcmp(HowProcessBtn, 'One') return; end;
mmmmm=0;
while 1

clear AbsFreqAxis % This is necessary because some of the flipping and transposes later on are no longer restricted to mmmmm = 0 
mStr=num2str(mmmmm);

if is2Dor3D == '2D'
    SIfname = 'SI0total1.dat';
elseif is2Dor3D == '3D'
    SIfname = strcat(FName2Ddata,'/SItotal',isPrepulseOrDelayT,mStr,'.dat');
    fid=fopen(strcat(DataPath,SIfname), 'r');
    if (fid==-1), break, end
    fclose(fid);
end
MatIntf1 = dlmread(strcat(DataPath, SIfname), '\t');
MatIntf1 = flipud(MatIntf1); % See beginning for when to comment this out!

% Select how many tau steps to include
if TauSteps
    MatIntf1 = MatIntf1(:,1:TauSteps);
end

NAbsDim0 = size(MatIntf1, 2);

if isResample
    [MatIntf1, FreqMin, FreqMax] = resampleWL2Freq(MatIntf1, WLMin, WLMax, NEmiDim);
else
    MatIntf1 = MatIntf1(end:(-1):1, :);
end

if isPhaseCycling == 0
    MatIntf1 = MatIntf1 - sum([FWMSpec1(1:10); FWMSpec1(1015:end)])/20; % subtract out dark counts (using initial FWM data)
    MatIntf1 = MatIntf1 - RefRatio * repmat(RefSpec_THz, 1, NAbsDim0); % subtract out ref.
end
MatIntf1 = transpose(MatIntf1);
% Inverse Fourier transform of MatIntf2 along EmisFreq Dim (to time-domain)
IFT_MatIntf = fft(fftshift(MatIntf1, 2), NEmiDim, 2);
clear MatIntf1;
    
    % Define 2D time cut window
    TimeWindow = ones(size(IFT_MatIntf));
    TimeWindow2 = TimeWindow;
    NRow = size(IFT_MatIntf,1);
    if Second
        if Echo
            % Follow photon echo
            for j = 1:NRow
                if j>c
                    TimeWindow(j,1:NEmiInitCut+floor(m*(j-c))) = 0;
                    TimeWindow(j,NEmiInitCut+floor(m*(j-c))+d+1:end) = 0;
                else
                    TimeWindow(j, 1 : NEmiInitCut) = 0;
                    TimeWindow(j,NEmiInitCut+d+1:end) = 0;
                end
            end
        else
            TimeWindow(:, 1 : NEmiInitCut) = 0;
            TimeWindow(:, NEmiFinCut : end) = 0;
        end
    else
        TimeWindow(:, 1 : NEmiInitCut) = 0;
        TimeWindow(:, NEmiDim/2+1: NEmiDim) = 0;
    end
    IFT_MatIntf = IFT_MatIntf.*TimeWindow;
    
    % Check time cut window
    for j = 1:NRow
        if j>c
            TimeWindow2(j,1:NEmiInitCut+floor(m*(j-c))) = 0;
            TimeWindow2(j,NEmiInitCut+floor(m*(j-c))+d+1:end) = 0;
        else
            TimeWindow2(j, 1 : NEmiInitCut) = 0;
            TimeWindow2(j,NEmiInitCut+d+1:end) = 0;
        end
    end
    
    if isDeconvolve
        DeconMat = exp(SpexG*RTimeAxis);
        DeconMat(NEmiDim/2+1:NEmiDim) = fliplr(DeconMat(1:NEmiDim/2));
%         figure(11);
%         plot(RTimeAxis,DeconMat);
%         xlabel('Real Time (ps)');
%         ylabel('Multiplication Factor');
        MultMat = repmat(DeconMat,NAbsDim0,1);
        IFT_MatIntf = IFT_MatIntf.*MultMat;
    end

    TauAxis1 = UdrSmplRatio / 2 / HeNeFreq * 1000 * (0:1:NAbsDim0-1)/1000;

    figure(6);
    set(gcf, 'Units', 'inch');
    set(gcf, 'position', [ 9.75 4 5 2.75 ]);
    contour(abs(IFT_MatIntf),100); % Plot 2D time domain
%     contour(RTimeAxis,TauAxis1,abs(IFT_MatIntf),100); % Plot 2D time domain
%     xlim([0 RTimeAxis(NEmiDim/2)]);
%     xlabel('Real Time (ps)');
%     ylabel('\tau (ps)');
%     title('Signal-Reference interference term in 2D time domain');
    hold on;
    contour(TimeWindow2,10);
    hold off;
    
% Perform FFT back to spectral-domain
MatFTShift = ifftshift( ifft(IFT_MatIntf, NEmiDim, 2), 2 );
clear IFT_MatIntf;
MatFT = MatFTShift ./ repmat(transpose(sqrt(RefSpec_THz)), NAbsDim0, 1);
clear MatFTShift;

% 1. remove linear term in the FWM phase,
% PhaseTerm1 is the linear phase exponential term, size: NEmiDimx1.
% 2. add the overall phase to FWM phase, OverallPhiFactor is scalar.
if isGDDcorrection
    PhaseTerm1 = transpose( exp(-i * 2*pi* FreqAxis .* Delay_t0 -i * 4*pi^2*GDD/2 * (FreqAxis - GDD0).^2) );
else
    PhaseTerm1 = transpose( exp(-i * 2*pi* FreqAxis .* Delay_t0) );
end
MatFT = MatFT .* repmat(PhaseTerm1, NAbsDim0, 1) .* exp(-i * OverallPhiFactor*pi);
    
if isOffsetCor == 1
    if TauIndex
        Offset = abs(MatFT(TauIndex,:));
    else
        Offset = abs(MatFT(NRow,:));
    end
    OffsetMat = repmat(Offset,NRow,1);
    AngleMat = angle(MatFT);
    OffsetMat = OffsetMat.*exp(1i*AngleMat);
    MatFT1 = MatFT - OffsetMat;
elseif isOffsetCor == 2 %Generalization of above. Neither of these really make sense.
    if TauIndex
        Offset = mean(mean(abs(MatFT(TauIndex-20:TauIndex,:))));
    else
        Offset = mean(mean(abs(MatFT(NRow-20:NRow,:))));
    end
    OffsetMat = Offset*ones(NRow,NEmiDim);
    AngleMat = angle(MatFT);
    OffsetMat = OffsetMat.*exp(i*AngleMat);
    MatFT1 = MatFT - OffsetMat;
elseif isOffsetCor == 3 %Good idea if the offset has a well-defined phase.
    if TauIndex
        Offset = mean(mean(MatFT(TauIndex-20:TauIndex,:)));
    else
        Offset = mean(mean(MatFT(NRow-20:NRow,:)));
    end
    OffsetMat = Offset*ones(NRow,NEmiDim);
    MatFT1 = MatFT - OffsetMat;      
else % Apply a window function along the tau direction
    % Hanning window
    %WindowFunc = 0.50 + 0.50*cos(pi.*transpose(linspace(0, 1, NAbsDim0)));
    % Hamming window
    %WindowFunc = 0.54 + 0.46*cos(pi.*transpose(linspace(0, 1, NAbsDim0)));
    % Arctan window
    %WindowFunc = atan(15*(1-transpose(linspace(0, 1, NAbsDim0))))./atan(15);
    % Tukey window (Reduces to Hanning when alpha=1)
    alpha = TukeyAlpha;
    alphaNAbsDim0 = round(alpha*NAbsDim0);
    WindowFunc(1:alphaNAbsDim0) = 0.5*(1+cos(pi*(2*(0:alphaNAbsDim0-1)/2/alphaNAbsDim0-1)));
    WindowFunc(alphaNAbsDim0+1:NAbsDim0) = 1;
    WindowFunc = flipud(transpose(WindowFunc));
    % Modified Tukey window, rounding off the initial heaviside.
%         alpha = TukeyAlpha; %Cannot get too big, or the function will fail.
%         TimeBandwidth = 80; %fs
%         TauAxis = UdrSmplRatio / 2 / HeNeFreq * 1000 * (0:1:(length(MatFT)-1));
%         alphares = TimeBandwidth*2/TauAxis(end);
%         alphaNAbsDim0 = round(alpha*NAbsDim0);
%         alpharesNAbsDim0 = round(alphares*NAbsDim0);      
%         WindowFunc(1:alphaNAbsDim0) = 0.5*(1+cos(pi*(2*(0:alphaNAbsDim0-1)/2/alphaNAbsDim0-1)));
%         WindowFunc(alphaNAbsDim0+1:NAbsDim0) = 1;
%         WindowFunc = fliplr(WindowFunc);   
%         WindowFuncc = ones(1,NAbsDim0);
%         WindowFuncc(1:alpharesNAbsDim0) = 0.5*(1+cos(pi*(2*(0:alpharesNAbsDim0-1)/2/alpharesNAbsDim0-1)));
%         WindowFunc = transpose(WindowFunc.*WindowFuncc);
    MatFT1 = MatFT .* repmat(WindowFunc, 1, size(MatFT, 2));      
end
clear WindowFunc; %CLS: This is very important for multiple spectra because otherwise WindowFunc will get transposed more than it should.
    
% Plot the retrieved FWM intensity vs. delay Tau at diff wavelength
fig4 = figure(4);
set(gcf, 'Units', 'inch');
% set(gcf, 'position', [ 4 0.75 4 4 ]);
set(gcf, 'position', [ 3.25 2 3 2.75 ]);
set(fig4, 'Name', 'Retrieved FWM intensity vs. Tau');
MatFTPrjInten = sum( (abs(MatFT)).^2, 2);
MatFT1PrjInten = sum( (abs(MatFT1)).^2, 2);
TauAxis = UdrSmplRatio / 2 / HeNeFreq * 1000 * (0:1:(length(MatFTPrjInten)-1));
% plot(TauAxis, MatFTPrjInten, 'b-', TauAxis, MatFT1PrjInten, 'r-');
semilogy(TauAxis, MatFTPrjInten, 'b-', TauAxis, MatFT1PrjInten, 'r-');
axis([0 max(TauAxis) min(MatFTPrjInten/10) max(MatFTPrjInten)]);
xlabel('Delay Tau (fs)', 'FontSize', 12);
ylabel('FWM Intensity (a.u.)', 'FontSize', 12);
title('Integrated FWM intensity vs. delay Tau', 'FontSize', 12);
% Write out retrieved TI-FWM to data file
% dlmwrite(strcat(fftDataPath, 'RtrFWMvsTau',num2str(fffff),'.dat'), [ TauAxis' MatFTPrjInten ], '\t');
dlmwrite(strcat(OutDataPath, 'RtrFWMvsTau.dat'), [ TauAxis' MatFTPrjInten ], '\t');
if isSaveImages
%     saveas(gcf, strcat(fftDataPath, 'Fig4',num2str(fffff)), 'pdf');
    saveas(gcf, strcat(OutDataPath,'Fig4',isPrepulseOrDelayT,num2str(mmmmm)), 'pdf');
end
clear MatFTPrjInten MatFT1PrjInten;
    
% Pad zeros at the end of matrix to increase length in tau
% F.T. FWM matrix with respect to tau (along the absorption dimension /row)
NAbsDim = (NAbsDim > NAbsDim0)*NAbsDim + (NAbsDim <= NAbsDim0)*NAbsDim0;
MatFTPadded = [ MatFT1; zeros(NAbsDim-size(MatFT1,1), size(MatFT1,2)) ];
M2D = fftshift(fft((MatFTPadded), NAbsDim, 1 ), 1);
clear MatFT MatFT1 MatFTPadded; 
    
% -----------------------------------------------------------------------
% Generate spectrum in the right range from the undersampled data
if (UdrSmplRatio == 1)||(UdrSmplRatio == 2)||(UdrSmplRatio == 4)||(UdrSmplRatio == 8)||(UdrSmplRatio == 16)||(UdrSmplRatio == 32)||(UdrSmplRatio == 64)
    EmiFreqStep = (FreqMax - FreqMin)/NEmiDim;
    AbsFreqStep = 2*NyquistFreq/NAbsDim;
    NSmplDim = round((NLastCol-NFirstCol)*EmiFreqStep/AbsFreqStep);
    %New frequency calculation (2016-05-31) starts here.
    %Non-rephasing version of this works for the 3D data. Hopefully also for S2 nonrephasing. 
    TauDelta = TauAxis(2)-TauAxis(1); %fs.
    FreqDelta = 1000/size(M2D,1)/TauDelta;
    AbsFreqAxisM2D = FreqDelta*(-size(M2D,1)/2:size(M2D,1)/2-1)';
    %Conventions below are a little twisted around because of the weird choice of using fft above instead of ifft, but this works anyway.
    %FreqOffset = HeNeFreq*2/UdrSmplRatio*round(FreqCenter/(HeNeFreq*2/UdrSmplRatio)); %Somehow seems more correct for HeNe steps of 633/2.
    if isRephasing == 1
        FreqOffset = 2*NyquistFreq*round(FreqCenter/(2*NyquistFreq)); %Same as previous line. See Chris Smallwood's (physical) Lab Notebook 2, p. 29.    
    elseif isRephasing == 0
        FreqOffset = -2*NyquistFreq*round(FreqCenter/(2*NyquistFreq));
    end
    AbsFreqAxisM2D = AbsFreqAxisM2D+FreqOffset;     
    % Now cut the spectrum w_tau axis down to relevant window size:
    for k = 1:NSmplDim
        FreqUnfolded = FreqMin + NFirstCol*EmiFreqStep + (k-1)*AbsFreqStep;
        NFolding = round(FreqUnfolded/(2*NyquistFreq));
        FreqFolded = (isRephasing==1)*(FreqUnfolded - 2*NFolding*NyquistFreq)...
            - (isRephasing==0)*(FreqUnfolded - 2*NFolding*NyquistFreq);
        if isRephasing == 1
            NAbsRow = round((FreqFolded + NyquistFreq)/AbsFreqStep); 
        elseif isRephasing==0
            NAbsRow = round((FreqFolded + NyquistFreq)/AbsFreqStep)+2; %This is a poor way to do things TKTK...
        end
        if NAbsRow < 1
            NAbsRow = NAbsRow + NAbsDim;
        elseif NAbsRow > NAbsDim
            NAbsRow = NAbsRow - NAbsDim;
        end
        M2DUdlSmpl(k,:) = M2D(NAbsRow, NFirstCol:NLastCol );
        if k == 1
            AbsFreqAxis(k) = AbsFreqAxisM2D(NAbsRow);
        elseif isRephasing == 1
            AbsFreqAxis(k) = AbsFreqAxis(k-1) + AbsFreqStep;
        else 
            AbsFreqAxis(k) = AbsFreqAxis(k-1) - AbsFreqStep;
        end
    end 
    AbsFreqAxisM2D = -AbsFreqAxisM2D; %Weird conventions get rectified here.
    AbsFreqAxis = -AbsFreqAxis; %Weird conventions get rectified here.
    if isRephasing == 1
        M2DUdlSmpl=flipud(M2DUdlSmpl); %Flipping makes some difference for MDim for S1.
        AbsFreqAxis =flipud(transpose(AbsFreqAxis));
    elseif isRephasing == 0
        AbsFreqAxis =transpose(AbsFreqAxis);
    end
    AbsEneAxis = Planck*AbsFreqAxis; 
else
    disp('Unexpected under-sampling ratio! Only 1, 2, 4 or 8 is valid.');
    return;
end

% ****************OPTIONAL********************************
clear M2D;    
fig5 = figure(5);
set(gcf, 'Units', 'inch');
% set(gcf, 'position', [ 8 0.75 4 4 ]);
set(gcf, 'position', [ 6.5 2 3 2.75 ]);
%plot(AbsEneAxis, sum(abs(M2DUdlSmpl).^2, 2), 'b-');
plot(AbsEneAxis, sum(abs(M2DUdlSmpl), 2), 'b-'); % CLS: Previous version seems strange. Already interacting with two fields, so it's already quadratic, isn't it?...
set(fig5, 'Name', 'Projection of 2D intensity to absorption frequency axis');
xlabel('Absorption Frequency (meV)', 'FontSize', 12);
ylabel('projected 2D intensity', 'FontSize', 12);
title( 'Projection of 2D intensity to absorption frequency axis', 'FontSize', 12);
if isSaveImages
    saveas(gcf, strcat(OutDataPath, 'Fig5',num2str(mmmmm)), 'pdf');
end
%dlmwrite(strcat(fftDataPath, 'MProjAbsRange',num2str(fffff),'.dat'), [AbsFreqAxis, sum(abs(M2DUdlSmpl).^2, 2)], '\t');
%dlmwrite(strcat(fftDataPath, 'MProjAbsRange.dat'), [AbsFreqAxis, sum(abs(M2DUdlSmpl).^2, 2)], '\t');
dlmwrite(strcat(OutDataPath, 'MProjAbsRange.dat'), [AbsFreqAxis, sum(abs(M2DUdlSmpl), 2)], '\t'); % CLS: See comment 10 lines up.
% ****************END OPTION********************************

% Code that existed for standard FTSA, not 3D
%     if mmmmm > 1
%         M2DUdlSmpl = flipud(M2DUdlSmpl);
%     end
     VMax = max(max(abs(M2DUdlSmpl)));

%% Normalize to the maximum value. 

% MaxAmp = max(max(abs(M2DUdlSmpl)))/2; %This is different from standard (2D version before merging) FTSA
% MaxAmp = max(max(abs(M2DUdlSmpl)));
MaxAmp = 1;

MAmpl = abs(M2DUdlSmpl) ./ MaxAmp;
MReal = real(M2DUdlSmpl)./ MaxAmp;
MImag = imag(M2DUdlSmpl)./ MaxAmp;
% Standard FTSA writes these to file at this point, but this version waits until the end.

%% Prepare the scales of 2D plot in energy (meV)
MDim(1:2) = [FreqAxis(NFirstCol), FreqAxis(NLastCol)];
MDim(3:4) = [AbsFreqAxis(1), AbsFreqAxis(end)];

% Export the frequency ranges of 2D matrix to a file, in the order of
% EmiFreq_min, EmiFreq_max, AbsFreq_min, and AbsFreq_max.
if mmmmm == 0
    fid = fopen( strcat(OutDataPath,'MFreqRange.txt'), 'w');
    fprintf(fid, 'NFirstCol \t NLastCol \t NSmplDim \r\n');
    fprintf(fid, '%10.0f \t %9.0f \t %10.0f \r\n', NFirstCol,NLastCol,NSmplDim);
    fprintf(fid, 'EmiFreq_min \t EmiFreq_max \t AbsFreq_min \t AbsFreq_max \r\n');
    fprintf(fid, '%11.4f \t %11.4f \t %11.4f \t %11.4f \r\n\r\n', MDim);
    MDimEngy = Planck.*MDim;
    fprintf(fid, 'EmiEngy_min \t EmiEngy_max \t AbsEngy_min \t AbsEngy_max \r\n');
    fprintf(fid, '%11.4f \t %11.4f \t %11.4f \t %11.4f \r\n\r\n', MDimEngy);
    % Output important parameters to the file
    fprintf(fid, 'TimeCut = %5.3f \r\n', TimeCut);
    fprintf(fid, 'OverallPhiFactor = %5.3f \r\n', OverallPhiFactor);
    fprintf(fid, 'PhiSI = %5.3f \r\n', PhiSI);
    fprintf(fid, 'PhiCam = %5.3f \r\n', PhiCam);
    fprintf(fid, 'PhiTR = %5.3f \r\n', PhiTR);
    fprintf(fid, 'Delay_t0 = %5.3f \r\n', Delay_t0);
    fprintf(fid, 'DataPath = %s \r\n', DataPath);
    fprintf(fid, 'Rephasing? = %g \r\n', isRephasing);
    fprintf(fid, 'Peak of SR-FWM ampl at tau=0: %7.4g \r\n', MaxSRFWMAmp);
    fprintf(fid, 'Peak of 2D amplitude matrix: %7.4g \r\n', MaxAmp);
    fprintf(fid, 'Normalization factor: %7.4g \r\n', MaxAmp/MaxSRFWMAmp);
    fclose(fid);
end 

% Plot Amp, real and imaginary parts of 2D data
[gEmiFreq, gAbsFreq] = meshgrid(linspace( MDimEngy(1), MDimEngy(2), NLastCol-NFirstCol+1 ),...
                                linspace( MDimEngy(3), MDimEngy(4), NSmplDim));
if isRephasing == 1   % Rephasing scan
    %FigName = strcat('Rephasing data of: ', DataPath);
    AxisRange = [ EneLowLmt EneHiLmt -EneHiLmt -EneLowLmt ];
    DiagnlLine = [linspace(EneLowLmt, EneHiLmt, 20); linspace(-EneLowLmt, -EneHiLmt, 20)];
elseif isRephasing == 0   % Non-Rephasing scan
    %FigName = strcat('Non-Rephasing data of: ', DataPath);
    AxisRange = [ EneLowLmt EneHiLmt EneLowLmt EneHiLmt ];
    DiagnlLine = [linspace(EneLowLmt, EneHiLmt, 20); linspace(EneLowLmt, EneHiLmt, 20)];
end
%dlmwrite(strcat(fftDataPath,'gEmiFreq.dat'),gEmiFreq,'\t'); %These get written in the original 2D version but not 3D version. TKTK.
%dlmwrite(strcat(fftDataPath,'gAbsFreq.dat'),gAbsFreq,'\t');

AxisRange1 = [ EneLowLmt EneHiLmt 0 1.1 ];
    
%% ************************ OPTION ************************** 
% Plot the 2D matrices (Optional)

    if Real
        fig7 = figure(7);
        set(gcf, 'Units', 'inch');
        set(gcf, 'position', [ 1 1 12 5 ]);
        subplot(231);
        plot(EneAxis, RefSpec_THz/max(RefSpec_THz), 'b-');
        %     plot(EneAxis, RefSpec_THz/max(RefSpec_THz)*max(linabs), 'b-', linenergy, linabs, 'k-');
        axis( AxisRange1 );
        set(gca, 'Units', 'inch');
        set(gca, 'position', [ 0.75 3.5 3 1 ]);
        subplot(232);
        plot(EneAxis, RefSpec_THz/max(RefSpec_THz), 'b-');
        %     plot(EneAxis, RefSpec_THz/max(RefSpec_THz)*max(linabs), 'b-', linenergy, linabs, 'k-');
        axis( AxisRange1 );
        set(gca, 'Units', 'inch');
        set(gca, 'position', [ 4.5 3.5 3 1 ]);
        subplot(233);
        plot(EneAxis, RefSpec_THz/max(RefSpec_THz), 'b-');
        %     plot(EneAxis, RefSpec_THz/max(RefSpec_THz)*max(linabs), 'b-', linenergy, linabs, 'k-');
        axis( AxisRange1 );
        set(gca, 'Units', 'inch');
        set(gca, 'position', [ 8.25 3.5 3 1 ]);
        subplot(234);
%         VMax = 1;
        hFig = contour(gEmiFreq, gAbsFreq, MAmpl, linspace(0, VMax, 1*NContourLevels), 'LineWidth', 1.5);
        axis( AxisRange );
        line(DiagnlLine(1,:), DiagnlLine(2,:), 'LineStyle', ':', 'Color', [0 0 0]);
        axis square;  colormap jet;
        set(gca, 'Units', 'inch');
        set(gca, 'position', [ 0.75 0.5 3 3 ]);
        colorbar('Location','East');
        xlabel('Emission Energy (meV)');
        ylabel('Excitation energy (meV)');
        title('Absolute value');
        subplot(235);
        MReal(end, end-1)= -1;  MReal(end, end)= 1;
        hFig = contour(gEmiFreq, gAbsFreq, MReal, linspace(-VMax, VMax, NContourLevels), 'LineWidth', 1.5);
        axis( AxisRange );
        caxis([-VMax,VMax]);
        line(DiagnlLine(1,:), DiagnlLine(2,:), 'LineStyle', ':', 'Color', [0 0 0]);
        axis square;  colormap jet;
        set(gca, 'Units', 'inch');
        set(gca, 'position', [ 4.5 0.5 3 3 ]);
        colorbar('Location','East');
        xlabel('Emission Energy (meV)');
        ylabel('Excitation energy (meV)');
        title('Real part');
        subplot(236);
        MImag(end, end-1)= -1;  MImag(end, end)= 1;
        hFig = contour(gEmiFreq, gAbsFreq, MImag, linspace(-VMax, VMax, NContourLevels), 'LineWidth', 1.5);
        axis( AxisRange );
        caxis([-VMax,VMax]);
        line(DiagnlLine(1,:), DiagnlLine(2,:), 'LineStyle', ':', 'Color', [0 0 0]);
        axis square;  colormap jet;
        set(gca, 'Units', 'inch');
        set(gca, 'position', [ 8.25 0.5 3 3 ]);
        colorbar('Location','East');
        xlabel('Emission Energy (meV)');
        ylabel('Excitation energy (meV)');
        title('Imaginary part');
        if isSaveImages
            if isRephasing
                saveas(gcf, strcat(OutDataPath,'2DS1',isPrepulseOrDelayT,'Fig7',num2str(mmmmm)), 'png');
            else
                saveas(gcf, strcat(OutDataPath,'2DS2',isPrepulseOrDelayT,'Fig7',num2str(mmmmm)), 'png');
            end
        end
            
        % Save 2D plot in 600dpi png format
        % print(gcf, '-r600', '-dpng', '-noui', strcat(fftDataPath, 'Fig7', '.png'));
    end

if Abs
    %Plot the absolute value matrix (This version has no Reference spectrum)
    % fig8 = figure(8);
    % set(gcf, 'Units', 'inch');
    % set(gcf, 'position', [ 2 2.2 7 7.5 ]);
    % newEmi = gEmiFreq(1,:)';
    % newAbs = gAbsFreq(:,1);
    % hFig = imagesc(newEmi, newAbs, MAmpl); set(gca,'Ydir','Normal');
    % line(DiagnlLine(1,:), DiagnlLine(2,:), 'LineStyle', ':', 'Color', [0 0 0]);
    % %colormap jet;
    % xlabel('Emission Energy (meV)');
    % ylabel('Excitation energy (meV)');
    % if isSaveImages
    %     if isRephasing
    %         saveas(gcf, strcat(DataPath, '2dfftmatrixS1/', '2DS1',isPrepulseOrDelayT, num2str(mmmmm)), 'png');
    %     else 
    %         saveas(gcf, strcat(DataPath, '2dfftmatrixS2/', '2DS2',isPrepulseOrDelayT, num2str(mmmmm)), 'png');
    %     end
    % end  
    
    %Plot the absolute value matrix (with Reference, possibly better alternate)
    fig8 = figure(8);
    set(gcf, 'Units', 'inch');
%     set(gcf, 'position', [ 2 2.2 7 7.5 ]);
    set(gcf, 'position', [ 2 2.25 5 5.5 ]);
    subplot(211);
    plot(EneAxis, RefSpec_THz/max(RefSpec_THz),'-k');
    axis( AxisRange1 );
    set(gca, 'Units', 'inch');
%     set(gca, 'position', [ 1 5.5 5 1.5 ]);
    set(gca, 'position', [ 1 4 3.5 1 ]);
    subplot(212);
        % VMax = 1;
        %      hFig = contour(gEmiFreq, gAbsFreq, MAmpl, linspace(0, VMax, NContourLevels), 'LineWidth', 1.5);
        %hFig = contourf(gEmiFreq, gAbsFreq, MAmpl, linspace(0, VMax, NContourLevels), 'LineStyle', 'none');
        %   hFig = contour(gEmiFreq, gAbsFreq, MReal, linspace(-0.75*VMax, 0.75*VMax, NContourLevels), 'LineWidth', 1.5);
    newEmi = gEmiFreq(1,:)';
    newAbs = gAbsFreq(:,1);
    hFig = imagesc(newEmi, newAbs, MAmpl); set(gca,'Ydir','Normal');
    %hFig = surf(newEmi,newAbs,MAmpl); view(2); set(hFig,'edgecolor','none');        
    axis( AxisRange );
    line(DiagnlLine(1,:), DiagnlLine(2,:), 'LineStyle', ':', 'Color', [0 0 0]);
    colormap jet;
    set(gca, 'Units', 'inch');
%     set(gca, 'position', [ 1 0.5 5 5 ]);
    set(gca, 'position', [ 1 0.5 3.5 3.5 ]);
    xlabel('Emission Energy (meV)');
    ylabel('Excitation energy (meV)');
    %set(8,'PaperPositionMode','auto'); % h is figure number
    if isSaveImages
        if isRephasing    
            saveas(gcf, strcat(OutDataPath,'2DS1',isPrepulseOrDelayT, num2str(mmmmm)),'png');
        else 
            saveas(gcf, strcat(OutDataPath,'2DS2',isPrepulseOrDelayT, num2str(mmmmm)),'png');
        end
    end
    
end
% % ***********************End OPTION *************************

%% Save the 2D matrices after FFT
dlmwrite(strcat(OutDataPath, 'MReal',isPrepulseOrDelayT, mStr,'.dat'), MReal, 'delimiter', '\t', 'precision', 6);
dlmwrite(strcat(OutDataPath, 'MImag',isPrepulseOrDelayT, mStr,'.dat'), MImag, 'delimiter', '\t', 'precision', 6);
dlmwrite(strcat(OutDataPath, 'MAmpl',isPrepulseOrDelayT, mStr,'.dat'), MAmpl, 'delimiter', '\t', 'precision', 6);
msgbox(strcat('The 2D FFT for ',isPrepulseOrDelayT,' ', mStr,' has been generated'),'replace');
mmmmm=mmmmm+1;

if strcmp(HowProcessBtn, 'One 2D spectrum'), break, end;
if is2Dor3D == '2D', break, end;



end

toc;
