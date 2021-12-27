%% 2021-11-19, PA1,
%% 2021-11-22, PA2, input: signal with band information and signle fLO!!!

clear all
close all

fnum = 111901

%% generate waveform based on RXIBW 491.52MHz
%% input: signal
config.bw_Channel = '20MHz'
config.fs = 491.52e6*1; % flag_fs_extend
config.MOD = '64QAM';
config.ratio_SamplesDMC = 'on'; % on:reduce length, off: no reduce
config.Carrier_Type = 'LTE';

%% input: channel filter
fir.Wtype = "kaiser"
fir.Ftype = "LPF"
fir.Order = NaN
fir.fTolerance = -0.1e6
fir.K_AttdB = 60
fir.K_fdelta = 0.5e6
fir.fcutoffL = 20e6/2 % based on config.bw_Channel = '20MHz'
fir.fcutoffH = 0
fir.fs = config.fs

%% input: signal output power
SNRdB_DL2UL = 60
LeakagedB_DL2UL = -60
PodB_UL = -15
PodB_DL = PodB_UL + SNRdB_DL2UL + LeakagedB_DL2UL % SNRdB_DL2UL = 60, Leakage_DL2UL = -60

%% input: signal with band information and signle fLO!!!
disp('Setting the signle LO frequency of IF architectue for B1 and B3 FDD mode, the LO freq. could be tunning for freq. plan')
fLO = 1900e6; %% 2021-11-22, PA2, input: signal with band information and signle fLO!!!

band1DL.band = 'B1'
band1DL.bw = 60e6;
band1DL.fLO = fLO;
band1DL.RFbwInband = [2110 2170]*1e6
band1DL.fIF = mean(band1DL.RFbwInband) - band1DL.fLO

band1UL.band = 'B1'
band1UL.bw = 60e6;
band1UL.fLO = fLO;
band1UL.RFbwInband = [1920 1980]*1e6
band1UL.fIF = mean(band1UL.RFbwInband) - band1UL.fLO

band3DL.band = 'B3'
band3DL.bw = 80e6;
band3DL.fLO = fLO;
band3DL.RFbwInband = [1800 1880]*1e6
band3DL.fIF = mean(band3DL.RFbwInband) - band3DL.fLO

band3UL.band = 'B3'
band3UL.bw = 80e6;
band3UL.fLO = fLO;
band3UL.RFbwInband = [1710 1790]*1e6
band3UL.fIF = mean(band3UL.RFbwInband) - band3UL.fLO

%% input: PA performance
PA.flag = 'on'

%% output: signal
[signal_DL1, Config_DL1] = OFDM_SG_SA(config, fir, PodB_DL, band1DL, PA, ['DL1'], fnum);
[signal_UL1, Config_UL1] = OFDM_SG_SA(config, fir, PodB_UL, band1UL, [], ['UL1'], fnum);

config.fs = Config_UL1.fs; % update the fs, make sure the fs of dualband are the same
[signal_DL3, Config_DL3] = OFDM_SG_SA(config, fir, PodB_DL, band3DL, PA, ['DL3'], fnum);
[signal_UL3, Config_UL3] = OFDM_SG_SA(config, fir, PodB_UL, band3UL, [], ['UL3'], fnum);
fnum = fnum+1

if 0 % save
    Config_DL3.signal_DL3 = signal_DL3;
    Config_UL3.signal_UL3 = signal_UL3;
    Config_DL1.signal_DL1 = signal_DL1;
    Config_UL1.signal_UL3 = signal_UL1;
    save('Waveform\waveform_OFDM_B3DL_7864p3MHz_struct.mat','Config_DL3')
    save('Waveform\waveform_OFDM_B3UL_7864p3MHz_struct.mat','Config_UL3')
    save('Waveform\waveform_OFDM_B1DL_7864p3MHz_struct.mat','Config_DL1')
    save('Waveform\waveform_OFDM_B1UL_7864p3MHz_struct.mat','Config_UL1')
end

fsOut = Config_UL1.fs
NsampsOut = length(signal_UL1)
df = fsOut/NsampsOut

%% Case1: Evaluate the DL leakage impact UL performance
%% signal combination
signal_combToUL1 = signal_DL3 + signal_DL1 + signal_UL1;
signal_combToUL3 = signal_DL3 + signal_DL1 + signal_UL3;

%% BPF for UL signal
%% input: BPF1
BPF1_fir_Wtype = "kaiser"
BPF1_fir_Ftype = "BPF"
BPF1_fir_Order = NaN
BPF1_fir_fTolerance = 10e6
BPF1_fir_K_AttdB = 40
BPF1_fir_K_fdelta = 10e6
BPF1_fir_fcutoffL = band1UL.RFbwInband(1)
BPF1_fir_fcutoffH = band1UL.RFbwInband(2)
BPF1_fir_Export = fsOut
NCarriers = 1

% input: BPF3
BPF3_fir_fcutoffL = band3UL.RFbwInband(1)
BPF3_fir_fcutoffH = band3UL.RFbwInband(2)

%% output: BPF
b_bpf1 = SYM_FIRApp(BPF1_fir_Wtype,BPF1_fir_Ftype,BPF1_fir_Order,BPF1_fir_K_AttdB,BPF1_fir_K_fdelta...
    ,BPF1_fir_fTolerance,BPF1_fir_fcutoffL,BPF1_fir_fcutoffH,df,band1UL.RFbwInband,NCarriers,BPF1_fir_Export,[]+2);
b_bpf1 = b_bpf1{:};

b_bpf3 = SYM_FIRApp(BPF1_fir_Wtype,BPF1_fir_Ftype,BPF1_fir_Order,BPF1_fir_K_AttdB,BPF1_fir_K_fdelta...
    ,BPF1_fir_fTolerance,BPF3_fir_fcutoffL,BPF3_fir_fcutoffH,df,band3UL.RFbwInband,NCarriers,BPF1_fir_Export,[]+2);
b_bpf3 = b_bpf3{:};

%% output: signal + BPF
signal_combToUL1_BPF = conv(signal_combToUL1, b_bpf1, 'same');
PLOT_FFT_dB_g(signal_combToUL1, fsOut, NsampsOut, ['signal'], 'df', 'full', 'pwr', [fnum]);
PLOT_FFT_dB_g(signal_combToUL1_BPF, fsOut, NsampsOut, ['signal + BPF'], 'df', 'full', 'pwr', [fnum]);
title('UL1 signal')

signal_combToUL3_BPF = conv(signal_combToUL3, b_bpf3, 'same');
PLOT_FFT_dB_g(signal_combToUL3, fsOut, NsampsOut, ['signal'], 'df', 'full', 'pwr', [fnum*10]);
PLOT_FFT_dB_g(signal_combToUL3_BPF, fsOut, NsampsOut, ['signal + BPF'], 'df', 'full', 'pwr', [fnum*10]);
title('UL3 signal')
fnum = fnum+1

%% output: UL demodulation
[~, ~, Demod_UL1] = OFDM_SG_SA(Config_UL1, fir, [], band1UL, [], ['UL1'], fnum, 'Demod', signal_combToUL1_BPF);
[~, ~, Demod_UL3] = OFDM_SG_SA(Config_UL3, fir, [], band3UL, [], ['UL3'], fnum*10, 'Demod', signal_combToUL3_BPF);

%% Case2: Evaluate the DL signal to ORX, ORX is RF sampling architecture
fnum = 112501

%% Load or generate waveform for ORX
if 0
    load('Waveform\waveform_OFDM_B3DL_7864p3MHz_struct.mat')
    load('Waveform\waveform_OFDM_B1DL_7864p3MHz_struct.mat')
    signal_DL1 = Config_DL1.signal_DL1;
    signal_DL3 = Config_DL3.signal_DL3;
 
else
    %% ORX is RF sampling architecture, and ORX input signal should incloud the 3*Fundmental freq.
    %% input: ORX ADC
    disp('Setting the ORX RF sampling rate of RF sampling architectue, the fsADC could be tunning for freq. plan ')
    fs_ADC = 5898.24e6
    fFundmax = max(band1DL.RFbwInband)
    NHD = 3+1
    fs_PAout = ceil(fFundmax*NHD/(fs_ADC/2))*fs_ADC
    
    %% output: signal
    config.fs = fs_PAout;
    [signal_DL1, Config_DL1] = OFDM_SG_SA(config, fir, PodB_DL, band1DL, PA, ['DL1'], fnum);
    [signal_DL3, Config_DL3] = OFDM_SG_SA(config, fir, PodB_DL, band3DL, PA, ['DL3'], fnum);
    fnum = fnum+1
end
fs = Config_DL3.fs;
Nsamps = numel(signal_DL3);

%% ORX path leakage/rejection/isolation...
LeakagedB_B3ToB1 = -30
Leakage_B3ToB1 = 10^(LeakagedB_B3ToB1/20);

%% signal DL to ORX
signal_DL1ToORX1 = signal_DL1 + signal_DL3*Leakage_B3ToB1;
signal_DL3ToORX3 = signal_DL1*Leakage_B3ToB1 + signal_DL3;
PLOT_FFT_dB_g(signal_DL1ToORX1, fs, Nsamps, ['signal input to ORX1'], 'df', 'full', 'pwr', [fnum]);
PLOT_FFT_dB_g(signal_DL3ToORX3, fs, Nsamps, ['signal input to ORX3'], 'df', 'full', 'pwr', [fnum]);
fnum = fnum+1

%% ORX: Oberservation Receiver
%% DSP, from 5898.24MHz to 491.52MHz
%% DSP, for B1: ADC1 + NCO1(-2123.732MHz) + DEC3(/3) + DEC2(/2) w/ FIR1 + HB1(/2) + NCO2(-16.268MHz)
%% DSP, for B3: ADC1 + NCO1(-1826.232MHz) + DEC3(/3) + DEC2(/2) w/ FIR1 + HB1(/2) + NCO2(-16.268MHz)

%% ORX path wo BPF: B1 and B3 share the same path
%% input: BPF
b_BPF = {1} 
if (iscell(b_BPF)&&cell2mat(b_BPF)==1) || (isnumeric(b_BPF)&&b_BPF==1)
    disp('No BPF in front of ORX ADC')
    disp('ORX ADC sharing for Dualband B1 and B3')
end

%% ORX, RF Sampling Receiver
%% ORX, ADC Sampling
%% input: ADC Quantization Error
Vref = 1.4
nbits = 10
LSB = Vref/2^nbits
flag_RemoveDC = 1
flag_Dither = 0

%% output: ADC Quantization
if 0
    signal_DL1ToORX1_ADCQuanErr = ADC_Quantizer_Dither(signal_DL1ToORX1, LSB, [], flag_RemoveDC, fs, fnum);
elseif flag_Dither
    [signal_DL1ToORX1_ADCQuanErrDH, disp_legend_1] = ADC_Quantizer_Dither(signal_DL1ToORX1, [Vref, nbits], [3], flag_RemoveDC, fs, fnum);
    [EVM_ADCQunaErrDH, ~] = dsp_evm_timexcorr_inband_g(signal_DL1ToORX1, signal_DL1ToORX1_ADCQuanErrDH, fs, [], [], [])
    disp_ADC = [disp_legend_1, ', EVM:',num2str(round(EVM_ADCQunaErr,2))];
elseif ~flag_Dither
    [signal_DL1ToORX1_ADCQuanErr, disp_legend_1] = ADC_Quantizer_Dither(signal_DL1ToORX1, [Vref, nbits], [], flag_RemoveDC, fs, fnum);
    [EVM_ADCQunaErr, ~] = dsp_evm_timexcorr_inband_g(signal_DL1ToORX1, signal_DL1ToORX1_ADCQuanErr, fs, [], [], [])
    disp_ADC = [disp_legend_1, ', EVM:',num2str(round(EVM_ADCQunaErr,2))];
end
PLOT_FFT_dB_g(signal_DL1ToORX1, fs, Nsamps, ['signal before ADC'], 'df', 'full', 'pwr', [fnum]);
fnum = fnum+1

%% input: ADC noise psd 
ratio_ADC1 = fs_ADC/fs
IpwrdB1Hz_Noise_ADC1 = -180

%% output: ADC sampling
if ~flag_Dither
    [signal_DL1ToORX1_ADCQuanErr_sampling, fs_ORX1_ADC, ~, ~] = FP_RFSampling_UpDownSampling_wiFIR_g(signal_DL1ToORX1_ADCQuanErr, fs, ratio_ADC1, b_BPF, IpwrdB1Hz_Noise_ADC1, [], [], fnum, {['ORX1 ',disp_ADC]}, []);
else
    [signal_DL1ToORX1_ADCQuanErr_sampling, fs_ORX1_ADC, ~, ~] = FP_RFSampling_UpDownSampling_wiFIR_g(signal_DL1ToORX1_ADCQuanErrDH, fs, ratio_ADC1, b_BPF, IpwrdB1Hz_Noise_ADC1, [], [], fnum, {['ORX1 ',disp_ADC]}, []);
end
fnum = fnum+1

%% input: ADC Interleaving Error
flag_ADCError = 'NA' % No Error or Not interleaving ADC
flag_ADCError = 'Level'
flag_ADCError = 'Gain'
flag_ADCError = 'Phase'
flag_ADCError = 'All'

error_ADCItl_Level = 1*0.1*[1, 3, -2, 0]
error_ADCItl_GaindB = 1*0.1*[-1, 0.8, -0.4, 0.5]
error_ADCItl_GaindB = 1*[0.5, 0, -0.5, 0]
error_ADCItl_PhaseDEG = 1*0.1*[3, 1, -1, -2]
    
%% output: ADC Interleaving Error
IpwrdBLimit = -80;
IpwrdBLimit = [];

flag_ADCError = 'NA' % No Error or Not interleaving ADC
fnum2 = [fnum, 5, 1, 1]
[signal_DL1ToORX1_ADCQuanErr_sampling_InterleavingErr,errorCell,freqHDkADC1ErCell] = ADC_Interleaving_Error_g(signal_DL1ToORX1_ADCQuanErr_sampling, fs_ORX1_ADC, error_ADCItl_Level, error_ADCItl_GaindB, error_ADCItl_PhaseDEG, flag_ADCError, IpwrdBLimit, fnum2, [], disp_legend_1);

flag_ADCError = 'Level'
fnum2 = [fnum, 5, 1, 2]
disp_legend_1 = ['Level Error: ',num2str(error_ADCItl_Level)]
[signal_DL1ToORX1_ADCQuanErr_sampling_InterleavingErr,errorCell,freqHDkADC1ErCell] = ADC_Interleaving_Error_g(signal_DL1ToORX1_ADCQuanErr_sampling, fs_ORX1_ADC, error_ADCItl_Level, error_ADCItl_GaindB, error_ADCItl_PhaseDEG, flag_ADCError, IpwrdBLimit, fnum2, [], disp_legend_1);

flag_ADCError = 'Gain'
fnum2 = [fnum, 5, 1, 3]
disp_legend_1 = ['GaindB Error: ',num2str(error_ADCItl_GaindB)]
[signal_DL1ToORX1_ADCQuanErr_sampling_InterleavingErr,errorCell,freqHDkADC1ErCell] = ADC_Interleaving_Error_g(signal_DL1ToORX1_ADCQuanErr_sampling, fs_ORX1_ADC, error_ADCItl_Level, error_ADCItl_GaindB, error_ADCItl_PhaseDEG, flag_ADCError, IpwrdBLimit, fnum2, [], disp_legend_1);

flag_ADCError = 'Phase'
fnum2 = [fnum, 5, 1, 4]
disp_legend_1 = ['PhaseDeg Error: ',num2str(error_ADCItl_PhaseDEG)]
[signal_DL1ToORX1_ADCQuanErr_sampling_InterleavingErr,errorCell,freqHDkADC1ErCell] = ADC_Interleaving_Error_g(signal_DL1ToORX1_ADCQuanErr_sampling, fs_ORX1_ADC, error_ADCItl_Level, error_ADCItl_GaindB, error_ADCItl_PhaseDEG, flag_ADCError, IpwrdBLimit, fnum2, [], disp_legend_1);

flag_ADCError = 'All'
fnum2 = [fnum, 5, 1, 5]
disp_legend_1 = ['Interleaving Level+Gain+Phase Error']
[signal_DL1ToORX1_ADCQuanErr_sampling_InterleavingErr,errorCell,freqHDkADC1ErCell] = ADC_Interleaving_Error_g(signal_DL1ToORX1_ADCQuanErr_sampling, fs_ORX1_ADC, error_ADCItl_Level, error_ADCItl_GaindB, error_ADCItl_PhaseDEG, flag_ADCError, IpwrdBLimit, fnum2, [], disp_legend_1);

%% export: signal_ORX1_ADC
signal_ORX1_ADC = signal_DL1ToORX1_ADCQuanErr_sampling_InterleavingErr;
fs_ORX1_ADC = fs_ADC
Nsamps_ORX1_ADC = numel(signal_ORX1_ADC)
fnum = fnum+1

%% ORX, NCO1+Decimation+NCO2
%% input: NCO1 of B1 and B3
fORX1_NCO1 = -2123.732e6
fORX3_NCO1 = -1826.232e6
t = (0:Nsamps_ORX1_ADC-1).'/fs_ORX1_ADC;
disp('Why fORX1_NCO1 = -2123.732MHz? fORX3_NCO1 = -1826.232MHz? to match for B1 and B3 with NCO2=-16.268MHz, there is Only One NCO2 freqs settings')

%% output: signal + NCO1
signal_ORX1_NCO1 =  signal_ORX1_ADC.*exp(1i*2*pi*fORX1_NCO1*t);
PLOT_FFT_dB_g(signal_ORX1_NCO1, fs_ORX1_ADC, Nsamps_ORX1_ADC, ['signal ORX1 NCO1'], 'df', 'full', 'pwr', [fnum, 4, 1, 1]);

%% input: Decimation FIR
b_DEC3 = [0.000976563, 0.000976563, -0.001953125, -0.012207031, -0.017578125, -0.006347656, 0.034179688, 0.065917969, 0.045898438, -0.063476563, -0.181152344, -0.176757813, 0.088867188, 0.55859375, 1.03515625, 1.2265625, 1.03515625, 0.55859375, 0.088867188, -0.176757813, -0.181152344, -0.063476563, 0.045898438, 0.065917969, 0.034179688, -0.006347656, -0.017578125, -0.012207031, -0.001953125, 0.000976563, 0.000976563];
b_DEC3_norm = b_DEC3/b_DEC3(ceil(numel(b_DEC3)/2))
PLOT_FFT_dB_g(b_DEC3_norm*Nsamps_ORX1_ADC, fs_ORX1_ADC, Nsamps_ORX1_ADC, ['FIRDEC3'], 'df', 'full', 'pwr', [fnum, 4, 1, 1]);

%% output: FIR+Decimation(ratio_DEC3 = 3)
signal_ORX1_FIR = DSP_filter_g(b_DEC3_norm, signal_ORX1_NCO1, 'FD');
PLOT_FFT_dB_g(signal_ORX1_FIR, fs_ORX1_ADC, Nsamps_ORX1_ADC, ['signal ORX1 + FIRDEC3'], 'df', 'full', 'pwr', [fnum, 4, 1, 1]);

ratio_DEC3 = 3
signal_ORX1_FIR_DEC3 = signal_ORX1_FIR(1:ratio_DEC3:end);
fs_ORX1_DEC3 = fs_ORX1_ADC/ratio_DEC3
Nsamps_ORX1_DEC3 = numel(signal_ORX1_FIR_DEC3)
PLOT_FFT_dB_g(signal_ORX1_FIR_DEC3, fs_ORX1_DEC3, Nsamps_ORX1_DEC3, ['signal ORX1 + FIR DEC3'], 'df', 'full', 'pwr', [fnum, 4, 1, 2]);

%% input: Decimation FIR1
b_FIR1 = [-0.000976563, 0, 0.005859375, 0, -0.021484375, 0, 0.0625, 0, -0.165039063, 0, 0.619140625, 1, 0.619140625, 0, -0.165039063, 0, 0.0625, 0, -0.021484375, 0, 0.005859375, 0, -0.000976563];
b_FIR1_norm = b_FIR1/b_FIR1(ceil(numel(b_FIR1)/2))
PLOT_FFT_dB_g(b_FIR1_norm*Nsamps_ORX1_DEC3, fs_ORX1_DEC3, Nsamps_ORX1_DEC3, ['FIR1'], 'df', 'full', 'pwr', [fnum, 4, 1, 2]);

%% output: FIR1+Decimation(ratio_FIR1 = 2)
signal_ORX1_FIR1 = DSP_filter_g(b_FIR1_norm, signal_ORX1_FIR_DEC3, 'FD');
PLOT_FFT_dB_g(signal_ORX1_FIR1, fs_ORX1_DEC3, Nsamps_ORX1_DEC3, ['signal ORX1 + FIR DEC3 + FIR1'], 'df', 'full', 'pwr', [fnum, 4, 1, 2]);

ratio_FIR1 = 2
signal_ORX1_FIR1DEC = signal_ORX1_FIR1(1:ratio_FIR1:end);
fs_ORX1_FIR1 = fs_ORX1_DEC3/ratio_FIR1
Nsamps_ORX1_FIR1 = numel(signal_ORX1_FIR1DEC)
PLOT_FFT_dB_g(signal_ORX1_FIR1DEC, fs_ORX1_FIR1, Nsamps_ORX1_FIR1, ['signal ORX1 + FIR DEC3 + FIR1 DEC2'], 'df', 'full', 'pwr', [fnum, 4, 1, 3]);

%% input: Decimation HB1
b_HB1 = [0.000610352, 0, -0.00100708, 0, 0.00177002, 0, -0.00289917, 0, 0.004455566, 0, -0.006561279, 0, 0.009368896, 0, -0.013061523, 0, 0.017791748, 0, -0.023956299, 0, 0.031982422, 0, -0.042755127, 0, 0.057769775, 0, -0.080413818, 0, 0.119384766, 0, -0.20690918, 0, 0.632843018, 0.996490479, 0.632843018, 0, -0.20690918, 0, 0.119384766, 0, -0.080413818, 0, 0.057769775, 0, -0.042755127, 0, 0.031982422, 0, -0.023956299, 0, 0.017791748, 0, -0.013061523, 0, 0.009368896, 0, -0.006561279, 0, 0.004455566, 0, -0.00289917, 0, 0.00177002, 0, -0.00100708, 0, 0.000610352 ];
b_HB1_norm = b_HB1/b_HB1(ceil(numel(b_HB1)/2))
PLOT_FFT_dB_g(b_HB1_norm*Nsamps_ORX1_DEC3, fs_ORX1_FIR1, Nsamps_ORX1_FIR1, ['HB1'], 'df', 'full', 'pwr', [fnum, 4, 1, 3]);

%% output: HB1+Decimation(ratio_HB1 = 2)
signal_ORX1_HB1 = DSP_filter_g(b_HB1_norm, signal_ORX1_FIR1DEC, 'FD');
PLOT_FFT_dB_g(signal_ORX1_HB1, fs_ORX1_FIR1, Nsamps_ORX1_FIR1, ...
    ['signal ORX1 + FIR DEC3 + FIR1 DEC2 + HB1'], 'df', 'full', 'pwr', [fnum, 4, 1, 3]);

ratio_HB1 = 2
signal_ORX1_HB1DEC = signal_ORX1_HB1(1:ratio_FIR1:end);
fs_ORX1_HB1 = fs_ORX1_FIR1/ratio_HB1
Nsamps_ORX1_HB1 = numel(signal_ORX1_HB1DEC)
PLOT_FFT_dB_g(signal_ORX1_HB1DEC, fs_ORX1_HB1, Nsamps_ORX1_HB1, ...
    ['signal ORX1 + FIR DEC3 + FIR1 DEC2 + HB1 DEC2'], 'df', 'full', 'pwr', [fnum, 4, 1, 4]);

fnum = fnum+1
%% input: NCO2
fORX1_NCO2 = -16.268e6
t = (0:Nsamps_ORX1_HB1-1).'/fs_ORX1_HB1;

%% output: signal + NCO2
signal_ORX1_NCO2 =  signal_ORX1_HB1DEC.*exp(1i*2*pi*fORX1_NCO2*t);
PLOT_FFT_dB_g(signal_ORX1_NCO2, fs_ORX1_HB1, Nsamps_ORX1_HB1, ['signal ORX1 NCO2'], 'df', 'full', 'pwr', [fnum]);
fnum = fnum+1

%% export:
signal_Demod = signal_ORX1_NCO2;
fs_Demod = fs_ORX1_HB1
Nsamps_Demod = numel(signal_ORX1_NCO2)

%% ORX1 DL1 Demodulation:
%% Test Mode Demodulation (change reference signal's sampling rate to fs_ORX1_HB1=491.52MHz)
Config_DL1_TM = Config_DL1
Config_DL1_TM.fs = fs_Demod
ratio = Config_DL1.fs/Config_DL1_TM.fs
Config_DL1_TM.cpLengths = Config_DL1_TM.cpLengths/ratio
Config_DL1_TM.Nfft = Config_DL1_TM.Nfft/ratio
signla_ref_TM = Config_DL1_TM.signal_ref(1:ratio:end);

% channel filter
chfir_Wtype = "kaiser"
chfir_Ftype = "LPF"
chfir_Order = NaN
chfir_fTolerance = 0.5e6
chfir_K_AttdB = 60
chfir_K_fdelta = 0.5e6
chfir_fcutoffL = Config_DL1_TM.bwCarrier/2
chfir_fcutoffH = 0
chfir_Export = fs_Demod
NCarriers = 1

% output: chfir
b_ch = SYM_FIRApp(chfir_Wtype,chfir_Ftype,chfir_Order,chfir_K_AttdB,chfir_K_fdelta...
    ,chfir_fTolerance,chfir_fcutoffL,chfir_fcutoffH,df,band1UL.RFbwInband,NCarriers,chfir_Export,[fnum]);
b_ch = b_ch{:};

%% output: signal + fnco + chfir and demodulation 
N = Config_DL1_TM.RF_Ncarriers
bwChannel = Config_DL1_TM.bwChannel
fnco = ([1:N]-mean([1:N]))*bwChannel+0;
for k=1:N
    signal_Demod_k = signal_Demod.*exp(1i*2*pi*fnco(k)*t);
    signal_Demod_k_lpf = conv(signal_Demod_k, b_ch, 'same');
    rxGrid_k = OFDM_Demod_g(Config_DL1_TM, signal_Demod_k_lpf, [fnum, 1, N, k]);
    [evm_k] = dsp_evm_g(signla_ref_TM, signal_Demod_k_lpf)
    title(['carrier', num2str(k),', evm:',num2str(round(evm_k,2))])
end

disp('Why carrier2 evm is worse than others? it could be the center carrier impact by the ACLR of the other two carrier')


