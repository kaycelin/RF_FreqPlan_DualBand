%% 2020-09-18, P1A
%% 2021-01-15, Add ConfigIn
%% 2021-02-26, 1RB signal
%% 2021-07-15, Support Sampling rate:fs setup
%% 2021-11-01, flag_gain_comp_ifft, gain compensation for ifft
%% 2021-11-03, MOD updated
%% 2021-11-04, export pwr_grid for demodulation and normalized

function [waveformOFDM, data_grid, Config] = OFDM_Mod_g(ConfigIn,MOD,ratio_SamplesDMC,Carrier_Type,fnum,fnum_save_dir)
% clear all
% close all

%% Inputs:
% Scs_kHz = 30;
% T = 10e-3;
% MOD = '256QAM'
% BW_MHz = 20
% Carrier_Type = 'NR'
% fnum=101

%% 2021-01-15, Add ConfigIn
if isstruct(ConfigIn)
    if isfield(ConfigIn,'bw_Channel')
        bw_Channel = ConfigIn.bw_Channel;
    end
    if isfield(ConfigIn,'NRB')
        NRB = ConfigIn.NRB;
    end
    if isfield(ConfigIn,'Scs_kHz')
        SCS_kHz = ConfigIn.Scs_kHz;
    end
    if isfield(ConfigIn,'fs')
        fs = ConfigIn.fs;
    end
    %     if ~exist('bw_Channel','var')||isempty(bw_Channel)
    %         error('check the input of ConfigIn!')
    %     end
else
    bw_Channel = ConfigIn;
end

% Default parameters
T_10ms = 10e-3;
if ~exist('T','var')||isempty(T)
    T = T_10ms;
    % else
    %     error('T is NOT based on 1ms?')
end
Nsubframes = T/1e-3; % 1ms/Subframe

%% Initial parameters
if ~exist('Carrier_Type','var')||isempty(Carrier_Type)
    Carrier_Type = 'LTE';
end

switch Carrier_Type
    case 'NR'
        if ~exist('SCS_kHz','var')||isempty(SCS_kHz)
            SCS_kHz = 30;
            Nslots = 2;
            Nsymbols = 14;
            CPLengths_NR100 = [352,repmat(288,[1,Nsymbols-1])];
            BWMHz = [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100].';
            fsMHz = [7.68, nan, nan, 23.04, nan, nan, nan, nan, nan, nan, 122.88].';
        else
            SCS_kHz = 15;
            Nslots = 1;
            Nsymbols = 14;
            CPLengths_NR100 = repmat([320,repmat(288,[1,Nsymbols/2-1])],[1,2]);
            BWMHz = [5, 10, 15, 20, 25, 30, 40, 50, nan, nan, nan].';
            fsMHz = [7.68, nan, nan, 23.04, nan, nan, nan, 61.44, nan, nan, nan].';
        end
        CPLengths_Base = CPLengths_NR100;
        
        %% generate table_NRB_SCSvsBW
        NRBofScs_15kHz = [25, 52, 79, 106, 133, 160, 216, 277, nan, nan, nan].';
        NRBofScs_30kHz = [11, 24, 38, 51, 65, 78, 106, 133, 162, 217, 273].';
        NRBofScs_60kHz = [nan, 11, 18, 24, 31, 38, 51, 65, 79, 107, 135].';
        table_BW_fs_NRBofScs = table(BWMHz,fsMHz, NRBofScs_15kHz,NRBofScs_30kHz,NRBofScs_60kHz)
        
    case 'LTE'
        SCS_kHz = 15;
        Nslots = 2;
        Nsymbols = 7;
        CPLengths_LTE20 = [160,repmat(144,[1,Nsymbols-1])];
        CPLengths_Base = CPLengths_LTE20; % LTE20, default

        
        %% generate table_LTE_SCSvsBW
        BWMHz = [1.4, 3, 5, 10, 15, 20].';
        fsMHz = [1.92, 3.84, 7.68, 15.36, 23.04, 30.72].';
        NRBofScs_15kHz = [6, 15, 25, 50, 75, 100].';
        table_BW_fs_NRBofScs = table(BWMHz,fsMHz, NRBofScs_15kHz)

end

if ~exist('bw_Channel','var')||isempty(bw_Channel)
    BW_MHz = max(BWMHz);
elseif ischar(bw_Channel)
    BW_MHz = str2double(erase(bw_Channel,'MHz'));
else
    BW_MHz=bw_Channel;
end

%% 2021-02-26, 1RB signal
if exist('NRB','var') && ~isempty(NRB) && NRB == 1
    disp('1RB signal')
    BW_MHz = 0.18;
else
    ind_BWMHz = find(BWMHz==BW_MHz);
    fs_carrier = fsMHz(ind_BWMHz)*1e6;
    if ~exist('fs','var')||isempty(fs)
        fs = fs_carrier;
    end
    if SCS_kHz == 30
        NRB = table_BW_fs_NRBofScs.NRBofScs_30kHz(ind_BWMHz);
    elseif SCS_kHz == 15
        NRB = table_BW_fs_NRBofScs.NRBofScs_15kHz(ind_BWMHz);
    elseif SCS_kHz == 60
        NRB = table_BW_fs_NRBofScs.NRBofScs_60kHz(ind_BWMHz);
    end
end

%% Export parameters
df_Scs = SCS_kHz*1e3;
fs_Max = max(fsMHz)*1e6;
Nfft = fs/df_Scs;
CPLengths = CPLengths_Base.*(fs/fs_Max);

NScsPerRB = 12;
NScs = NRB*NScsPerRB;
bw_Carrier = NScs*df_Scs;
NOFDMsymbols = Nsymbols*Nslots*Nsubframes;
Nsamps = fs*T;

if fs == fs_carrier
    if Nsamps ~= sum(repmat(Nfft, [1,Nsymbols])+CPLengths)*Nslots*Nsubframes
        error('Check the Nfft and CPLengths!')
    end
end

%% 2020-03-09, Reduce the Lengths of data
if  ~exist('ratio_SamplesDMC','var')||isempty(ratio_SamplesDMC)||strcmp(ratio_SamplesDMC,'off')
    ratio_R=1;
elseif strcmp(ratio_SamplesDMC,'on')
    ratio_R = gcd(CPLengths(1),CPLengths(2));
    
    while any(rem([Nsamps/ratio_R, Nfft/ratio_R, NScs/ratio_R],1)~=0)
        ratio_R = ratio_R/2;
    end
elseif isnumeric(ratio_SamplesDMC)
    ratio_R=ratio_SamplesDMC;
end

if ratio_R~=1
    %     data_binary = data_binary(1:ratio_R:end,:);
    df_Scs = df_Scs*ratio_R;
    Nsamps = Nsamps/ratio_R;
    Nfft = Nfft/ratio_R;
    NScs = NScs/ratio_R;
    if mod(NScs,1)~=0
        NScs = 1;
    end
    CPLengths = CPLengths/ratio_R;
end

%% Generate data_binary from random matrix [NScs x NOFDMsymbols]
if ~exist('MOD','var')||isempty(MOD)
    MOD = '64QAM'; % default
end
switch MOD
    case 'QPSK'
        M=4;
    case '64QAM'
        M=64;
    case '256QAM'
        M=256;
end

data = randi([0, M-1], NScs, NOFDMsymbols);
data_binary = cellstr(dec2bin(data,log2(M))); %*****
data_binary = reshape(data_binary,NScs,NOFDMsymbols);

%% Constellation Mapping for data_binary
data_grid = zeros(NScs,NOFDMsymbols);
M_bit = dec2bin([0:M-1],log2(max(M)));
for i=1:M
    ind = find(ismember(data_binary,M_bit(i,:))==1);
    
    switch MOD
        case '256QAM'
            Bit0=str2num(M_bit(i,log2(M)-0));
            Bit1=str2num(M_bit(i,log2(M)-1));
            Bit2=str2num(M_bit(i,log2(M)-2));
            Bit3=str2num(M_bit(i,log2(M)-3));
            Bit4=str2num(M_bit(i,log2(M)-4));
            Bit5=str2num(M_bit(i,log2(M)-5));
            Bit6=str2num(M_bit(i,log2(M)-6));
            Bit7=str2num(M_bit(i,log2(M)-7));
            d = 1/sqrt(170)*[(1-2*Bit0)*[8-(1-2*Bit2)*[4-(1-2*Bit4)*[2-(1-2*Bit6)]]]+ ...,
                1j*(1-2*Bit1)*[8-(1-2*Bit3)*[4-(1-2*Bit5)*[2-(1-2*Bit7)]]]];
        case '64QAM'
            Bit0=str2num(M_bit(i,log2(M)-0));
            Bit1=str2num(M_bit(i,log2(M)-1));
            Bit2=str2num(M_bit(i,log2(M)-2));
            Bit3=str2num(M_bit(i,log2(M)-3));
            Bit4=str2num(M_bit(i,log2(M)-4));
            Bit5=str2num(M_bit(i,log2(M)-5));
            d = 1/sqrt(42)*[(1-2*Bit0)*[4-(1-2*Bit2)*[2-(1-2*Bit4)]]+1j*(1-2*Bit1)*[4-(1-2*Bit3)*[2-(1-2*Bit5)]]];
        case 'QPSK'
            Bit0=str2num(M_bit(i,log2(M)-0));
            Bit1=str2num(M_bit(i,log2(M)-1));
            d = 1/sqrt(2)*[(1-2*Bit0)+1j*(1-2*Bit1)];
    end
    
    data_grid(ind) = d;
end
pwr_grid = mean(abs(data_grid(:,1)).^2); %% 2021-11-04, export pwr_grid for demodulation and normalized

% plot Constellation
if exist('fnum','var')&&~isempty(fnum)
    %     PLOT_Constellation(data_grid, [MOD],200*fnum)
    if length(fnum)==4
        %         PLOT_Constellation(data_grid, [bw_Channel,', ',MOD], [fnum(1) fnum(2) 2 2*fnum(4)-1])
        PLOT_Constellation(data_grid, [MOD], [fnum(1) fnum(2) 2 2*fnum(4)-1])
    else
        %         PLOT_Constellation(data_grid, [bw_Channel,', ',MOD], [fnum(1) 1 2 1])
        PLOT_Constellation(data_grid, [MOD], [fnum(1) 1 2 1])
    end
    hold off
end

%% OFDM Modulation & Constellation Mapping
if mod(Nfft,2)==0
    f = (-Nfft/2:Nfft/2-1)*df_Scs;
elseif mod(Nfft,2)==1
    f = [-(Nfft-1)/2:(Nfft-1)/2]*df_Scs;
end

if mod(NScs,2)==0
    f_Scs = [(-NScs/2:-1),(1:NScs/2)]*df_Scs; % exclude DC
elseif mod(NScs,2)==1
    f_Scs = [(-NScs+1)/2:(NScs-1)/2]*df_Scs; % include DC
end
ind_Scs = find(ismember(f,f_Scs)==1);
if length(f_Scs)~=size(data_grid,1)
    error('Check the f_Scs or data_grid!')
else
    ind_Scs = find(ismember(f,f_Scs)==1);
end

NAnts = 1;
ind_NextSymbol = 1;
waveformOFDM = [];
for n=1:NOFDMsymbols
    
    % Mapping data_grid into subcarrier
    data_Map2Nfft = zeros(Nfft,NAnts);
    data_Map2Nfft(ind_Scs,1) = data_grid(:,n,NAnts);
    
    % ifft
    data_ifft = [];
    data_ifft = ifft(fftshift(data_Map2Nfft))/ratio_R;
    flag_gain_comp_ifft = 1; %% 2021-11-01, flag_gain_comp_ifft, gain compensation for ifft
    if flag_gain_comp_ifft
        data_ifft = data_ifft*sqrt(size(data_Map2Nfft,1));
    end
    
    % Extend data with CPLength
    CPLengths_n = CPLengths(mod(n-1,length(CPLengths))+1);
    data_CP = data_ifft(end-CPLengths_n+1:end);
    data_CPExtend = [data_CP; data_ifft];
    dataLength = length(data_CPExtend);
    
    % IQdata Export
    waveformOFDM(ind_NextSymbol:(ind_NextSymbol-1)+dataLength,NAnts) = data_CPExtend;
    ind_NextSymbol = ind_NextSymbol+dataLength;
end
if length(waveformOFDM)~=Nsamps
    error('Check the Nsamps!')
end

% plot
if exist('fnum','var')&&~isempty(fnum)
    
    if length(fnum)==4
        PLOT_FFT_dB_g(waveformOFDM, fs, length(waveformOFDM), [Carrier_Type,num2str(BW_MHz),'MHz'], 'df', 'full', 'pwr', [fnum(1) fnum(2) 2 2*fnum(4)]);
    else
        %     PLOT_Constellation(data_grid, [fnum(1) 1 2 1])
        PLOT_FFT_dB_g(waveformOFDM, fs, length(waveformOFDM), [Carrier_Type,num2str(BW_MHz),'MHz'], 'df', 'full', 'pwr', [fnum(1) 1 2 2]);
    end
    %% 2020-10-31, fnum_save_dir: save picture to folder
    if exist('fnum_save_dir')&&~isempty(fnum_save_dir)
        fnum_save_file = [fnum_save_dir,'\',num2str(fnum),'.fig']
        saveas(gcf,[fnum_save_file])
    end
end
PARdB = CCDF_g(waveformOFDM, Nsamps, []);

%% 2020-01-22, Export configuration -> 2021-11-03, MOD updated
Config = struct('MOD',{MOD},'PwrGrid',{pwr_grid},'bw_Channel',{BW_MHz},'bwChannel',{BW_MHz*1e6},'fs',{fs},'bwCarrier',{bw_Carrier},'cpLengths',{CPLengths},'Nfft',{Nfft},'NRB',{NRB},'NScs',{NScs},'dfScs',{df_Scs},'SamplesDecimation',{ratio_R},'PARdB',{PARdB});
end