% [waveformOFDM, data_grid, Config] = genOFDM_MOD_g(bw_Channel,MOD,ratio_SamplesDMC,fnum,Carrier_Type)
%% 2020-09-18, P1A
clear all
close all

%% Inputs:
Scs_kHz = 30;
T = 10e-3;
MOD = '256QAM'
BW_MHz = 20
Carrier_Type = 'NR'
fnum=101

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
        end
        df_Scs = SCS_kHz*1e3;
        Nslots = 2;
        Nsymbols = 14;
        % NR100, default
        fs_NR100 = 122.88e6;
        CPLengths_NR100 = [352,repmat(288,[1,Nsymbols-1])];
        
        %         if SCS_kHz == 30
        %             if ~exist('BW_MHz','var')||isempty(BW_MHz)||BW_MHz == 100||strcmp(BW_MHz,'100MHz')
        %                 % NR100, default
        %                 BW_MHz = 100;
        %                 fs = fs_NR100;
        %                 NRB = 273;
        %             elseif BW_MHz == 20||strcmp(BW_MHz,'20MHz')
        %                 fs = 23.04e6;
        %                 NRB = 51;
        %             elseif BW_MHz == 40||strcmp(BW_MHz,'40MHz')
        %                 fs = 30.72e6;
        %                 NRB = 106;
        %             end
        %         end
        
        %% generate table_NRB_SCSvsBW
        BWMHz = [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100].';
        fsMHz = [7.68, nan, nan, 23.04, nan, nan, nan, nan, nan, nan, 122.88].';
        NRBofScs_15kHz = [25, 52, 79, 106, 133, 160, 216, 277, nan, nan, nan].';
        NRBofScs_30kHz = [11, 24, 38, 51, 65, 78, 106, 133, 162, 217, 273].';
        NRBofScs_60kHz = [nan, 11, 18, 24, 31, 38, 51, 65, 79, 107, 135].';
        table_BW_fs_NRBofScs = table(BWMHz,fsMHz, NRBofScs_15kHz,NRBofScs_30kHz,NRBofScs_60kHz)
        
        if ischar(BW_MHz)
            BW_MHz = str2double(erase(BW_MHz,'MHz'));
        end
        ind_BWMHz = find(BWMHz==BW_MHz);
        fs = fsMHz(ind_BWMHz)*1e6;
        if SCS_kHz == 30
            NRB = table_BW_fs_NRBofScs.NRBofScs_30kHz(ind_BWMHz);
        elseif SCS_kHz == 15
            NRB = table_BW_fs_NRBofScs.NRBofScs_15kHz(ind_BWMHz);
        elseif SCS_kHz == 60
            NRB = table_BW_fs_NRBofScs.NRBofScs_60kHz(ind_BWMHz);
        end
        
        NFFT = fs/df_Scs;
        CPLengths = CPLengths_NR100.*(fs/fs_NR100);
        
    case 'LTE'
        % LTE20
        if ~exist('SCS_kHz','var')||isempty(SCS_kHz)
            SCS_kHz = 15;
        end
        df_Scs = SCS_kHz*1e3;
        Nslots = 2;
        Nsymbols = 7;
        fs_LTE20 = 30.72e6;
        NFFT_LTE20 = fs_LTE20/df_Scs;
        CPLengths_NR100 = [160,repmat(160,[1,Nsymbols-1])];
        BW_MHz = 20;
        Nslots = 2;
        Nsymbols = 7;
        
        if ~exist('BW_MHz','var')||isempty(BW_MHz)||BW_MHz == 20||strcmp(BW_MHz,'20MHz')
            fs = fs_LTE100;
            NRB = 100;
        elseif BW_MHz == 5||strcmp(BW_MHz,'5MHz')
            fs = 7.68e6;
            NRB = 25;
        end
        NFFT = NFFT_NR100*(fs/fs_NR100);
        CPLengths = CPLengths_NR100.*(fs/fs_NR100);
end

NScsPerRB = 12;
NScs = NRB*NScsPerRB;
bw_Carrier = NScs*df_Scs;
NOFDMsymbols = Nsymbols*Nslots*Nsubframes;
Nsamps = fs*T;
if Nsamps ~= sum(repmat(NFFT, [1,Nsymbols])+CPLengths)*Nslots*Nsubframes
    error('Check the NFFT and CPLengths!')
end

%% 2020-03-09, Reduce the Lengths of data
if  ~exist('ratio_SamplesDMC','var')||isempty(ratio_SamplesDMC)||strcmp(ratio_SamplesDMC,'off')
    ratio_R=1;
elseif strcmp(ratio_SamplesDMC,'on')
    ratio_R = gcd(CPLengths(1),CPLengths(2));
elseif isnumeric(ratio_SamplesDMC)
    ratio_R=ratio_SamplesDMC;
end

if ratio_R~=1
    %     data_binary = data_binary(1:ratio_R:end,:);
    df_Scs = df_Scs*ratio_R;
    Nsamps = Nsamps/ratio_R;
    NFFT = NFFT/ratio_R;
    NScs = NScs/ratio_R;
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
if mod(NFFT,2)==0
    f = (-NFFT/2:NFFT/2-1)*df_Scs;
elseif mod(NFFT,2)==1
    f = [-(NFFT-1)/2:(NFFT-1)/2]*df_Scs;
end

if mod(NScs,2)==0
    f_Scs = [(-NScs/2:-1),(1:NScs/2)]*df_Scs; % exclude DC
elseif mod(NScs,2)==1
    f_Scs = [(-NScs-1)/2:(NScs-1)/2]*df_Scs; % include DC
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
    data_Map2NFFT = zeros(NFFT,NAnts);
    data_Map2NFFT(ind_Scs,1) = data_grid(:,n,NAnts);
    
    % ifft
    data_ifft = [];
    data_ifft = ifft(fftshift(data_Map2NFFT));
    
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
if length(fnum)==4
    PLOT_FFT_dB_g(waveformOFDM, fs, length(waveformOFDM), [Carrier_Type,num2str(BW_MHz),'MHz'], 'df', 'full', 'pwr', [fnum(1) fnum(2) 2 2*fnum(4)]);
else
    %     PLOT_Constellation(data_grid, [fnum(1) 1 2 1])
    PLOT_FFT_dB_g(waveformOFDM, fs, length(waveformOFDM), [Carrier_Type,num2str(BW_MHz),'MHz'], 'df', 'full', 'pwr', [fnum(1) 1 2 2]);
end
PARdB = CCDF_g(waveformOFDM, Nsamps, []);

%% 2020-01-22, Export configuration
Config = struct('bw_Channel',{BW_MHz},'bwChannel',{BW_MHz*1e6},'fs',{fs},'bwCarrier',{bw_Carrier},'cpLengths',{CPLengths},'nFFT',{NFFT},'nRB',{NRB},'nSubcarriers',{NScs},'dfHzSubcarrier',{df_Scs},'SamplesDecimation',{ratio_R},'PARdB',{PARdB});
