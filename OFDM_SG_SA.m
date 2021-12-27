%% 2021-11-19, PA1
%% 2021-11-22, Add Band, flag_BW, flag_NCO, flag_LO, flag_fs_extend
%% 2021-11-22, Add PA for unlinearity performance
%% 2021-11-23, Change function name from OFDM_SignalGenerator to OFDM_SG_SA
%% 2021-11-24, Add flag_Demod and sigIn for demodualtion
%% 2021-11-24, Add sigGrid for evm calculation at demodulation
%% 2021-11-24, Add flag_Decimation=0, default
%% 2021

function [signal, config, demod] = OFDM_SG_SA(config, fir, PodB, Band, PA, fnum_legend, fnum, flag_Demod, sigIn)

if ~exist('flag_Demod')||isempty(flag_Demod)||any(flag_Demod==0)||strcmpi(flag_Demod,'MOD')
    flag_Demod = 0;
elseif (any(flag_Demod==1) || strcmpi(flag_Demod,'Demod')) && (exist('sigIn','var')&&~isempty(sigIn))
    flag_Demod = 1;
else
    error('vargin for demodulation?!')
end

if 0
    %% input: signal
    config.bw_Channel = '20MHz'
    config.fs = 491.52e6; % flag_fs_extend
    config.MOD = '64QAM';
    config.ratio_SamplesDMC = 'on'; % on:reduce length, off: no reduce
    config.Carrier_Type = 'NR';
    
    %% input: channel filter
    fir.Wtype = "kaiser"
    fir.Ftype = "LPF"
    fir.Order = NaN
    fir.fTolerance = -0.1e6
    fir.K_AttdB = 60
    fir.K_fdelta = 0.5e6
    fir.fcutoffL = bwChannel/2
    fir.fcutoffH = 0
    fir.Export = fs
    
    %% input: signal output power
    PodB = -15
end

if ~exist('fnum_legend','var')||isempty(fnum_legend)
    fnum_legend = 'signal';
end

if exist('PA','var')&&~isempty(PA)
    flag_PA = 1;
else
    flag_PA = 0;
end

if exist('Band','var')&&~isempty(Band)
    if isfield(Band,'bw')&&~isempty(Band.bw)
        flag_BW = 1;
        bwRF = Band.bw;
    else
        flag_BW = 0;
    end
    if isfield(Band,'fIF')&&~isempty(Band.fIF)
        flag_NCO = 1;
        fIF = Band.fIF;
    else
        flag_NCO = 0;
        fIF = 0;
    end
    if isfield(Band,'fLO')&&~isempty(Band.fLO)
        flag_LO = 1;
        fLO = Band.fLO;
        fRF = fLO+fIF;
    else
        flag_LO = 0;
    end
    
else
    flag_BW = 0;
    flag_NCO = 0;
    flag_LO = 0;
    flag_fs_extend = 0;
end

if flag_LO
    fRF_max = (bwRF/2+fIF+fLO);
    if fRF_max > config.fs/2
        flag_fs_extend = 1;
        %         config.fs = 491.52e6
        config.fs = 2^ceil(log2(fRF_max/(config.fs/2)))*config.fs; % Nyquist rule: fs/2>fRF_max
    end
end

if ~flag_Demod % SG
    
    %% generate signal
    MOD = config.MOD;
    ratio_SamplesDMC = config.ratio_SamplesDMC; % on:reduce length, off: no reduce
    Carrier_Type = config.Carrier_Type;
    [signal, sigGrid, config] = OFDM_Mod_g(config,MOD,ratio_SamplesDMC,Carrier_Type,[],[]);
    config.txGrid = sigGrid; %% 2021-11-24, Add sigGrid for evm calculation at demodulation
    
    %% output: configuration
    fs = config.fs;
    bwInband = (config.bwCarrier/2+0.5e6)*[-1 1];
    bwChannel = config.bwChannel;
    Nsamps = length(signal);
    Nbr = 1;
    df = fs/Nsamps;
    
    %% channel filter
    fir_Wtype = fir.Wtype;
    fir_Ftype = fir.Ftype;
    fir_Order = fir.Order;
    fir_fTolerance = fir.fTolerance;
    fir_K_AttdB = fir.K_AttdB;
    fir_K_fdelta = fir.K_fdelta;
    fir_fcutoffL = bwChannel/2;
    fir_fcutoffH = fir.fcutoffH;
    fir_fs = config.fs;
    NCarriers = 1;
    
    %% output: channel filter
    b_ch = SYM_FIRApp(fir_Wtype,fir_Ftype,fir_Order,fir_K_AttdB,fir_K_fdelta,fir_fTolerance,fir_fcutoffL,fir_fcutoffH,df,bwInband,NCarriers,fir_fs,[]+2);
    b_ch = b_ch{:};
    if 0
        signal_ch = conv(signal(:,1), b_ch, 'same');
    else % replace time domain convlution by fft for reducing calculation time
        [signal_ch] = DSP_filter_g(b_ch, signal, 'FD');
    end
    PLOT_FFT_dB_g(signal_ch, fs, Nsamps, ['signal + BPF'], 'df', 'full', 'pwr', []);
    
    %% gain
    % if flag_PA && isfield(PA,'gaindB')
    %    GdB = PA.gaindB;
    % end
    GdB = 0;
    vo = 10.^((PodB+GdB)/20);
    vi = ( sqrt(mean(abs(signal_ch).^2)) );
    
    %% output: gain
    signal_gain = signal_ch.*vo/vi;
    
    %% export: reference signal for demodulation
    config.signal_ref = signal_gain;
    
    %% Band: 2021-11-22
    if flag_BW||flag_NCO
        
        %% BW and NCO
        N = bwRF/bwChannel;
        config.RF_Ncarriers = N;
        if mod(N ,1)~=0
            error('BW should be N times of bwChannel')
        end
        
        fnco = ([1:N]-mean([1:N]))*bwChannel+fIF;
        signal_nco_comb = zeros(size(signal_gain));
        t = (0:Nsamps-1).'/fs;
        for n=1:N
            signal_nco = signal_gain.*exp(1i*2*pi*fnco(n)*t);
            signal_nco_comb = signal_nco+signal_nco_comb;
        end
        
        config.RF_bwChannel = bwRF;
        fnum_legend = [fnum_legend, ', bwRFMHz:',num2str(bwRF/1e6)];
    else
        signal_nco_comb = signal_gain;
    end
    
    if flag_LO
        
        %% generate LO
        %% input: LO
        LOU.leveldB = 20;
        LOU.fLO = fLO;
        
        LOU.PN = [];
        LOU.PN.ThetaDeg = 1; % %% 2020-3-14, if PN_LOU.PN_ThetaDeg==0, no-assignment
        LOU.PN.MagDriftdB1Hz = 1;
        
        LOU.PN.f_offset_Hz = 1e9*[0    0.0001    0.0010    0.0050    0.0100    0.0500    0.1000    1.0000]; % phase noise spectrum, frequencies
        LOU.PN.g_offset_dBc1Hz = [0   -80  -100  -120  -130  -130  -130  -130]; % phase noise spectrum, magnitude
        
        % ********** LO IQ Imbalance input: **********
        LOU.IMB = [];
        LOU.IMB.PhsDeg = -1.0;
        LOU.IMB.MagdB = 1;
        
        % ********** LO SPURS input: **********
        LOU.SPURS = [];
        LOU.SPURS.foffset_spurs_Hz = 1e3*[10         100        1000        1500]; % discrete spurs, freq relative to fLO
        LOU.SPURS.g_spurs_dBc1Hz = [-60   -60   -60   -60]; % discrete spurs, power relative to fLO
        
        flagT4_LOU.PhsNoise = 0; % LO wo PN
        flagT4_LOU.IMB = 1; % LO wo IQ imbalance
        flagT4_LOU.SPURS = 0; % LO wo SPURS
        flagT4_LOU.QEC = 'off';
        
        if flagT4_LOU.IMB
            disp('Add Mixer Imbalance')
        end
        %% output: LO
        [~,loU_realistic,tableT4_LOUInput,tableT4_LOU] = SYM_LOgenApp(LOU, fs, Nsamps, flagT4_LOU, [], 'half', 'LOU', []);
        
        %% upconversion
        [signal_upconversion, ~, ~] = Mixer_Up_Down_Convert_g(signal_nco_comb, loU_realistic(1,:), loU_realistic(2,:), fs, ['upconversion'],'Up',[],[]);
        config.RF_fRF = fRF;
        fnum_legend = [fnum_legend, ', fcRFMHz:',num2str(fRF/1e6)];
        
    else
        signal_upconversion = signal_nco_comb;
    end
    
    %% PA: unlinearity performance, 2021-11-22
    if flag_PA
        
        %% calculate Ipwr(Average Eneragy)
        bwInbandRF = fRF+[-bwRF bwRF]/2;
        [IpwrdBm, ~, ~] = Pwr_Inband_g(fft(signal_upconversion), fs, bwInbandRF, 5e6, 'half', [], 'dBm');
        
        %% input: PA unlinearity performance
        if 0
            paIIP3dBm = -15+50;
            paAMPMdeg = 1;
            paLinearGaindB = 0;
            paPowerUpperLimit = 45;
            paRippledB = 0; %% 2021-10-22, Add Ripple to ORX
            paSNRdBin = []; %% 2021-11-08, Add AWGN before pa
        elseif 1
            paIIP3dBm = IpwrdBm+20;
            paAMPMdeg = 5;
            paLinearGaindB = 0;
            paPowerUpperLimit = IpwrdBm+10; % PAR=10dB
            paRippledB = 0; %% 2021-10-22, Add Ripple to ORX
            paSNRdBin = []; %% 2021-11-08, Add AWGN before pa
        else
            
        end
        
        pa = DPD_PA_MemorylessNonlinearity_g('IIP3dBm',paIIP3dBm,'AMPMdeg',paAMPMdeg,...
            'LinearGaindB',paLinearGaindB,'PwrdBmLimitUpper',paPowerUpperLimit,...
            'Ripple',paRippledB,'SNRdB',paSNRdBin);
        
        signal_PA = pa(signal_upconversion);
        [IpwrdBm_PA, ~, ~] = Pwr_Inband_g(fft(signal_PA), fs, bwInbandRF, 5e6, 'half', [], 'dBm');
        config.PA_IpwrdB_Compression = IpwrdBm_PA-IpwrdBm;
        
        if 0 % debug
            signal_compression = signal_PA./signal_upconversion;
            y_shifter_PMdeg = angle(signal_compression)*180/pi;
            y_shifter_AMdeg = 20*log10(signal_compression);
            x_PindBm = 20*log10(abs(signal_upconversion));
            figure('Name','AM2PM')
            plot(x_PindBm,y_shifter_PMdeg)
            xlabel('Pin(dBm)'), ylabel('PM shifter (deg)')
            
            figure('Name','AM2AM')
            plot(x_PindBm,y_shifter_AMdeg)
            xlabel('Pin(dBm)'), ylabel('AM compression (dB)')
        end
    else
        signal_PA = signal_upconversion;
    end
    
    %% export
    signal = signal_PA;
    PLOT_FFT_dB_g(signal, fs, Nsamps, [fnum_legend], 'df', 'full', 'pwr', [fnum]);
    
else % SA
    %% 2021-11-24, Signal demodulation and analyizer process: 1.Downconversion, 2.NCO, 3.Channel filter
    %% 2021-11-24, vargin: 1.config: signal input, 2.fir: channel filter, 3.
    %% 2021-11-24, change config to signal for demodulation
    
    %% input:
    fs = config.fs;
    Nsamps = numel(sigIn);
    bwChannel = config.bwChannel;
    df = fs/Nsamps;
    
    %% Downconversion:
    if flag_LO
        % generate LO, downcoversion, perfect LO performance
        % input: LO
        LOD.leveldB = 20;
        LOD.fLO = fLO;
        
        LOD.PN = [];
        if 0
            LOD.PN.ThetaDeg = 0; % %% 2020-3-14, if PN_LOU.PN_ThetaDeg==0, no-assignment
            LOD.PN.MagDriftdB1Hz = 0;
            
            LOD.PN.f_offset_Hz = 1e9*[0    0.0001    0.0010    0.0050    0.0100    0.0500    0.1000    1.0000]; % phase noise spectrum, frequencies
            LOU.PN.g_offset_dBc1Hz = [0   -80  -100  -120  -130  -130  -130  -130]-100; % phase noise spectrum, magnitude
        end
        
        % ********** LO IQ Imbalance input: **********
        LOU.IMB = [];
        if 0
            LOU.IMB.PhsDeg = 0;
            LOU.IMB.MagdB = 0;
        end
        
        % ********** LO SPURS input: **********
        LOU.SPURS = [];
        if 0
            LOU.SPURS.foffset_spurs_Hz = 1e3*[10         100        1000        1500]; % discrete spurs, freq relative to fLO
            LOU.SPURS.g_spurs_dBc1Hz = [-60   -60   -60   -60]-100; % discrete spurs, power relative to fLO
        end
        
        flag_LOD.PhsNoise = 0; % LO wo PN
        flag_LOD.IMB = 0; % LO wo IQ imbalance
        flag_LOD.SPURS = 0; % LO wo SPURS
        flag_LOD.QEC = 'off';
        
        %% output: LO
        [~,loD_realistic,table_LODInput,table_LOD] = SYM_LOgenApp(LOD, fs, Nsamps, flag_LOD, [], 'half', 'LOD', []);
        
        %% Downconversion, BPF ?
        [signal_downconversion, ~, ~] = Mixer_Up_Down_Convert_g(sigIn, loD_realistic(1,:), loD_realistic(2,:), fs, ['downconversion'],'Down',[],[]);
    else
        signal_downconversion = sigIn;
    end
    
    %% LPF + NCO + Channel filter + Decimation
    if flag_NCO
        
        % input: LPF or channel filter
        fir_Wtype = fir.Wtype;
        fir_Ftype = fir.Ftype;
        fir_Order = fir.Order;
        fir_fTolerance = fir.fTolerance;
        fir_K_AttdB = fir.K_AttdB;
        fir_K_fdelta = fir.K_fdelta;
        fir_fcutoffL = fir.fcutoffL;
        fir_fcutoffH = fir.fcutoffH;
        fir_Export = fs;
        NCarriers = 1;
        
        flag_Decimation = 1; %% 2021-11-24, Add flag_Decimation=0, default
        
        % output: channel filter
        b_ch = SYM_FIRApp(fir_Wtype,fir_Ftype,fir_Order,fir_K_AttdB,fir_K_fdelta,fir_fTolerance,fir_fcutoffL,fir_fcutoffH,df,[],NCarriers,fir_Export,[]+2);
        b_ch = b_ch{:};
        
        % output: LPF
        disp('LPF Wide')
        fir_fTolerance = 0e6;
        fir_K_fdelta = 100e6;
        fir_fcutoffL = bwRF+abs(fIF)+fir_K_fdelta;
        b_LPF = SYM_FIRApp(fir_Wtype,fir_Ftype,fir_Order,fir_K_AttdB,fir_K_fdelta,fir_fTolerance,fir_fcutoffL,fir_fcutoffH,df,[],NCarriers,fir_Export,[]+2);
        b_LPF = b_LPF{:};
        
        %% LPF
        if 0
            signal_lpf = conv(signal_downconversion, b_LPF, 'same');
        else % replace time domain convlution by fft for reducing calculation time
            signal_lpf = DSP_filter_g(b_LPF, signal_downconversion, 'FD');
        end
        PLOT_FFT_dB_g(signal_lpf, fs, Nsamps, ['signal + LPF '], 'df', 'full', 'pwr', []);
        
        % input: NCO
        N = bwRF/bwChannel;
        if mod(N ,1)~=0
            error('BW should be N times of bwChannel')
        end
        
        fnco = -1*( ([1:N]-mean([1:N]))*bwChannel+fIF );
        signal_lpf_nco = zeros(Nsamps,N);
        t = (0:Nsamps-1).'/fs;
        for n=1:N
            
            %% NCO
            signal_lpf_nco(:,n) =  signal_lpf.*exp(1i*2*pi*fnco(n)*t);
            %% Channel filter
            disp('dmc filter or channel filter')
            if 0
                signal_lpf_nco_chfir(:,n) = conv(signal_lpf_nco(:,n), b_ch, 'same');
            else % replace time domain convlution by fft for reducing calculation time
                signal_lpf_nco_chfir(:,n) = DSP_filter_g(b_ch, signal_lpf_nco(:,n), 'FD');
            end
            PLOT_FFT_dB_g(signal_lpf_nco_chfir, fs, Nsamps, ['signal + LPF + NCO + ch fir ',num2str(n)], 'df', 'full', 'pwr', []);
        end
        
        flag_Decimation = 1; %% 2021-11-24, Add flag_Decimation=0, default
        if flag_Decimation %% 2021-11-24, bypass the decimation and demodulation directly based wide sampling rate
            % input: decimation
            switch config.bw_Channel
                case 20
                    fsOut = 30.72e6 % LTE
            end
            Ndmc = fs/fsOut;
            
            %% Decimation
            signal_lpf_nco_chfir_dmc = signal_lpf_nco_chfir(1:Ndmc:end,:);
            
            sigOut = signal_lpf_nco_chfir_dmc;
            NsampsOut = length(sigOut);
            PLOT_FFT_dB_g(signal_lpf_nco_chfir_dmc, fsOut, NsampsOut, ['signal + LPF + NCO + ch fir + DMC'], 'df', 'full', 'pwr', []);
            
            if 1 % generate configureation for demodulation
                MOD = config.MOD;
                ratio_SamplesDMC = 'on'; % on:reduce length, off: no reduce
                Carrier_Type = 'LTE';
                config2=config;
                config2.fs=[]
                [~, ~, configDMOD] = OFDM_Mod_g(config2,MOD,ratio_SamplesDMC,Carrier_Type,[],[]);
                configDMOD.txGrid = config.txGrid;
            end
            
        else
            sigOut = signal_nco_div_fir;
            fsOut = config.fs;
            NsampsOut = Nsamps;
            configDMOD = config;
        end
    end
    
    %% Demodulation
    for n=1:N
        rxGrid = OFDM_Demod_g(configDMOD, sigOut(:,n));
        fnum2 = [fnum,2,N,n];
        PLOT_FFT_dB_g(sigOut(:,n), fsOut, NsampsOut, ['signal demod.'], 'df', 'full', 'pwr', [fnum2]);
        title(['carrier ',num2str(n)])
        
        ref = configDMOD.txGrid;
        mea = rxGrid;
        [evm(n)] = dsp_evm_g(ref, mea)
        fnum2 = [fnum,2,N,N+n];
        PLOT_Constellation(rxGrid,['carrier ',num2str(n),', evm:',num2str(round(evm(n),2))],[fnum2],config.MOD,1)
    end
    
    %% export
    demod.evm = evm;
    signal = [];
    config = [];
end
end