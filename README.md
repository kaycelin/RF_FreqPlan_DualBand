# RF_FreqPlan_DualBand

DRAFT
Background: 
- Freq plan for IF architecture with 
- Only One LO freq. for TX and RX used, used NCO to shift to zero freq
- The ADC BW is limited
- Dualband: B1 and B3 FDD mode

# Generate waveform for B1 and B3 FDD mode
  generate waveform based on RXIBW 491.52MHz
  - input: signal
    - config.bw_Channel = '20MHz'
    - config.fs = 491.52e6*1; % flag_fs_extend
    - config.MOD = '64QAM';
    - config.ratio_SamplesDMC = 'on'; % on:reduce length, off: no reduce
    - config.Carrier_Type = 'LTE'

 - input: channel filter
    - fir.Wtype = "kaiser"
    - fir.Ftype = "LPF"
    - fir.Order = NaN
    - fir.fTolerance = -0.1e6
    - fir.K_AttdB = 60
    - fir.K_fdelta = 0.5e6
    - fir.fcutoffL = 20e6/2 % based on config.bw_Channel = '20MHz'
    - fir.fcutoffH = 0
    - fir.fs = config.fs

  - input: signal output power
    - SNRdB_DL2UL = 60
    - LeakagedB_DL2UL = -60
    - PodB_UL = -15
    - PodB_DL = PodB_UL + SNRdB_DL2UL + LeakagedB_DL2UL % SNRdB_DL2UL = 60, Leakage_DL2UL = -60

**_Setting the signle LO frequency of IF architectue for B1 and B3 FDD mode, the LO freq. could be tunning for freq. plan_**
  - input: signal with band information and signle fLO!!!
    - fLO = 1900e6; %% 2021-11-22, PA2, input: signal with band information and signle fLO!!!
    
    - band3DL.band = 'B3'
    - band3DL.bw = 60e6;
    - band3DL.fLO = fLO;
    - band3DL.RFbwInband = [2110 2170]*1e6
    - band3DL.fIF = mean(band3DL.RFbwInband) - band3DL.fLO

    - band3UL.band = 'B3'
    - band3UL.bw = 60e6;
    - band3UL.fLO = fLO;
    - band3UL.RFbwInband = [1920 1980]*1e6
    - band3UL.fIF = mean(band3UL.RFbwInband) - band3UL.fLO

    - band1DL.band = 'B1'
    - band1DL.bw = 80e6;
    - band1DL.fLO = fLO;
    - band1DL.RFbwInband = [1800 1880]*1e6
    - band1DL.fIF = mean(band1DL.RFbwInband) - band1DL.fLO

    - band1UL.band = 'B1'
    - band1UL.bw = 80e6;
    - band1UL.fLO = fLO;
    - band1UL.RFbwInband = [1710 1790]*1e6
    - band1UL.fIF = mean(band1UL.RFbwInband) - band1UL.fLO

  - input: PA performance
     - PA.flag = 'on'
 
  - output: signal (B1+B3 FDD waveform)
![image](https://user-images.githubusercontent.com/87049112/143376793-9d1fefa2-5987-4e5e-810e-6fc84f19c938.png)

# Case1: Evaluate the DL leakage impact UL performance
  1. signal combination
    - signal_combToUL1
    - signal_combToUL3

  2. BPF for UL signal
  - input: BPF1
    - BPF1_fir_Wtype = "kaiser"
    - BPF1_fir_Ftype = "BPF"
    - BPF1_fir_Order = NaN
    - BPF1_fir_fTolerance = 10e6
    - BPF1_fir_K_AttdB = 40
    - BPF1_fir_K_fdelta = 10e6
    - BPF1_fir_fcutoffL = band1UL.RFbwInband(1)
    - BPF1_fir_fcutoffH = band1UL.RFbwInband(2)
    - BPF1_fir_Export = fsOut
    - NCarriers = 1
    
- output: signal + BPF
![image](https://user-images.githubusercontent.com/87049112/143377144-dc7a4b84-c4ec-4bae-a10e-9574716ae0cc.png)
![image](https://user-images.githubusercontent.com/87049112/143377591-4879e9ce-de6d-4bbc-9817-eea7a5ed7b1c.png)

- output: UL demodulation
  - SNR:10dB (UL1 to DL3's image)
  - UL1 demodulation
![image](https://user-images.githubusercontent.com/87049112/143378064-f706f010-5fe0-4001-a483-eb9416cc0711.png)
  - UL3 demodulation
![image](https://user-images.githubusercontent.com/87049112/143378182-3d5ddde0-607f-48ec-9c61-7a45442cf7d3.png)
  - SNR:27dB (UL1 to DL3's image)
  - UL1 demodulation
![image](https://user-images.githubusercontent.com/87049112/143387597-b744a83d-46a2-4ae8-9d99-c3bf510cee8d.png)
  - UL3 demodulation
![image](https://user-images.githubusercontent.com/87049112/143387626-03088424-e85a-4795-b2aa-e1b800943190.png)

**Summary:**        
- **1. The parameters of _DL leakage(LeakagedB_DL2UL)_ will impact the UL1 demodulation performance, the DL3's image overlape with UL1 signal, the EVM of UL1 is up to 18%(need to improvement)**       
- **1b. SNR 10dB introudce to UL1 EVM 18%, SNR 27dB introduce to UL1 EVM 5.21%**      
- **2. The parameters of _DL leakage(LeakagedB_DL2UL) and UL BPF_ will impact the UL3 demodulation performance, the DL3's aclr overlape with UL3 signal, the EVM of UL3 is 4.68% at UL3's high freq. part**     

# Case2: Evaluate the DL signal to ORX, ORX is RF sampling architecture     
on-goin...
