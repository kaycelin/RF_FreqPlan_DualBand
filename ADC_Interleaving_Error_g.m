%% 2020-03-19, PA!
%% 2020-03-25, Add ADC Interleaving Offset and Mismatch impact
%% 2020-03-26, Process by ROW and Export to COLUMN

function [xo_ADC, ErrorCell, freqCell] = ADC_Interleaving_Error_g(xi, fs, ErrorVLevel, ErrorGaindB, ErrorPhaseDEG, flag_Errors, IpwrdBMin, fnum, flag_IpwrMark, fnum_legend)
%% Notice: Error = Offset = Mismatch
%% M: Numbers of ADC

% ROW
if size(xi,1)>size(xi,2) % COLUMN
    xi=xi.'; % switch to ROW
    flag_xi_original = 'COLUMN';
else
    flag_xi_original = 'ROW';
end
NCarriers=size(xi,1);

if 1
    %     disp_Err = ['Sampling by ', num2str(round(fs/1e6,2)),'MHz', ' with Interleaving Error ',flag_Errors]
    disp_title_Err = ['ADC Interleaving ', flag_Errors, ' Error']
end

if strcmp(flag_Errors,'Level')||strcmp(flag_Errors,'All')
    M = length(ErrorVLevel);
elseif strcmp(flag_Errors,'Gain')
    M = length(ErrorGaindB);
elseif strcmp(flag_Errors,'Phase')
    M = length(ErrorPhaseDEG);
elseif strcmp(flag_Errors,'NA')
    M = 1;
    disp_title_Err = ['ADC without Interleaving Error']
else
    error('check interleaving error parameters?!')
end

% if strcmp(flag_Errors,'All')
for k=1:NCarriers
    %% 2020-06-20, xo_Er is based on M-ADC
    [xo_Er(1:M,:)] = ADC_Error_g(xi(k,:), M, ErrorVLevel, ErrorGaindB, ErrorPhaseDEG, flag_Errors);
    [xo_ADC(k,:)] = ADC_interleaving_g(xo_Er, M);
    
    flag_debug=0;
    if flag_debug==1;
        PLOT_FFT_dB_g(xo_ADC, fs, [], [], 'df', 'full', 'pwr', 20200630);
    end
    
    xo_Er=[];
    % plot
    if exist('fnum','var')&&~isempty(fnum) && exist('fs','var')&&~isempty(fs)
        flag_debug='off';
        if strcmp(flag_debug,'on')
            if isrow(fnum)&&length(fnum)==4
                fnum(4)=k;
                fnum(2)=NCarriers;
            end
        end
        if ~exist('IpwrdBMin','var')&&~isempty(IpwrdBMin)
            IpwrdBMin=-174+5+40*log10(fs/length(xo_ADC));
        end
        %         PLOT_FFT_dB_g(xo_ADC(k,:), fs, length(xo_ADC), ['HD',num2str(k),', Mismatch ',flag_Errors], 'df', 'full', 'pwr',fnum);
        if exist('fnum_legend','var')&&~isempty(fnum_legend)
            freqCell(k,:) = FP_IpwrFreqCapturePlot_g(xo_ADC(k,:), fs, length(xo_ADC), [IpwrdBMin], fs/length(xo_ADC), ['min'], [fnum_legend], fnum, flag_IpwrMark);
        else
            freqCell(k,:) = FP_IpwrFreqCapturePlot_g(xo_ADC(k,:), fs, length(xo_ADC), [IpwrdBMin], fs/length(xo_ADC), ['min'], ['HD',num2str(k)], fnum, flag_IpwrMark);
        end
        
        title(disp_title_Err)
    end
end

% export and recover
if strcmp(flag_xi_original,'COLUMN')
    xo_ADC=xo_ADC.'; % switch to COLUMN
end

flag_detection = 'off';
if strcmp(flag_detection,'on')
    if strcmp(flag_Errors,'Level')||strcmp(flag_Errors,'All')
        IpwrdB_DC = 20*log10(sum((ErrorVLevel))/M);
        for j=1:M/2
            diff_ErrorVLevelj(j)=diff(ErrorVLevel(2*j-1:2*j));
        end
        IpwrdB_fsPer2 = 20*log10(abs(sum(diff_ErrorVLevelj))/M);
    elseif strcmp(flag_Errors,'Gain')
    elseif strcmp(flag_Errors,'Phase')
    end
end

for m=1:M
    if all(ErrorVLevel==0)&&all([ErrorGaindB, ErrorPhaseDEG]==0)
        ErrorCell = [];
    elseif ~all(ErrorVLevel==0)&&all([ErrorGaindB, ErrorPhaseDEG]==0)
        LevelError_cell{m,1} = [num2str(m-1),'*fs/',num2str(M)];
    elseif all(ErrorVLevel==0)&&~all([ErrorGaindB, ErrorPhaseDEG]==0)
        LevelError_cell{m,1} = [num2str(m-1),'*fs/',num2str(M)];
        GainPhaseError_cell{m,1} = [LevelError_cell{m,1},'+/-HD'];
        LevelError_cell = [];
    elseif ~all(ErrorVLevel==0)&&~all([ErrorGaindB, ErrorPhaseDEG]==0)
        LevelError_cell{m,1} = [num2str(m-1),'*fs/',num2str(M)];
        GainPhaseError_cell{m,1} = [LevelError_cell{m,1},'+/-HD'];
    end
end
% export
ErrorCell.LevelError_cell=LevelError_cell;
ErrorCell.GainPhaseError_cell=GainPhaseError_cell;

end


function [xo_Er] = ADC_Error_g(xi, M,ErrorVLevel, ErrorGaindB, ErrorPhaseDEG, flag_Errors)

% ROW
if size(xi,1)>size(xi,2) % COLUMN
    xi=xi.'; % switch to ROW
    flag_xi_original = 'COLUMN';
else
    flag_xi_original = 'ROW';
end

if size(xi,1)>size(xi,2) % COLUMN
    xi=xi'; % switch to ROW
end

xo_Er=zeros(M,length(xi));
xin = repmat(xi,M,1); % ROW

% zeros assignment
ErLevel = zeros(M,1);
ErGain = 10.^(zeros(M,1)./20);
ErPhase = exp(1i*(zeros(M,1)./180*pi));

% Error assignment
if strcmp(flag_Errors,'All')
    ErLevel = ErrorVLevel(:);
    ErGain = 10.^(ErrorGaindB(:)./10);
    ErPhase = exp(1i*(ErrorPhaseDEG(:)./180*pi));
elseif strcmp(flag_Errors,'Level')||isempty(ErrorGaindB)||isempty(ErrorPhaseDEG)
    ErLevel = ErrorVLevel(:);
elseif strcmp(flag_Errors,'Gain')||isempty(ErrorVLevel)||isempty(ErrorGaindB)
    ErGain = 10.^(ErrorGaindB(:)./10);
elseif strcmp(flag_Errors,'Phase')||isempty(ErrorVLevel)||isempty(ErrorPhaseDEG)
    ErPhase = exp(1i*(ErrorPhaseDEG(:)./180*pi));
    ErPhase = ErrorPhaseDEG(:)./180*pi;
end

% Mismatch
if 0
    xo_Er(:,:)=((xin+ErLevel).*ErGain).*exp(1i*(ErPhase));
else % 2021-12-23, updated
    xo_Er(:,:)=((xin+ErLevel).*ErGain).*ErPhase;
end

% export and recover
if strcmp(flag_xi_original,'COLUMN')
    xo_Er=xo_Er.'; % switch to COLUMN
end
end


function [xo] = ADC_interleaving_g(xi, M)

% ROW
if size(xi,1)>size(xi,2) % COLUMN
    xi=xi.'; % switch to ROW
    flag_xi_original = 'COLUMN';
else
    flag_xi_original = 'ROW';
end

Nsamps = length(xi);
% n=[1:Nsamps];
x_interleaving = zeros(M,Nsamps);
% ind_n_interleaving = zeros(M,Nsamps);
if mod(log2(M),1)~=0
    error('M')
end
for m=1:M
    x_interleaving(m,m:M:end) = xi(m,m:M:end);
end
xo = sum(x_interleaving,1); % ROW
% export and recover
if strcmp(flag_xi_original,'COLUMN')
    xo=xo.'; % switch to COLUMN
end
end

