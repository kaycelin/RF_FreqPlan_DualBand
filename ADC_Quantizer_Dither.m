function [xout, disp_legend_export] = quant(xi,LSB,scale_DitherLSB,flag_DCRomve,fs,fnum,fnum_legend)
% close all

if ~exist('fnum_legend','var')||isempty(fnum_legend)
    disp_legend = [];
else
    disp_legend = [', ',fnum_legend];
end

if isrow(xi)
    flag_row = 1;
    xi = xi.';
else
    flag_row = 0;
end

if ~exist('LSB','var')||isempty(LSB)
    Vref = 2.5;
    nbits = 4;
    LSB = Vref/2^nbits;
    disp_LSB = ['Vref:',num2str(Vref),', nbits:',num2str(nbits)];
elseif numel(LSB)==2
    Vref = LSB(1);
    nbits = LSB(2);
    disp_LSB = ['Vref:',num2str(Vref),', nbits:',num2str(nbits)];
    LSB = Vref/2^nbits;
elseif numel(LSB)==1
    disp_LSB = ['LSB:',num2str(LSB)];
else
    error('LSB?!')
end

if ~exist('scale_DitherLSB','var')||isempty(scale_DitherLSB)
    Dither = (1+1i).*zeros(size(xi));
    DitherLSBstr = 'off';
    
else 
    DitherLSB = LSB/scale_DitherLSB;
    Dither = (1+1i).*ones(size(xi))*(DitherLSB);
    DitherLSBstr = num2str(DitherLSB);

end

method_Quantize = '1'
switch method_Quantize
    case '1'
        xiD = xi;
        err_quan = (0.5+0.5i).*ones(size(xiD));
        NDither = 10;
        Dither_k = zeros(size(xiD,1), size(xiD,2), NDither);
        xout_k = zeros(size(xiD,1), size(xiD,2), NDither);
        for k=1:NDither
            Dither_k(:,:,k) = randi([-1,1],size(xiD)).*Dither;     
            
            xout_k(:,:,k) = LSB * floor(1/LSB * (xiD+Dither_k(:,:,k)) + err_quan);            

            % Dither Noise Remove
           flag_DitherNoiseRemove = 1;
            if flag_DitherNoiseRemove
                xiD = xout_k(:,:,k)-1*Dither_k(:,:,k);
            else
                xiD = xout_k(:,:,k);
            end
        end
%         xout_k(:,:,k) = xout_k(:,k)-1*Dither_k(:,:,k);

        flag_debug = 0;
        if flag_debug
            figure(20210417)
            plot(real(xi(:,1))), hold on
            plot(real(xout_k(:,1,10))), hold on
        end
        
        % Noise Cancellation
        if ~flag_DitherNoiseRemove
            xout = mean(xout_k-Dither_k,3);
        else
            % Average
            xout = mean(xout_k,3);
        end
        % DC Remove
        for idBR=1:size(xi,2)
            [xout_RemoveDC(:,idBR), gcorr(idBR)] = dsp_normalize_RemoveDC_krait_g(xi(:,idBR), xout(:,idBR), 1);
        end
        xout = xout_RemoveDC;
        stop=1;

    case '2'
        
        VoutLSB = 0:LSB:Vref;

        xiD = xi + Dither;
        xi_I = real(xiD);
        xi_Q = imag(xiD);
        xoADC_I = zeros(size(xiD));
        xoADC_Q = zeros(size(xiD));
        
        for n=1:2^nbits
            ind_I_n = find(abs(xi_I) >= VoutLSB(n)-LSB/2 & abs(xi_I) < VoutLSB(n)+LSB/2);
            xoADC_I(ind_I_n) = sign(xi_I(ind_I_n))*VoutLSB(n);
        end
        
        if ~isreal(xi)
            for n=1:2^nbits
                ind_Q_n = find(abs(xi_Q) >= VoutLSB(n)-LSB/2 & abs(xi_Q) < VoutLSB(n)+LSB/2);
                xoADC_Q(ind_Q_n) = sign(xi_I(ind_Q_n))*VoutLSB(n);
            end
        end
        
        xout2 = xoADC_I + 1i*xoADC_Q;
        
end

if flag_row
    xout = xout.';
end

if exist('fs','var')&&~isempty(fs)
    if ~exist('fnum','var')||isempty(fnum)
        fnum = 20010416
    end
%     PLOT_FFT_dB_g(xi, fs, length(xi), ['Input'], 'df', 'full', 'pwr', [20010416], [], []);
%     PLOT_FFT_dB_g(xout, fs, length(xout), ['Quantiztion, LSB:',num2str(LSB),', DitherLSB:',(DitherLSBstr)], 'df', 'full', 'pwr', [fnum], [], []);
disp_legend_export = ['ADC Quantiztion, ',disp_LSB,', DitherLSB:',(DitherLSBstr),disp_legend];
PLOT_FFT_dB_g(xout, fs, length(xout), disp_legend_export, 'df', 'full', 'pwr', [fnum], [], []);
end

end