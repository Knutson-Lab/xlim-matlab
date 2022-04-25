function [taucomps,thetacomps,dtaucomps,dthetacomps,ddtaucomps,ddthetacomps] = ...
    tcspc_dk_conv_fft(tcal,ndata,irf,param,ntau,ntheta)
    %TCSPC_DK_CONV_FFT Convolves IRF with parameters to generate new curve

    %% Make time axis
    tax = (tcal / 2):tcal:(ndata * tcal);
    %% Isolate parameters from param
    if ntau > 0
        tau = param(2:2:(2 * ntau));
    else
        tau = [];
    end
    if ntheta > 0
        offs = (2 * ntau) + 1;
        theta = param((offs + 1):2:(offs + (2 * ntheta)));
    else
        theta = [];
    end
    %% Generate basic exponential curve
    if ntau > 0
        taucomps = exp(-tax ./ tau)';
        dtaucomps = (tax ./ (tau .^ 2))' .* taucomps;
        ddtaucomps = ((tax - tau) ./ (tau .^ 3))' .* taucomps;
    else
        taucomps = [];
        dtaucomps = [];
        ddtaucomps = [];
    end
    if ntheta > 0
        if ntau > 0
            efftau = permute(((1 ./ tau) + (1 ./ theta')) .^ (-1),[3 1 2]);
            thetacomps = exp(-tax' ./ efftau);
            dthetacomps = (tax ./ (efftau .^ 2))' .* thetacomps;
            ddthetacomps = ((tax - efftau) ./ (efftau .^ 3))' .* thetacomps;
        else
            thetacomps = [];
            dthetacomps = [];
            ddthetacomps = [];
        end
    else
        thetacomps = [];
        dthetacomps = [];
        ddthetacomps = [];
    end
    
    %% Convolve with IRF
    tempirf = fft(irf,[],1);
    if ntau > 0
        taucomps = abs(ifft(fftshift(fft(taucomps,[],1) .* tempirf,1),[],1));
        dtaucomps = abs(ifft(fftshift(fft(dtaucomps,[],1) .* tempirf,1),[],1));
        ddtaucomps = abs(ifft(fftshift(fft(ddtaucomps,[],1) .* tempirf,1),[],1));
    end
    if ntheta > 0 && ntau > 0
        thetacomps = abs(ifft(fftshift(fft(thetacomps,[],1) .* tempirf,1),[],1));
        dthetacomps = abs(ifft(fftshift(fft(dthetacomps,[],1) .* tempirf,1),[],1));
        ddthetacomps = abs(ifft(fftshift(fft(ddthetacomps,[],1) .* tempirf,1),[],1));
    end
    
end

