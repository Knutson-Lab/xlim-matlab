function [taucomps,thetacomps,dtaucomps,dthetacomps] = ...
    tcspc_dk_delta_function_conv(tcal,ndata,irf,param,ntau,ntheta,taur,convflag)
    %TCSPC_DK_DELTA_FUNCTION_CONV Convolves function using delta function
    %% First generate convolution components
    switch convflag
        case "FFT"
            [taucomps,thetacomps,dtaucomps,dthetacomps,ddtaucomps,ddthetacomps] = ...
                tcspc_dk_conv_fft(tcal,ndata,irf,param,ntau,ntheta);
    end
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
    %% Generate effective lifetime curves
    efftau = permute(((1 ./ tau) + (1 ./ theta')) .^ (-1),[3 1 2]);
    %% Scale convolved components appropriately
    taucomps = irf + ((1 / taur) - (1 ./ tau)) .* taucomps;
    dtaucomps = (dtaucomps / taur) + ddtaucomps;
    thetacomps = irf + ((1 / taur) - (1 ./ efftau)) .* thetacomps;
    dthetacomps = (dthetacomps / taur) + ddthetacomps;
end

