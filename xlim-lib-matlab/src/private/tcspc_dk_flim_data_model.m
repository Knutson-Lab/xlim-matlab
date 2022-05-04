function [I,Jalpha,Jtau] = tcspc_dk_flim_data_model...
    (tcal,ndata,irf,param,ntau,taur,convflag)
    %TCSPC_DK_FLIM_DATA_MODEL Model for fitting basic FLIM curves
    %% Convolve parameters with IRF
    if ~isempty(taur)
        [D,~,Jtau] = tcspc_dk_delta_function_conv...
            (tcal,ndata,irf,param,ntau,taur,convflag);
    else
        switch convflag
            case "FFT"
                [D,~,Jtau] = tcspc_dk_conv_fft(tcal,ndata,irf,param,ntau,0);
        end
    end
    %% Weight by amplitude appropriately
    alpha = param(1:2:(2 * ntau));
    I = sum(alpha' .* D,2);
    Jalpha = D;
    Jtau = alpha' .* Jtau;
end

