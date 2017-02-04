function [FWM_Ic, FWMfq, FWM_I] = processFWMSF(pumpS, pumpWL, SF_S, SF_WL, FWM_S, FWM_WL, plotting) 

FWM_Ic = [];
FWM_I  = [];
FWMfq  = [];


[pumpfit, pumpfitq] = fit(pumpWL, pumpS, fittype('gauss1'));
for i = 1:size(FWM_S,2);
   
    
    [FWMfit,  FWMfitq]  = fit(FWM_WL, FWM_S(:,i), fittype('gauss1'));
    [SFfit,   SFfitq]   = fit(SF_WL,  SF_S (:,i), fittype('gauss1'));
    
    if (plotting)
        subplot(1,2,1); 
        plot(FWM_WL, FWM_S(:,i), FWM_WL, FWMfit(FWM_WL));
        if (FWMfitq.adjrsquare < 0.7 || SFfitq.adjrsquare < 0.7 || pumpfitq.adjrsquare < 0.7)
            title(sprintf('FWM signal, %d  - excluded',i));
            pause;
        else
            title(sprintf('FWM signal, %d is %d', i, size(FWM_S,2)));
        end;
        subplot(1,2,2);
        plot(SF_WL, SF_S(:,i), SF_WL, SFfit(SF_WL));
        title('Sum-frequency signal');
        
        %subplot(2,2,3);
        %plot(pumpWL, pumpS(:,i), pumpWL, pumpfit(pumpWL));
        %title('Pump signal');
        
        drawnow;
    end;
    
    if (FWMfitq.adjrsquare < 0.7 || SFfitq.adjrsquare < 0.7 || pumpfitq.adjrsquare < 0.7)
        disp(sprintf('Number %d - bad fit, excluded',i));
    else
        FWM_I    = [FWM_I, FWMfit.a1]; 
        FWM_Ic   = [FWM_Ic, FWMfit.a1/SFfit.a1];
        FWMfq    = [FWMfq, (2./pumpfit.b1 - 1./SFfit.b1)./1e-7];
    end;
    
    
end;

    
    