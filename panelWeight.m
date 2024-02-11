function panWt = panelWeight( solDATA,lamDATA,panXY,psiI,etaI,gl_wts,plot_flag )
%% calculates panel Weight 
%
if(lamDATA.tDeg > 0)
    panWt = panelWeightVthk2( lamDATA,panXY,psiI,etaI,gl_wts,solDATA.intMethod,plot_flag );
else
    panWt = panelWeightVthk(solDATA,lamDATA,panXY );
end