function eigError=getEigErr(panXY,solDATA,lamDATA,ldsDATA,Eg,EgN,plot_flag)
%% get error in eigenvalue solution
%disp('getEigError called');

eigError=0;
if(plot_flag >= 1)
    if(solDATA.isoltype == 1)
        eigError=eigErr(panXY,solDATA,lamDATA,EgN);
    elseif(solDATA.isoltype == 2)
        eigError=eigErrBuckl(panXY,lamDATA,ldsDATA,Eg);
    end
    %disp({'Ritz EgN=',EgN'}); EgN'
end
