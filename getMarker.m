function marker = getMarker()
%% marker = getMarker() 
%  returns a plot marker

lines   = {'-','--',':','-.'};
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
colors  = {'r','g','b','c','m','y','k','w'};

marker = markers{mod(i,numel(markers))+1};

end