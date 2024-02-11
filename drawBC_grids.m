function drawBC_grids(plot_flag,psiBCMesh,etaBCMesh,xBCMesh,yBCMesh)
%if(plot_flag == 1 || plot_flag == 3)
if(plot_flag >= 1)
    disp('Drawing Grids...')
    figure
    subplot(1,2,1)
    h_grid = plot(xBCMesh(:),yBCMesh(:),'r.','linewidth',3);
    axis equal; axis tight;
    text(xBCMesh(:),yBCMesh(:),int2str((1:numel(xBCMesh))'),'color','b','FontSize',12)
    title({'\bfBC Grid IDs'});
    
    subplot(1,2,2)
    h_NonGrid = plot(psiBCMesh(:),etaBCMesh(:),'r.','linewidth',3);
    axis equal; axis tight;
    text(psiBCMesh(:),etaBCMesh(:),int2str((1:numel(psiBCMesh))'),'color','b','FontSize',12)
    title({'\bfNon-Dimensional BC Grid IDs'});
end
