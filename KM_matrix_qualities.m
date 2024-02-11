function KM_matrix_qualities(K,M,plot_flag)
if(plot_flag>0)
    
    disp({'Size of K = ',int2str(size(K))})
    disp({'Rank of K = ',int2str(rank(K))})

    maxKERR=max(max(abs(K-K')));
    maxMERR=max(max(abs(M-M')));
    
    disp({'maxKError=',num2str(maxKERR)})
    disp({'maxMError=',num2str(maxMERR)})
    
    if( plot_flag == 3)
        disp('K = ');K
    end
    
    if( plot_flag == 3)
        Ktrunc=round(K)
    end
    
    if( plot_flag == 3)
        disp('M = ');M
    end
    
    if(plot_flag >= 2)
        %format long
        format compact
        disp('K matrix eigenvalues'); eig_K=eig(K)
        disp('M matrix eigenvalues'); eig_M=eig(M)
    end
    
    disp({'Kww: ',matrix_pos_definite( K )})
end