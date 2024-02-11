function bc = findAdjPanEdge(panID, panDATA)
% returns array of panel edge adjacencies

adj=[];
bc='';

for i=1:4
    if( i < 4 )
        edge1=sortrows([panDATA(panID).panXY(i,1) panDATA(panID).panXY(i,2)
                  panDATA(panID).panXY(i+1,1) panDATA(panID).panXY(i+1,2)]);
    else
        edge1=sortrows([panDATA(panID).panXY(i,1) panDATA(panID).panXY(i,2)
                  panDATA(panID).panXY(1,1) panDATA(panID).panXY(1,2)]);
    end
    %disp(['Main PANEL edge ' num2str(i)]);edge1
    for j=1: length(panDATA)
        if( j ~= panID )
            %disp('adj panel ID = '); j
            for k=1:4
                if( k < 4 )
                    edge2=sortrows([panDATA(j).panXY(k,1) panDATA(j).panXY(k,2)
                          panDATA(j).panXY(k+1,1) panDATA(j).panXY(k+1,2)]);
                else
                    edge2=sortrows([panDATA(j).panXY(k,1) panDATA(j).panXY(k,2)
                          panDATA(j).panXY(1,1) panDATA(j).panXY(1,2)]);
                end

                %disp('edge2=');edge2
                if( isequal(edge1,edge2) >= 1 )
                    adj=[adj; j k];
                    break;
                end
            end
        end
    end
    if(length(adj) > 0)
        if(i<4)
            bc=strcat(bc,'c');
        else
            bc=strcat('c',bc);
        end
        adj=[];
    else
        if(i<4)
            bc=strcat(bc,'f');
        else
            bc=strcat('f',bc);
        end
    end
end

