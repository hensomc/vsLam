function [e1,e2,nu12,g12,rho,mtype] = get_mat_props(mat_names)
%get_mat_props(mat_names);
%    mat_names = string array of materials
%    e1   = array of E11
%    e2   = array of E22
%    nu12 = array of nu12
%    g12  = array of G11

num_mat = length(mat_names);
for(i=1:num_mat)
    
    %switch mat_names(i)    
    if(strcmp(mat_names(i),'carbon_ep'))
        %case 'carbon_ep'  % Carbon/Ep Tape
%             e1(i)  = 30.0e6;
%             e2(i)  = 1.5e6;
%             g12(i) = 0.91e6;
%             nu12(i)= 0.3;
%             rho(i) = 0.058;
            e1(i)  = 22.15e6;
            e2(i)  = 1.38e6;
            g12(i) = 0.86e6;
            nu12(i)= 0.321;
            rho(i) = 0.058;
            mtype(i) = {'ortho'};
            
    elseif(strcmp(mat_names(i),'aluminum'))
        %case 'aluminum' % Aluminum           
            e1(i)  = 10.0e6;
            e2(i)  = 10.0e6;
            nu12(i)= 0.3;
            g12(i) = e1(i)/(2.*(1+nu12(i)));
            rho(i) = 0.1;
            mtype(i) = {'iso'};
          
    elseif(strcmp(mat_names(i),'steel'))
       %case 'steel' % Steel           
            e1(i)  = 30.0e6;
            e2(i)  = 30.0e6;
            nu12(i)= 0.3;
            g12(i) = e1(i)/(2.*(1+nu12(i)));
            rho(i) = 0.3;
            mtype(i) = {'iso'};
            
    elseif(strcmp(mat_names(i),'hi_ortho_carbon_ep'))
        %case 'high orthotropy carbon_ep'  % Carbon/Ep Tape
%             e1(i)  = 30.0e6;
%             e2(i)  = 1.5e6;
%             g12(i) = 0.91e6;
%             nu12(i)= 0.3;
%             rho(i) = 0.058;
            e1(i)  = 40.0e6;
            e2(i)  = 1.00e6;
            g12(i) = 0.50e6;
            nu12(i)= 0.25;
            rho(i) = 0.058;
            mtype(i) = {'ortho'};
            
    elseif(strcmp(mat_names(i),'syncore'))  % syncore_9872: ref: http://www.loctite.co.th/tht/content_data/LT3729_TT_Aerospace_Syncore_Design_Guide.pdf
            e1(i)  = 0.39e6;
            e2(i)  = 0.39e6;
            g12(i) = 0.15e6;
            nu12(i)= 0.30;
            rho(i) = 0.0243;
            mtype(i) = {'iso'};
    end
end
end