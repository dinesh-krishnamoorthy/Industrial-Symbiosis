function [GOR_real,Constraints] = Disturbances(sim_k,nIter,par)

    GOR_real = par.GOR;
    Constraints.w_in = 9.5;
    Constraints.w_pg = 27.5;
    
%% GOR - 8 changes
    if sim_k >= nIter*3/24 && sim_k < nIter*6/24
        GOR_real = [par.GOR(1)+par.GOR_var(1);par.GOR(2)+par.GOR_var(2); ....
                    par.GOR(3); par.GOR(4);par.GOR(5);par.GOR(6)];    
    elseif sim_k>= nIter*6/24 && sim_k < nIter*9/24
        GOR_real = [par.GOR(1);par.GOR(2); ....
                    par.GOR(3)+par.GOR_var(3); par.GOR(4)-par.GOR_var(4);par.GOR(5);par.GOR(6)];    
    elseif sim_k>= nIter*9/24 && sim_k < nIter*11/24
        GOR_real = [par.GOR(1)-par.GOR_var(1);par.GOR(2)-par.GOR_var(2); ....
                    par.GOR(3); par.GOR(4);par.GOR(5)+par.GOR_var(5);par.GOR(6)+par.GOR_var(6)];    
    elseif sim_k>= nIter*11/24 && sim_k < nIter*13/24
        GOR_real = [par.GOR(1)+par.GOR_var(1);par.GOR(2)+par.GOR_var(2); ....
                    par.GOR(3)+ par.GOR_var(3)*((sim_k-nIter*11/24)/(nIter*13/24-nIter*11/24));...
                    par.GOR(4)- par.GOR_var(4)*((sim_k-nIter*11/24)/(nIter*13/24-nIter*11/24));...
                    par.GOR(5)-par.GOR_var(5);par.GOR(6)];    
    elseif sim_k>= nIter*13/24 && sim_k < nIter*16/24
            GOR_real = [par.GOR(1);par.GOR(2); ....
                        par.GOR(3); par.GOR(4);par.GOR(5);par.GOR(6)];    
    elseif sim_k>= nIter*16/24 && sim_k < nIter*20/24
    GOR_real = [par.GOR(1)+par.GOR_var(1)*((sim_k-nIter*16/24)/(nIter*20/24-nIter*16/24)); ...
                par.GOR(2)+par.GOR_var(2)*((sim_k-nIter*16/24)/(nIter*20/24-nIter*16/24)); ...
                par.GOR(3)-par.GOR_var(3)*((sim_k-nIter*16/24)/(nIter*20/24-nIter*16/24)); ...
                par.GOR(4)-par.GOR_var(4)*((sim_k-nIter*16/24)/(nIter*20/24-nIter*16/24)); ...
                par.GOR(5)-par.GOR_var(5)*((sim_k-nIter*16/24)/(nIter*20/24-nIter*16/24)); ...
                par.GOR(6)+par.GOR_var(6)*((sim_k-nIter*16/24)/(nIter*20/24-nIter*16/24))];    
    elseif sim_k>= nIter*20/24 && sim_k < nIter*22/24
        GOR_real = [par.GOR(1)+par.GOR_var(1);par.GOR(2)+par.GOR_var(2); ....
                    par.GOR(3)+par.GOR_var(3); par.GOR(4)-par.GOR_var(4);...
                    par.GOR(5)-par.GOR_var(5);par.GOR(6)-par.GOR_var(6)];    
    elseif sim_k>= nIter*22/24
        GOR_real = [par.GOR(1)+par.GOR_var(1);par.GOR(2)-2*par.GOR_var(2); ....
                    par.GOR(3)+par.GOR_var(3); par.GOR(4)+2*par.GOR_var(4);...
                    par.GOR(5)+par.GOR_var(5);par.GOR(6)+par.GOR_var(6)];    
    end         
      
        %% GOR - 3 changes

        %     if sim_k >= nIter/2 && sim_k < nIter*2/3
        %     GOR_real = [par.GOR(1)+par.GOR_var(1);par.GOR(2)+par.GOR_var(2); ....
        %                 par.GOR(3); par.GOR(4);par.GOR(5);par.GOR(6)];    
        % 
        %         else if sim_k >= nIter*2/3 && sim_k < nIter*5/6                 
        %                     GOR_real = [par.GOR(1);par.GOR(2); par.GOR(3) + par.GOR_var(3)*((sim_k-nIter*2/3)/(nIter*5/6-nIter*2/3)); ...
        %                 par.GOR(4) + par.GOR_var(4)*((sim_k-nIter*2/3)/(nIter*5/6-nIter*2/3)); par.GOR(5);par.GOR(6)];
        % 
        %             else if sim_k >= nIter*5/6 && sim_k <= nIter
        %             GOR_real = [par.GOR(1);par.GOR(2);par.GOR(3); par.GOR(4); ...
        %                     par.GOR(5) + par.GOR_var(5);par.GOR(6)+ par.GOR_var(6)];
        %                 end
        %             end
        %     end
        %     
%% Gas produced constraint
    
    if sim_k >= nIter/2 && sim_k < nIter*16/24
        Constraints.w_pg = 24;
    elseif sim_k >= nIter*16/24 && sim_k < nIter*21/24
        Constraints.w_pg = 27.5;
    elseif sim_k >= nIter*21/24
        Constraints.w_pg = 24;
    end
    
%% Gas Lift constraint

    if sim_k >= nIter*6/24 && sim_k < nIter*8/24
        Constraints.w_in = 7.5;
    elseif sim_k >= nIter*8/24 && sim_k < nIter*14/24 
        Constraints.w_in = 12;
    elseif sim_k >= nIter*14/24 && sim_k
        Constraints.w_in = 9.5;
    end
        
end