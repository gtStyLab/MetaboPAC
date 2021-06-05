function penalty = pen_twoContMetCorr_chass(concMatrix,fluxMatrix,timeVec,fluxTimeVec)
% Calculate correlation penalty between metabolite and flux when holding
% one metabolite constant.

% First two columns are controller metabolites, second column is the target
% flux that the controller metabolites interact with.
twoContMet_interactions = {[2 18],4;
                           [3 10],6;
                           [5 7],12;
                           [5 6],13;
                           [7 8],16;
                           [8 9],18;
                           [9 10],19;
                           [4 10],20;
                           [4 10],21;
                           [10 17],24;
                           [13 16],28;
                           [13 14],29;
                           [4 18],31};

for interaction_iter = 1:length(twoContMet_interactions)
    met1 = twoContMet_interactions{interaction_iter,1}(1);
    met2 = twoContMet_interactions{interaction_iter,1}(2);
    flux = twoContMet_interactions{interaction_iter,2};
    
    x1 = concMatrix(:,met1);
    x2 = concMatrix(:,met2);
    v = fluxMatrix(:,flux);

    new_timeVec = linspace(timeVec(1),timeVec(end),201)';
    new_fluxTimeVec = linspace(fluxTimeVec(1),fluxTimeVec(end),200)';

    new_x1 = interp1(timeVec,x1,new_timeVec(1:end-1));
    new_x2 = interp1(timeVec,x2,new_timeVec(1:end-1));
    new_x1_x2 = [new_x1 new_x2];
    new_v = interp1(fluxTimeVec,v,new_fluxTimeVec);

    num_intersect_lines = 25;
    y1 = linspace(min(new_x1),max(new_x1),num_intersect_lines)';
    y2 = linspace(min(new_x2),max(new_x2),num_intersect_lines)';
    y_v = linspace(min(new_v),max(new_v),num_intersect_lines)';

    corr_x1 = [];
    corr_x2 = [];
    for iter = 2:num_intersect_lines-1
        P1 = InterX([new_x1'; new_x2'],[y1(iter)*ones(size(y2))'; y2']);
        if size(P1,2) > 1
            intersect_y2 = P1(2,:);
            F1 = scatteredInterpolant(new_x1,new_x2,new_v);
            intersect_y2_flux = F1(P1(1,:),P1(2,:));
            corr_x1 = [corr_x1; corr(intersect_y2',intersect_y2_flux','Type','Spearman')];
        end
        P2 = InterX([new_x1'; new_x2'],[y1'; y2(iter)*ones(size(y2))']);
        if size(P2,2) > 1
            intersect_y1 = P2(1,:);
            F2 = scatteredInterpolant(new_x1,new_x2,new_v);
            intersect_y1_flux = F2(P2(1,:),P2(2,:));
            corr_x2 = [corr_x2; corr(intersect_y1',intersect_y1_flux','Type','Spearman')];
        end
    end

    if flux == 4
        penalty(1) = abs(nanmean(corr_x1) + 1) + abs(nanmean(corr_x2) - 1);
    elseif flux == 6
        penalty(2) = abs(nanmean(corr_x1) + 1) + abs(nanmean(corr_x2) - 1);
    elseif flux == 12
        penalty(3) = abs(nanmean(corr_x1) + 1) + abs(nanmean(corr_x2) - 1);
    elseif flux == 13
        penalty(4) = abs(nanmean(corr_x1) - 1) + abs(nanmean(corr_x2) + 1);
    elseif flux == 16
        penalty(5) = abs(nanmean(corr_x1) + 1) + abs(nanmean(corr_x2) - 1);
    elseif flux == 18
        penalty(6) = abs(nanmean(corr_x1) + 1) + abs(nanmean(corr_x2) - 1);
    elseif flux == 19
        penalty(7) = abs(nanmean(corr_x1) + 1) + abs(nanmean(corr_x2) - 1);
    elseif flux == 20
        penalty(8) = abs(nanmean(corr_x1) - 1) + abs(nanmean(corr_x2) - 1);
    elseif flux == 21
        penalty(9) = abs(nanmean(corr_x1) - 1) + abs(nanmean(corr_x2) - 1);
    elseif flux == 24
        %penalty(10) = abs(nanmean(corr_x1) - 1) + abs(nanmean(corr_x2) + 1);
        penalty(10) = 0;
    elseif flux == 28
        penalty(11) = abs(nanmean(corr_x1) + 1) + abs(nanmean(corr_x2) - 1);
    elseif flux == 29
        penalty(12) = abs(nanmean(corr_x1) + 1) + abs(nanmean(corr_x2) - 1);
    elseif flux == 31
        penalty(13) = abs(nanmean(corr_x1) - 1) + abs(nanmean(corr_x2) - 1);
    end
end

penalty = nansum(penalty);