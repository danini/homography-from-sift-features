%% The proposed 2SIFT solver
function sols = solverHomography2SIFT(data)
    C = setupEliminationTemplates(data);
    C0 = C(:,1:6);
    C1 = C(:,7:end);
    C1 = C0 \ C1;
    RR = [-C1(end-1:end,:);eye(4)];
    AM_ind = [5,1,6,2];
    AM = RR(AM_ind,:);
    [V,D] = eig(AM);
    V = V ./ (ones(size(V,1),1)*V(1,:));
    sols(1,:) = V(2,:);
    sols(2,:) = diag(D).';
end