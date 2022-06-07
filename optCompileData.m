for i = 1 :numel(genAll)
    tic;
    gen = genAll{i};
    try
        % load data
        load([pwd '\Data\DataGen\' gen '_Nov13.mat'],'Data')
        load([pwd '\Data\DataSTCW\' gen ' STCW.mat'],'curvPks','curvWalks','stopCond','boundCond')
    catch
        load([pwd '\Data\DataGen\' gen '_Nov13.mat'],'Data','curvPks','curvWalks','stopCond','boundCond')
    end
    save([pwd '\Data\DataGen\' gen '_March2022.mat'],'Data','curvPks','curvWalks','stopCond','boundCond')
end