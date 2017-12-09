function BMA_CMBHC_ClusterNum
    addpath export_fig
    theoutfile = Select_BMA(1, '.mat');
    load(theoutfile{1})
    clear theoutfile

    % Gives you options to select cluster number
    fig = dendrogram(thelink, 0);
    
    temp = thelink(:,3);
    temp(1:(length(temp)-14)) = [];
    temp(1:length(temp)) = temp(length(temp):-1:1);
    dc = (temp(1:length(temp)-1)-temp(2:length(temp)))./temp(1:length(temp)-1);
    figure, plot(3:1:length(temp)+1, dc)
    title('Relative Difference in Cophenetic Distances with Increasing Cluster Number')
    xlabel('Number of Clusters Transitioning To')
    ylabel('Difference in Cophenetic Distances')
    clear temp
    
    for a = 2:15
        t = cluster(thelink, 'MaxClust', a);
        for b = 1:a
            check{b} = find(t==b);
        end
        clear t
        fin{a-1} = check;
        clear check
    end
                
    for a = 1:13
        first = fin{a};
        second = fin{a+1};
        for b = 1:numel(first)
            for c = 1:numel(second)
                check(b,c) = isequal(first{b}, second{c});
            end
        end
        sumrow = sum(check);
        sumcol = sum(check,2);
        low = first(find(sumcol==0));
        high = second(find(sumrow==0));
        clear first second check sumcol sumrow
        low1 = numel(low{1});
        high1 = numel(high{1});
        high2 = numel(high{2});
        clear low high
        if high1>high2
            expsd(a) = high1/low1;
        else
            expsd(a) = high2/low1;
        end
        clear high1 high2 low1
    end  
    figure, plot(3:15, expsd, 'r')
    title('Experiment Separation Density with Increasing Cluster Number')
    xlabel('Number of Clusters Transitioning To')
    ylabel('Maximum Proportion of Experiment Separation')
    
    clear expsd dc fin
    
    clustnum = chooseclusters();
        
    close all
    clear fig

    if clustnum == 0
        return
    else
        t = cluster(thelink, 'MaxClust', clustnum);
    end
    
    tempdir = strcat(temptask, '_', num2str(clustnum), '_Cluster_Solution');
    mkdir(tempdir)
    
    % Begins formatting a nice black background dendrogram
        den = dendrogram(thelink, 0, 'ColorThreshold', thelink(size(thelink,1)-clustnum+2,3));
        set(gcf, 'Position', get(0,'screensize'))
        set(den, 'LineWidth', 2)
        getden = get(den);

        set(gca, 'YColor', [1 1 1]);
        set(gca, 'XColor', [1 1 1]);
        h = findobj('Color', [0 0 0]);
        gg = get(h);
        set(h, 'Color', [1 1 1]);
        set(gca, 'XColor', [0 0 0]);
        set(gca, 'Color', [0 0 0]);
        set(gcf, 'Color', [0 0 0]);

        clear h gg

        for a = 1:numel(getden)
            cols(a,:) = getden(a).Color;
        end

        uncols = unique(cols, 'rows');

        if length(find(ismember(uncols, [0 0 0], 'rows'))) > 0
            uncols(find(ismember(uncols, [0 0 0], 'rows')),:) = [];
        end

        clear cols den getden

        for a = 1:size(uncols,1)
            h = findobj('Color', uncols(a,:));
            if numel(h) < numel(finalexps)*0.05
                set(h, 'Color', [1 1 1])
            end
            clear h
        end

        clear uncols

        set(gca, 'YLim', [0 0.9])
        
        export_fig(strcat(tempdir, '/', 'Dendrogram_1'), '-tiff')

        set(gca, 'YLim', [0.85 (thelink(size(thelink,1),3)+0.1)])
        
        export_fig(strcat(tempdir, '/', 'Dendrogram_2'), '-tiff')
    
    %Exports text file with coordinates to run through GingerALE
    
    for a = 1:clustnum
        exps{a} = finalexps(find(t==a));
        
        templateline = '//Reference=Talairach';
        thefilname = strcat(tempdir, '/', temptask, '_Cluster_', num2str(a), '_of_', num2str(clustnum), '.txt');
        fid = fopen(thefilname, 'w');
        fprintf(fid, '%s\n', templateline);
        fclose(fid);
        secondstarters = exps{a};
        
        for b = 1:numel(secondstarters)
            fid = fopen(thefilname, 'a');
            subnum = secondstarters(b).Subject_Total;
            subnum = num2str(subnum);
            fprintf(fid, '%s\n', strcat('//Subjects=',subnum));
            clear subnum
            peaknum = secondstarters(b).Peaks;
            for c = 1:peaknum
                checkpts = secondstarters(b).XYZmm_Tal(c,:);
                fprintf(fid, '%f\t%f\t%f\n', [checkpts(1) checkpts(2) checkpts(3)]);
                clear checkpts
            end
            clear peaknum
            fprintf(fid, '%s\n', '');
            fclose(fid);
        end
            clear secondstarters
            
        if numel(exps) > 1
            fid = fopen(thefilname);
            a = fscanf(fid, '%s');
            a(1:2) = [];
            space = a((find(a=='=', 1, 'first')+1):(find(a=='/', 1, 'first')-1));
            clear a fid
            namefile = thefilname;
            namefile(find(namefile=='.', 1, 'last'):length(namefile)) = [];
            %GingerALE parameters
            maskTal = '-mask=Tal_wb.nii';
            maskMNI = '-mask=MNI_wb.nii';
            aleimage = strcat('-ale=', namefile, '_ALE.nii');
            pimage = strcat('-pval=', namefile, '_PVal.nii');
            threshimage = strcat('-thresh=', namefile, '_ALE_C05_1k.nii');
            pvalent = '-p=0.001';
            perment = '-perm=1000';
            clustent = '-clust=0.05';
            nonadditive = '-nonAdd';

            switch space
                case 'Talairach'
                    systemcommand = ['java -cp GingerALE.jar org.brainmap.meta.getALE2 ' thefilname ' ' maskTal ' ' perment ' ' clustent ' ' pvalent ' ' nonadditive ' ' aleimage ' ' pimage ' ' threshimage];

                case 'MNI'
                    systemcommand = ['java -cp GingerALE.jar org.brainmap.meta.getALE2 ' thefilname ' ' maskMNI ' ' pvalent ' ' perment ' ' clustent ' ' nonadditive ' ' aleimage ' ' pimage ' ' threshimage];
            end

            [status results] = system(systemcommand);
        end
            
    end
    
    [OutFin Labels Clust] = BMA_ForRevInf(exps);
    Labels(find(OutFin==0)) = [];
    Clust(find(OutFin==0)) = [];
    OutFin(find(OutFin==0)) = [];
    
    %Runs a leave-one-out metadata forward and reverse inference analysis
    %to determine if one paper is having a significant effect on the
    %metadata results
    
    for a = 1:numel(exps)
        a
        temp = exps{a};
        for b = 1:numel(temp)
            temp1(b) = temp(b).Citation.BrainMap_ID;
        end
        unbmapid = unique(temp1);
        for b = 1:length(unbmapid)
            tempexps = temp;
            check = find(unbmapid(b)==temp1);
            tempexps(check) = [];
            clear check
            expiter{b} = tempexps;
            clear tempexps
        end
        clear temp temp1 unbmapid
    
        [tempoutfin templabels tempclust] = BMA_ForRevInf(expiter);
        templabels(find(tempoutfin==0)) = [];
        tempclust(find(tempoutfin==0)) = [];
        clear expiter tempoutfin
        
        untemplabs = unique(templabels);
        
        for b = 1:numel(untemplabs)
            check(b) = isequal(unique(tempclust(find(strcmp(untemplabs{b}, templabels)))), [1:max(tempclust)]');
        end
        
        if exist('check')==1
            finlabs{a} = untemplabs(find(check));
        else
            finlabs{a} = [];
        end
        
        clear check untemplabs templabels tempclust
    end
    
    for a = 1:max(Clust)
        temp1 = find(Clust==a);
        somelabs = Labels(temp1);
        tempfinlabs = finlabs{a};
        for b = 1:numel(somelabs)
            check(b) = length(find(strcmp(somelabs(b), tempfinlabs)));
        end
        if exist('check')==1
            temp2 = find(check==0);
            Labels(temp1(temp2)) = [];
            OutFin(temp1(temp2)) = [];
            Clust(temp1(temp2)) = [];
        end
        clear temp2 check tempfinlabs somelabs temp1
    end
    
    clear finlabs
    
    SigMetaData = dataset(Labels, OutFin, Clust);
    
    save(strcat(tempdir, '/', temptask, '_', clustnum, '_Cluster_Solution.mat'), 'corrmat', 'exps', 'finalexps', 'thelink', 't', 'OutFin', 'Labels', 'Clust')
    export(SigMetaData, 'File', strcat(tempdir, '/', temptask, '_', clustnum, '_Cluster_Solution_SigMetaData.txt'))
        
    clear SigMetaData Labels Clust OutFin
    
end

function clustnum = chooseclusters()

    figure
    set(gcf, 'Color', [0 0 0.25], 'NumberTitle', 'off', 'Menubar', 'none', 'Toolbar', 'none', 'Position', [650 300 200 200])
    
    currentnum = '0';
    thistext = uicontrol('Style', 'Text', 'ForegroundColor', [0.7 0.7 0], 'String', 'Enter Number of Clusters', 'FontSize', [18], 'Position', [0 125 200 50], 'BackgroundColor', [0 0 0.25]); 
    theeditbox = uicontrol('Style', 'Edit', 'FontSize', [16], 'Position', [75 75 50 50], 'UserData', currentnum, 'Callback', @theeditbox2);
    okbox = uicontrol('Style', 'pushbutton', 'Fontsize', [12], 'Position', [25 25 50 25], 'String', 'OK', 'Callback', @okbox2);
    exitbox = uicontrol('Style', 'pushbutton', 'Fontsize', [12], 'Position', [125 25 50 25], 'String', 'Exit', 'Callback', @exitbox2);
    
    uiwait(gcf)
    clustnumscheck = get(theeditbox, 'UserData');
    clustnumscheck = str2num(clustnumscheck);
    if clustnumscheck == 0
        clustnum = 0;
    else
        clustnum = clustnumscheck;
    end
    close

end

function theeditbox2(~,~)


end

function okbox2(~,~)

    chillins = get(gcf, 'Children');
    thechoicenum = get(chillins(3), 'string');
    set(chillins(3), 'UserData', thechoicenum);
    uiresume

end