function [Labels OutFinFor OutFinRev Clust] = BMA_ForRevInf(someexps)

load('BrainMapMetaData_NormalActs.mat', 'Experiments');

OutFin = [];
Labels = [];
Clust = [];
OutFinFor = [];
OutFinRev= [];

finalexps = Experiments;

finalexpbmapid = [];
finalexpexpid = [];

for a = 1:numel(finalexps)
    finalexpbmapid = [finalexpbmapid; finalexps(a).Citation.BrainMap_ID];
    finalexpexpid = [finalexpexpid; finalexps(a).Citation.Experiment_ID];
end

BD = 'dummy';
PC = 'dummy';
SM = 'dummy';
ST = 'dummy';
RM = 'dummy';
RT = 'dummy';
Instruct = 'dummy';
Contrast = 'dummy';

for a = 1:numel(finalexps)
    temp = finalexps(a).Experiments.Behavioral_Domain;
    for b = 1:numel(temp)
        BD = [BD;temp(b)];
    end
    clear temp
    temp = finalexps(a).Experiments.Paradigm_Class;
    for b = 1:numel(temp)
        PC = [PC;temp(b)];
    end
    clear temp
    temp = finalexps(a).Conditions.Stimulus_Modality;
    for b = 1:numel(temp)
        SM = [SM;temp(b)];
    end
    clear temp
    temp = finalexps(a).Conditions.Stimulus_Type;
    for b = 1:numel(temp)
        ST = [ST;temp(b)];
    end
    clear temp
    temp = finalexps(a).Conditions.Response_Modality;
    for b = 1:numel(temp)
        RM = [RM;temp(b)];
    end
    clear temp
    temp = finalexps(a).Conditions.Response_Type;
    for b = 1:numel(temp)
        RT = [RT;temp(b)];
    end
    clear temp
    temp = finalexps(a).Conditions.Instructions;
    for b = 1:numel(temp)
        Instruct = [Instruct;temp(b)];
    end
    clear temp
    temp = finalexps(a).Experiments.Contrast;
    for b = 1:numel(temp)
        Contrast = [Contrast;temp(b)];
    end
    clear temp
end


BD(1,:) = [];
uniqueBD = unique(BD);
PC(1,:) = [];
uniquePC = unique(PC);
SM(1,:) = [];
uniqueSM = unique(SM);
ST(1,:) = [];
uniqueST = unique(ST);
RM(1,:) = [];
uniqueRM = unique(RM);
RT(1,:) = [];
uniqueRT = unique(RT);
Instruct(1,:) = [];
uniqueInstruct = unique(Instruct);
Contrast(1,:) = [];
uniqueContrast = unique(Contrast);

clear BD PC SM ST RM RT Instruct Contrast

for a = 1:numel(uniqueBD)
    for b = 1:numel(finalexps)
        temp = finalexps(b).Experiments.Behavioral_Domain;
        for c = 1:numel(temp)
            checkBD(c) = strcmp(temp(c), uniqueBD(a));
        end
        if length(find(checkBD)) > 0
            isBD(b) = 1;
        else
            isBD(b) = 0;
        end
        clear checkBD temp
    end
    AllBD{a} = finalexps(find(isBD));
    BDExps(a,:) = isBD;
    clear isBD
end

for a = 1:numel(uniquePC)
    for b = 1:numel(finalexps)
        temp = finalexps(b).Experiments.Paradigm_Class;
        for c = 1:numel(temp)
            checkPC(c) = strcmp(temp(c), uniquePC(a));
        end
        if length(find(checkPC)) > 0
            isPC(b) = 1;
        else
            isPC(b) = 0;
        end
        clear checkPC temp
    end
    AllPC{a} = finalexps(find(isPC));
    PCExps(a,:) = isPC;
    clear isPC
end

for a = 1:numel(uniqueSM)
    for b = 1:numel(finalexps)
        temp = finalexps(b).Conditions.Stimulus_Modality;
        for c = 1:numel(temp)
            checkSM(c) = strcmp(temp(c), uniqueSM(a));
        end
        if length(find(checkSM)) > 0
            isSM(b) = 1;
        else
            isSM(b) = 0;
        end
        clear checkSM temp
    end
    AllSM{a} = finalexps(find(isSM));
    SMExps(a,:) = isSM;
    clear isSM
end

for a = 1:numel(uniqueST)
    for b = 1:numel(finalexps)
        temp = finalexps(b).Conditions.Stimulus_Type;
        for c = 1:numel(temp)
            checkST(c) = strcmp(temp(c), uniqueST(a));
        end
        if length(find(checkST)) > 0
            isST(b) = 1;
        else
            isST(b) = 0;
        end
        clear checkST temp
    end
    AllST{a} = finalexps(find(isST));
    STExps(a,:) = isST;
    clear isST
end

for a = 1:numel(uniqueRM)
    for b = 1:numel(finalexps)
        temp = finalexps(b).Conditions.Response_Modality;
        for c = 1:numel(temp)
            checkRM(c) = strcmp(temp(c), uniqueRM(a));
        end
        if length(find(checkRM)) > 0
            isRM(b) = 1;
        else
            isRM(b) = 0;
        end
        clear checkRM temp
    end
    AllRM{a} = finalexps(find(isRM));
    RMExps(a,:) = isRM;
    clear isRM
end

for a = 1:numel(uniqueRT)
    for b = 1:numel(finalexps)
        temp = finalexps(b).Conditions.Response_Type;
        for c = 1:numel(temp)
            checkRT(c) = strcmp(temp(c), uniqueRT(a));
        end
        if length(find(checkRT)) > 0
            isRT(b) = 1;
        else
            isRT(b) = 0;
        end
        clear checkRT temp
    end
    AllRT{a} = finalexps(find(isRT));
    RTExps(a,:) = isRT;
    clear isRT
end

for a = 1:numel(uniqueInstruct)
    for b = 1:numel(finalexps)
        temp = finalexps(b).Conditions.Instructions;
        for c = 1:numel(temp)
            checkInstruct(c) = strcmp(temp(c), uniqueInstruct(a));
        end
        if length(find(checkInstruct)) > 0
            isInstruct(b) = 1;
        else
            isInstruct(b) = 0;
        end
        clear checkInstruct temp
    end
    AllInstruct{a} = finalexps(find(isInstruct));
    InstructExps(a,:) = isInstruct;
    clear isInstruct
end

for a = 1:numel(uniqueContrast)
    for b = 1:numel(finalexps)
        temp = finalexps(b).Experiments.Contrast;
        for c = 1:numel(temp)
            checkContrast(c) = strcmp(temp(c), uniqueContrast(a));
        end
        if length(find(checkContrast)) > 0
            isContrast(b) = 1;
        else
            isContrast(b) = 0;
        end
        clear checkContrast temp
    end
    AllContrast{a} = finalexps(find(isContrast));
    ContrastExps(a,:) = isContrast;
    clear isContrast
end

%The number of experiments in each Behavioral Domain
for a = 1:numel(AllBD)
    AllBDAvail(a) = numel(AllBD{a});
end

%The number of experiments in each Paradigm Class
for a = 1:numel(AllPC)
    AllPCAvail(a) = numel(AllPC{a});
end
% 
%The number of experiments in each Stimulus Modality
for a = 1:numel(AllSM)
    AllSMAvail(a) = numel(AllSM{a});
end

%The number of experiments in each Stimulus Type
for a = 1:numel(AllST)
    AllSTAvail(a) = numel(AllST{a});
end

%The number of experiments in each Response Modality
for a = 1:numel(AllRM)
    AllRMAvail(a) = numel(AllRM{a});
end

%The number of experiments in each Response Type
for a = 1:numel(AllRT)
    AllRTAvail(a) = numel(AllRT{a});
end

%The number of experiments in each Instruction
for a = 1:numel(AllInstruct)
    AllInstructAvail(a) = numel(AllInstruct{a});
end

%The number of experiments in each Contrast
for a = 1:numel(AllContrast)
    AllContrastAvail(a) = numel(AllContrast{a});
end
% 
allcoords = 0;

%Find Total Number of Coordinates in Database
for a = 1:numel(finalexps)
    allcoords = allcoords+size(finalexps(a).XYZ_Tal,1);
end


    
    
    
    
    
for count2 = 1:numel(someexps)
    count2
    tempexps = someexps{count2};

    tempexpbmapid = [];
    tempexpexpid = [];

    for a = 1:numel(tempexps)
        tempexpbmapid = [tempexpbmapid; tempexps(a).Citation.BrainMap_ID];
        tempexpexpid = [tempexpexpid; tempexps(a).Citation.Experiment_ID];
    end

    expcheck = ismember([finalexpbmapid finalexpexpid], [tempexpbmapid tempexpexpid], 'rows');

    clear check newexps

    

    for count = 1:2 % Change to 1:8 for all Metadata fields
        switch count
            case 1
                AllTask = AllBD;
                TaskAvail = AllBDAvail;
                TaskExps = BDExps;
                uniqueTask = uniqueBD;
            case 2
                AllTask = AllPC;
                TaskAvail = AllPCAvail;
                TaskExps = PCExps;
                uniqueTask = uniquePC;
            case 3
                AllTask = AllSM;
                TaskAvail = AllSMAvail;
                TaskExps = SMExps;
                uniqueTask = uniqueSM;
            case 4
                AllTask = AllST;
                TaskAvail = AllSTAvail;
                TaskExps = STExps;
                uniqueTask = uniqueST;
            case 5
                AllTask = AllRM;
                TaskAvail = AllRMAvail;
                TaskExps = RMExps; 
                uniqueTask = uniqueRM;
            case 6
                AllTask = AllRT;
                TaskAvail = AllRTAvail;
                TaskExps = RTExps;
                uniqueTask = uniqueRT;
            case 7
                AllTask = AllInstruct;
                TaskAvail = AllInstructAvail;
                TaskExps = InstructExps;
                uniqueTask = uniqueInstruct;
            case 8
                AllTask = AllContrast;
                TaskAvail = AllContrastAvail;
                TaskExps = ContrastExps;   
                uniqueTask = uniqueContrast;
        end

        for a = 1:numel(AllTask)
            PTask = numel(AllTask{a})/sum(TaskAvail);
            checkexps = AllTask{a};
            sumcoords = 0;

            checkexpbmapid = [];
            checkexpexpid = [];

            for b = 1:numel(checkexps)
                checkexpbmapid = [checkexpbmapid; checkexps(b).Citation.BrainMap_ID];
                checkexpexpid = [checkexpexpid; checkexps(b).Citation.Experiment_ID];
                sumcoords = sumcoords+size(checkexps(b).XYZ_MNI,1);
            end

            thistask = ismember([checkexpbmapid checkexpexpid], [tempexpbmapid tempexpexpid], 'rows');

            clear checkexpbmapid checkexpexpid

            PActITask = sum(thistask)/sumcoords;
            PAct = numel(tempexps)/allcoords;

            ForInf(a,1) = PActITask;
            ForInf(a,2) = sum(thistask);
            ForInf(a,3) = sumcoords*numel(tempexps)/allcoords;

            if sum(thistask)>4 & ForInf(a,2)>ForInf(a,3)%sum(thistask) should be > 4
                [table, chi2, taskp(a,1)] = crosstab(expcheck, TaskExps(a,:));
                taskp(a,2) = 1-binocdf(sum(thistask), sumcoords, numel(tempexps)/allcoords);
            else
                taskp(a,1:2) = [1 1];
            end

            clear thistask
            RevInf(a) = (PActITask*PTask)/PAct;
            clear PActITask PTask PAct checkexps table chi2
        end

        RevInf = RevInf/sum(RevInf);

        tempRevInf = zeros(size(RevInf,1), size(RevInf,2));
        tempForInf = zeros(size(ForInf,1), 1);

        PAct = numel(tempexps)/allcoords;

        for run = 1:2

            if run == 1
                xQ = taskp(:,2) < 0.05;%0.05
            elseif run == 2
                Ps = sort(taskp(:,2), 'ascend');
                S = length(Ps);
                Fi = (1:S)'/S*0.05/1;%0.05
                I = find(Ps<=Fi, 1, 'last');
                if isempty(I)
                    u = 0;
                else
                    u = Ps(I);
                end
                xQ = taskp(:,2) <= max((0.05/sum(ForInf(:,1)>0)), u);%0.05
            end

            if run == 2
                tempForInf(xQ) = ForInf(xQ,1)./PAct;
            end

            clear xQ Data Names Ps S Fi I u

            if run == 1
                xQ = taskp(:,1) < 0.05;%0.05        
            elseif run == 2
                Ps = sort(taskp(:,1), 'ascend');
                S = length(Ps);
                Fi = (1:S)'/S*0.05/1;%0.05
                I = find(Ps<=Fi, 1, 'last');
                if isempty(I)
                    u = 0;
                else
                    u = Ps(I);
                end
                xQ = taskp(:,1) <= max((0.05/sum(ForInf(:,1)>0)), u);%0.05
            end

            if run == 2
                tempRevInf(xQ) = RevInf(xQ);
            end

            clear xQ Data Names Ps S Fi I u
        end

        clear AllTask TaskAvail TaskExps taskp a ans b c finalexps run sumcoords temptask

        FinForInf(:,1) = ForInf(:,1)./PAct;
        FinRevInf(:,1) = RevInf(1,:);
        FinForInfsig(:,1) = tempForInf(:,1);
        FinRevInfsig(:,1) = tempRevInf(1,:);
        OutFinForInf(:,1) = dataset(uniqueTask);
        OutFinRevInf(:,1) = dataset(uniqueTask);
        OutFinForInf_sig(:,1) = dataset(uniqueTask);
        OutFinRevInf_sig(:,1) = dataset(uniqueTask);
        OutFinForInf(:,2) = dataset(FinForInf);
        OutFinRevInf(:,2) = dataset(FinRevInf);
        OutFinForInf_sig(:,2) = dataset(FinForInfsig);
        OutFinRevInf_sig(:,2) = dataset(FinRevInfsig);
        
        clear PAct

        %OutFin = [OutFin; FinForInfsig; FinRevInfsig];
        OutFinFor = [OutFinFor; FinForInf];
        OutFinRev = [OutFinRev; FinRevInf];
        %Labels = [Labels; uniqueTask; uniqueTask];
        Labels = [Labels; uniqueTask];
        %Clust = [Clust; ones(numel(uniqueTask),1)*count2; ones(numel(uniqueTask),1)*count2];
        Clust = [Clust; ones(numel(uniqueTask),1)*count2];

    %     switch count
    %         case 1
    %             mkdir('BehavioralDomain')
    %             cd('BehavioralDomain')
    %         case 2
    %             mkdir('ParadigmClass')
    %             cd('ParadigmClass')
    %         case 3
    %             mkdir('StimulusModality')
    %             cd('StimulusModality')
    %         case 4
    %             mkdir('StimulusType')
    %             cd('StimulusType')
    %         case 5
    %             mkdir('ResponseModality')
    %             cd('ResponseModality')
    %         case 6
    %             mkdir('ResponseType')
    %             cd('ResponseType')
    %         case 7
    %             mkdir('Instructions')
    %             cd('Instructions')
    %         case 8
    %             mkdir('Contrast')
    %             cd('Contrast')
    %     end
    %     
    %     export(OutFinForInf)
    %     export(OutFinRevInf)
    %     export(OutFinForInf_sig)
    %     export(OutFinRevInf_sig)
    %             
    %     cd ..

        clear FinForInf FinRevInf FinForInfsig FinRevInfsig OutFinForInf OutFinRevInf OutFinRevInf_sig ForInf RevInf tempForInf tempRevInf uniqueTask OutFinForInf_sig
    end
    
    clear tempexps
        
end

    %Labels(find(OutFin==0)) = [];
    %Clust(find(OutFin==0)) = [];
    %OutFin(find(OutFin==0)) = [];
    
    %SigMetaData = dataset(Labels, OutFin, Clust);