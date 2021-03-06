clear
River='Sacramento'; 
%River='PoDS'; %downstream section of the Po River
%River='Po'; %upstream section of the Po River
pathtodata=['./RawData/' River '/'];
MinReachLen=5; 
tcritReach=2;
numbregresspts=10; %should be similar to the MinReachLen/Node Spacing
numbregressptsDams=25;
tcrit=1.5165; %for dams 85% CL for 23 DF
tcritslope=1.5165; %for slopes used in the detection of dams 85% CL for 23 DF

%sets default configurations for Po and Sacramento rivers

if strfind(River,'Sacramento')
    SWATHboundaries= [0 116927.0128 151838.8427]/1000; %division by 1000 to translate to km
    ReferenceDay='170'; %day used for the definition of reach boundaries and for the estimation of A0
    %Actual day user wants to compute reach-averaged values
    Day='23'; %23 is the lowest flow from the series. 86 and 170 are intermediate flow and 65 is the highest flow, 128 is high flow, 
    %available days are: 2 23 44 53 65 86 107 128 149 170
    %DaysUsedForDams={'2', '23', '44', '86', '107', '149', '170'};
    DaysUsedForDams={'2','23','44','86','107','149'};
    %DaysForReaches={'2', '23', '44', '86', '107', '149', '170'};
    DaysForReaches={'2','23','44','86','107','149'};
    filenameReference=[pathtodata River 'Day' ReferenceDay '.mat'];
    filenameDataset=[pathtodata River 'Day' Day '.mat'];
end
if strfind(River,'Po')
    ReferenceDay='220'; %day used to trace the center line and to estimate A0
    RefOverpass='560';
    SWATHboundaries=[];
    %For overpass 560: available days are: '136' '157' '178' '199' '220' '241' '262' '283' '304' '325' '346' '367' '388' '409'
    %For overpass 211: available days are: '145' '166' '187' '208' '229' '250' '271' '292' '313' '334' '355' '376' '397' '418' '439' '460' '481'
    Day='157';
    Overpass='560';
    DaysUsedForDams={'220', '241', '304', '283', '136', '346'};
    PassForDams='560';
    DaysForReaches={'220', '241','304'};
    PassForReaches='560';
    filenameReference=[pathtodata River 'Day' ReferenceDay '-' RefOverpass '.mat'];
    filenameDataset=[pathtodata River 'Day' Day '-' Overpass '.mat'];
end

Makeplots=1;
OutputPath='./output/';
OutFileName=[River 'HC' Day '.mat'];
SaveResults=0;
SmoothData=1;
VariableSmoothingWindow=1;
%builds average profile for dam finding.
for count=1:length(DaysUsedForDams)
    if strfind(River,'Sacramento')
        fnamedams=[pathtodata River 'Day' DaysUsedForDams{count} '.mat'];
    else
        if strfind(River,'Po')
           fnamedams=[pathtodata River 'Day' DaysUsedForDams{count} '-' PassForDams '.mat']; 
        end
    end
    DamNodes=load(fnamedams,'RiverObs');
    DamNodes=DamNodes.RiverObs;
    DamH(:,count)=DamNodes.H;        
end
if length(DaysUsedForDams)>1
    DamNodes.H=mean(DamH,2);    
end
Dams=DetectDamsFastCheckSlope(DamNodes.x,DamNodes.H, numbregressptsDams,tcrit,tcritslope,0);
%builds average profile for reach definition
for count=1:length(DaysForReaches)
    if strfind(River,'Sacramento')
        fnameReach=[pathtodata River 'Day' DaysForReaches{count} '.mat'];
    else
        if strfind(River,'Po')
           fnameReach=[pathtodata River 'Day' DaysForReaches{count} '-' PassForReaches '.mat'];
        end
    end
    ReachNodes=load(fnameReach,'RiverObs');
    ReachNodes=ReachNodes.RiverObs;
    ReachH(:,count)=ReachNodes.H;        
end
ReachNodes.H=mean(ReachH,2);
AveProfile=ReachNodes;
[ReachBoundaries,StructureFlag,Concavity]=FindReachesHydControl(AveProfile,Dams,SWATHboundaries,MinReachLen,tcritReach,numbregresspts,Makeplots);
load(filenameReference);
RefRiverObs=RiverObs;
RefTrue=True;
load(filenameDataset,'RiverObs','True');

% %add missing nodes:
% RiverObs.H(40:68)=nan;
% RiverObs.W(40:68)=nan;
% 
% RiverObs.H(325:330)=nan;
% RiverObs.W(325:330)=nan;

[Reach,RiverData,Metadata,ReachTrue,Nodes,NodesTrue]=ReachAveraging(ReachBoundaries, Dams, StructureFlag, RiverObs,True,RefRiverObs,RefTrue,Day,SaveResults,SmoothData,VariableSmoothingWindow,OutputPath, OutFileName,Makeplots);