function [ ReachBoundaries ] = MergeShortReaches(ReachBoundaries,FlowDist,Concavity,MinReachLen)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    ReachLength=zeros(length(ReachBoundaries)-1,1);
    for count=1:length(ReachBoundaries)-1;
        ReachLength(count)=FlowDist(ReachBoundaries(count+1))-FlowDist(ReachBoundaries(count));
    end
    [MinLength,ReachID ]=min(ReachLength);
    while MinLength(1)<MinReachLen
        %merge reaches until all are larger than the minimum length
        %find the smaller reach. merge it with the surroundings until it is
        %larger than the minimum length. Once that reach is large enough, than
        %sort again and proceed with merging    
        if ReachID>1
        %the reach is not the first in the series, so check if its
        %sinuisity is closer to the up or downstream's
            if ReachID+2<=length(ReachBoundaries)
                %then there is an upstream reach, otherwhise, ReachID is
                %the last reach
                AveSinCurr=mean(Concavity(ReachBoundaries(ReachID):ReachBoundaries(ReachID+1)));
                AveSdownstr=mean(Concavity(ReachBoundaries(ReachID+1):ReachBoundaries(ReachID+2)));
                AveSupstr=mean(Concavity(ReachBoundaries(ReachID-1):ReachBoundaries(ReachID)));
                if abs(AveSdownstr-AveSinCurr)<abs(AveSupstr-AveSinCurr)
                    %current is more similar to the downstream
                    ReachBoundaries=[ReachBoundaries(1:ReachID);ReachBoundaries(ReachID+2:length(ReachBoundaries))];
                    ReachLength(ReachID)=FlowDist(ReachBoundaries(ReachID+1))-FlowDist(ReachBoundaries(ReachID));
                    ReachLength=[ReachLength(1:ReachID);ReachLength(ReachID+2:length(ReachLength))];
                else
                    %current is more similar to the upstream
                    ReachBoundaries=[ReachBoundaries(1:ReachID-1);ReachBoundaries(ReachID+1:length(ReachBoundaries))];
                    ReachLength=[ReachLength(1:ReachID-1);ReachLength(ReachID+1:length(ReachLength))];
                    ReachID=ReachID-1;
                    ReachLength(ReachID)=FlowDist(ReachBoundaries(ReachID+1))-FlowDist(ReachBoundaries(ReachID));                                                            
                end
            else
                %ReachID is the last reach, so merge with upstream
                ReachBoundaries=[ReachBoundaries(1:ReachID-1);ReachBoundaries(ReachID+1)];
                ReachID=ReachID-1;
                ReachLength=ReachLength(1:ReachID);
                ReachLength(ReachID)=FlowDist(ReachBoundaries(ReachID+1))-FlowDist(ReachBoundaries(ReachID));
            end
        else
            %Since it is the first reach, merge first and second reaches
            ReachBoundaries=[ReachBoundaries(1);ReachBoundaries(3:length(ReachBoundaries))];
            ReachLength=ReachLength(2:length(ReachLength));
            ReachLength(1)=FlowDist(ReachBoundaries(2))-FlowDist(ReachBoundaries(1));
        end
        [MinLength,ReachID ]=min(ReachLength);
    end


end

