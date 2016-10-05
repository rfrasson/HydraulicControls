function [ReachBoundaries,Concavity]=FindReachesHydControl(RiverObs,DamLocations,SWATHboundaries,MinReachLen,tcritReach,numbregresspts,Makeplots)
%This function defines reach boundaries based on inflection points detected in
%the water surface profile, swath boundaries, and detected 

%List of inputs
%RiverObs nodes    : Structure containing RiverObs nodes (see explanation
%                    below)
%RiverObsDams      : Structure containing RiverObs nodes (see explanation
%                    below) used for detection of dams. Best when using
%                    average of a few "low flow" overpasses
%DamLocations      : Node ids that contain dams
%SWATHboundaries   : Location of the intersections between SWATH boundaries
%                    and the river centerline as flow distances in km.
%MinReachLen       : Minimum reach length in km
%tcritReach        : Number os Standard Deviations used to create
%                    dH/dx confidense intervals. Tested value = 2
%numbregresspts    : Number of regression points used to estimate dS/dx
%                    tested value = 10


%Makeplots         : Generate reach boundary plots if Makeplots =1

%list of outputs:
%ReachBoundaries   : Indices of the nodes that correspond to reach
%                    boundaries
%Concavity         : Concavity of the WSP at computed at each node


%Clarification on the format of RiverObs nodes.
%The variable passed as RiverObs must be a structure containing the
%following elements (as vectors with dimenstion eaqual to Number of nodes,1:
%x          : flow distance measured in meters
%H          : Node elevation in meters
%W          : River width in meters
%Easting    : Projected coordinate of the node in m
%Northing   : Projected coordinate of hte node in m
%Lat        : Latitude of the node in decimal degrees
%Lon        : Longitude of the node in decimal degrees
%xtrack     : cross-track position of the node

%Author: Renato Frasson August 31, 2015

    x=RiverObs.x/1000;
    y=RiverObs.H;

    %Find closest node to the Swath boundary
    SwathBoundIDs=nan(size(SWATHboundaries));
    NodeID=1;
    for countSwath=1:length(SWATHboundaries);
        countnode=NodeID;
        while countnode<length(x)&& SWATHboundaries(countSwath)>x(countnode)
            countnode=countnode+1;
        end
        if countnode>1&&countnode<length(x)
            SwathBoundIDs(countSwath)=countnode;
        end
        NodeID=countnode;
    end
    SwathBoundIDs=SwathBoundIDs(~isnan(SwathBoundIDs));
    
    %Search for more reaches. 
    %Step 1: break up river at dams
    %Step 2: search for inflection points (DetectCurvatureLinear)
    %Step 3: Combine ReachBoundaries in a single, concise variable
    %Step 4: add boundaries at the swath's edges
    if isempty(DamLocations)
        %treat the river as one
        [ReachBoundaries,Concavity,~,~]=DetectCurvatureLinear(x,y,1,MinReachLen,numbregresspts,tcritReach);
    else
        %Step 2:
        NumberBoundaries=0;
        Concavity=nan(size(x));
        beg=1;
        endn=DamLocations(1)-1;
        countsections=1;
        [Section(countsections).ReachBoundaries,Concavity(beg:endn),~,~]=DetectCurvatureLinear(x(beg:endn),y(beg:endn),1,MinReachLen,numbregresspts,tcritReach);
        NumberBoundaries=NumberBoundaries+length(Section(countsections).ReachBoundaries);
        countsections=countsections+1;
        for count=1:length(DamLocations)-1
            beg=DamLocations(count);
            endn=DamLocations(count+1)-1;
            [Section(countsections).ReachBoundaries,Concavity(beg:endn),~,~]=DetectCurvatureLinear(x(beg:endn),y(beg:endn),1,MinReachLen,numbregresspts,tcritReach);
            offset=beg-Section(countsections).ReachBoundaries(1); %Section(count).ReachBoundaries comes out with indices starting in 1, so we need to apply an offset
            %to place them with respect to "global" indices
            Section(countsections).ReachBoundaries=Section(countsections).ReachBoundaries + offset;
            NumberBoundaries=NumberBoundaries+length(Section(countsections).ReachBoundaries);
            countsections=countsections+1;
        end
        beg=DamLocations(length(DamLocations));
        endn=length(x);
        [Section(countsections).ReachBoundaries,Concavity(beg:endn),~,~]=DetectCurvatureLinear(x(beg:endn),y(beg:endn),1,MinReachLen,numbregresspts,tcritReach);
        offset=beg-Section(countsections).ReachBoundaries(1); %Section(count).ReachBoundaries comes out with indices starting in 1, so we need to apply an offset
        %to place them with respect to "global" indices
        Section(countsections).ReachBoundaries=Section(countsections).ReachBoundaries + offset;
        NumberBoundaries=NumberBoundaries+length(Section(countsections).ReachBoundaries);     
        %Step 3:
        ReachBoundaries=nan(NumberBoundaries,1);
        endn=length(Section(1).ReachBoundaries);
        ReachBoundaries(1:endn)=Section(1).ReachBoundaries(1:endn);
        for count=2:length(Section)
            beg=endn+1;
            endn=endn+length(Section(count).ReachBoundaries);
            ReachBoundaries(beg:endn)=Section(count).ReachBoundaries; 
        end
    end
    %add Swath boundaries
    ReachBoundaries = [ReachBoundaries; SwathBoundIDs];
    ReachBoundaries = sort(ReachBoundaries);
    if Makeplots==1
        figure
        plot(RiverObs.Easting,RiverObs.Northing)
        hold on
        plot(RiverObs.Easting(ReachBoundaries),RiverObs.Northing(ReachBoundaries),...,
                'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
        legend('Centerline', 'Reach boundaries');
        xlabel('Easting (km)');
        ylabel('Northing (km)');
        
        figure
        for count=1:length(ReachBoundaries)-1
            if ReachBoundaries(count+1)-ReachBoundaries(count)>1
                %this if skips considering dams as a separate reach
                if mean(Concavity(ReachBoundaries(count):ReachBoundaries(count+1)))>0
                    %plot in blue
                    plot(x(ReachBoundaries(count):ReachBoundaries(count+1)),y(ReachBoundaries(count):ReachBoundaries(count+1)),'o',...
                           'MarkerEdgeColor','b',...
                           'MarkerFaceColor','b',...
                           'MarkerSize',6)
                       hold on
                    
                else
                    %plot in red
                    plot(x(ReachBoundaries(count):ReachBoundaries(count+1)),y(ReachBoundaries(count):ReachBoundaries(count+1)),'o',...
                           'MarkerEdgeColor','r',...
                           'MarkerFaceColor','r',...
                           'MarkerSize',6)
                    hold on
                end
            end
        end
        xlabel('Flow distance, km')
        ylabel('Water elevation, m')
    end
end

