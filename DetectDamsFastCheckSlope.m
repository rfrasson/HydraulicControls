function [ DamIndex ] = DetectDamsFastCheckSlope( X, Y, NumbRegPoints,tcrit,tcritslope,debug)
%This function looks for dams given the water surface elevation (Y) at
%several points along the centerline.
%The location of Dams are defined as points where there is a mismatch in
%height estimated using the points downstream and upstream from the current
%location

%this version does not search for influential points to increase
%performance

%inputs:
    %Y=Water surface elevation X km downstream from the first river station
    %NumbRegPoints: Number of points used in the regression
    %tcrit=Critical T distribution values for the desired level of
    %confidence
    %debug 1 = print candidates
%outputs:
    %position in the y vector where the funciton detected a possible Dam
    order =1; %linear regression
    n=length(X);
    DamIndex=NaN(n,1);
    %debug=1; %now a parameter
    %vector containing the predictions of the WSE using the upstream(1) and donwstream (2) points
    %I have more than 1 prediction, as I generate predictions excluding one
    %point from the subset at a time, to account for uncertainty on the
    %Observed elevation
    Coef1=zeros(NumbRegPoints+1,order+1); %upsteam subset 
    Coef2=zeros(NumbRegPoints+1,order+1); %downstream subset
    WaterDrop=zeros(NumbRegPoints,1); %vector that holds discontinuities on the WSP for a subinterval
    %I will use it to decide where exactly to place a dam, once one is located
    Locations=zeros(NumbRegPoints,1);
    dropindex=1;
    NumbDams=1;
    count=NumbRegPoints+1;%index of the current examined node
    while count<n-NumbRegPoints+1 %cycles until count=n-NumbRegPoints
        
        x1=X(count-NumbRegPoints:count);
        y1=Y(count-NumbRegPoints:count);
        %select the subset of points downstream from the current node
%         x2=X(count+1:count+NumbRegPoints);
%         y2=Y(count+1:count+NumbRegPoints);
        x2=X(count:count+NumbRegPoints);
        y2=Y(count:count+NumbRegPoints);
        Coef1=polyfit(x1,y1,order); %Linear regression
        Coef2=polyfit(x2,y2,order);
        
        slope1=mean(Coef1(1));
        if slope1>0
            %we cannot resolve hydraulic jumps, so I will set slope to 0
            slope1=0;
            intercept1=mean(y1);
        else
            intercept1=mean(Coef1(2));
        end
        slope2=mean(Coef2(1));
        if slope2>0 %same thing for the downstream subinterval
            slope2=0;
            intercept2=mean(y2);
        else
            intercept2=mean(Coef2(2));
        end
        %calculate residuals
        Residual1=y1-(slope1*x1+intercept1);
        Residual2=y2-(slope2*x2+intercept2);
        s1=sqrt(1/(NumbRegPoints-2)*sum(Residual1.^2)); %regression standard error for the first line
        %coming from pag 566 of Moore's Basic practice of statistics
        s2=sqrt(1/(NumbRegPoints-2)*sum(Residual2.^2));%regression standard error for the second line
        yhat1=(slope1*X(count)+intercept1); %predicted water surface elevation at the evaluated point using curve 1
        yhat2=(slope2*X(count)+intercept2); %predicted water surface elevation at the evaluated point using curve 2
        
        %if the confidence intervals for yhat2 and yhat1 do not overlap,
        %than might be a dam. Since if there is a dam, yhat2 must be smaller than yhat1 as water cannot climb a dam
        %coefficient 2.365 comes from t distribution considering 8 degrees of freedom and 95% CI
        %tcrit=2.365; %now an input
      
        lower1=yhat1-s1*sqrt(1+1/NumbRegPoints+((X(count)-mean(x1))^2)/sum((x1-mean(x1)).^2))*tcrit; %formula from Moore's Basic practice of statistics page 579
        upper2=yhat2+s2*sqrt(1+1/NumbRegPoints+((X(count)-mean(x2))^2)/sum((x2-mean(x2)).^2))*tcrit;   
        
        %I'm assuming that the slope upstream from the dam must be milder
        %than the slope downstream, so I calculate the confidence interval
        %for slopes and see if the slope downstream is (statistically)
        %significantly milder than the slope upstream.
        
        SEb=s1/sqrt(sum((x1-mean(x1)).^2));
        %slope1up=slope1+tcritslope*SEb;
        slope1down=slope1-tcritslope*SEb;
        
        SEb=s2/sqrt(sum((x2-mean(x2)).^2));
        slope2up=slope2+tcritslope*SEb;
        %slope2down=slope2-tcritslope*SEb;        
        
        if upper2<lower1 && slope1down>slope2up %slope upstream must be less negative than slope downstream           
            %dam detected in this sub-interval
            %WaterDrop(dropindex)=yhat1-yhat2;
            WaterDrop(dropindex)=lower1-upper2;
            Locations(dropindex)=count;
            dropindex=dropindex+1;                
            if debug==1
                disp('lower bound using curve1=')
                disp(lower1)
                disp('upper bound using curve2=')
                disp(upper2)
                disp(count)
                ypred1=(slope1*x1+intercept1); %predicted water surface elevation at the evaluated point using curve 1
                ypred2=(slope2*x2+intercept2); %predicted water surface elevation at the evaluated point using curve 2
                figure
                %plot(x1,y1,'o',x2,y2,'o',x1,ypred1,'-r',x2,ypred2,'-b')
                plot(x1,y1,'o','MarkerEdgeColor','b','MarkerFaceColor','b')
                hold on
                plot(x2,y2,'o','MarkerEdgeColor','r','MarkerFaceColor','r')
                plot(x1,ypred1,'-b',x2,ypred2,'-r')
                title(['Possible dam at point: ', num2str(count)])
            end
            count=count+1;
        else           
            count=count+1;
            if dropindex~=1
                %there was an occurence of a potential dam, pinpoint location
                [~,index]=max(WaterDrop);
                DamIndex(NumbDams)=Locations(index);%location of the maximum drop in water elevation
                NumbDams=NumbDams+1;
                WaterDrop(:)=0;%reset drops
            end
            dropindex=1;
        end
    end
    DamIndex=DamIndex(~isnan(DamIndex));
end

