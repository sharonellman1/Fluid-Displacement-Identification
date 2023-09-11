%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to determine displacements and calculate event based curvatures
% and contact angles

% By Tom Bultreys, Arjen Mascini and Sharon Ellman
% Adapted for Schluter et al. data no pad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run setUpData.m
indexSafe = 1;

% create sparse connection matrix
C=sparse(linksClean(:,1), linksClean(:,2), linksClean(:,1), numPores, numPores );
C(C>0) = 1;
C = C+C';

% create vector to keep track if a pore has been filled before in the time series,
% to reduce jitter artifacts
indexHasBeenFilled = zeros(numPores, 1);
indexHasBeenFilled(indexSafe(numTimeSteps == 0)) = 1; %no jittering in pores filled at the beginning of the sequence

% remove jittering in throats for snap-off detection
indexHasBeenFilledSnapOff = zeros(numThroats, 1); 

nEvents = zeros(max(numTimeSteps),1); % events per time step
fillingVol = zeros(max(numTimeSteps),1); % total filling event volume per time step in voxels
cumulFillingVol =  zeros(max(numTimeSteps),1); % number of multi-pore events per time step
saturation =  zeros(max(numTimeSteps),1); % number of multi-pore events per time step

nMultiPoreEvents = zeros(max(numTimeSteps),1); % number of multi-pore events per time step

eventNumber = 1;
eventSizes = []; 
eventVolumesIn = [];
eventVolumesAll = [];
eventT = [];

inEventPoreNumbers = {};
inEventThroatNumbers = {};
outEventPoreNumbers = {};
outEventThroatNumbers = {};
allEventPoreNumbers = {};
allEventThroatNumbers = {};

noNeighborEvent = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avEventCurvatureIn = [];
maxEventCurvatureIn = [];
avEventStdIn = [];

avEventCurvatureAll = [];
maxEventCurvatureAll = [];
avEventStdAll = [];

avEventCaIn = [];                                 
avEventCaStdIn = [];

avEventCaAll = [];
avEventCaStdAll = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxInvasionThroatRadius = [];
minInternalThroatRadius = [];
minInvasionThroatRadius = [];

MaxEventPoreRadius = [];
MinEventPoreRadius = [];
eventShapeFactor_node = [];
eventNodeID = [];
eventNodeID1 = [];
eventAvgAspectRatio_node = [];
eventMaxAspectRatio_node = [];

eventLinkID = [];
eventShapeFactor = []; 
eventAvgAspectRatio = [];
eventMaxAspectRatio = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fillingTPerPore = NaN(numPores,1);
fillingEventPerPore = NaN(numPores,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSnapOff = [0]; % number of snap-off's per time step
AllStandAloneSnapOff = {}; % throat numbers of throats where snap-off occurred
PoreI_nType = NaN(numPores,1); % gives the n-value of the event type (e.g. for I_1 event, the value at that pore is 1) for every pore. NaN occurs when the pore was not filled or when it was an internal pore and the event type could not be determined. 
Event_I_n_type = [];
snapOffPc = {}; % same format as AllStandAloneSnapOff
SnapCurvPrs = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extracting all CAs and Kcs for every timestep per pore %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numPores = 604;
numTimeSteps = 52;
Ca_dyIMB = nan(numPores,numTimeSteps,3); % this will be the matrix with the contact angle values
display(length(inEventPoreNumbers))

%folder with the input files
prefix_CaName1 ='./data/contactAngle_data_perPore/'; 

% last bit of the name of the files
postfix_CaName1 = 'CaOutputFile.txt';


for t = 0:51 %t here refers to the numbers on the file names
    n_strPaddedBefore = sprintf( '%03d', (t));
    CaData = dlmread([prefix_CaName1 n_strPaddedBefore  '_' postfix_CaName1], '\t', 2,0); %Sharon: '\t' means horizontal tab acts as the delimiter. Start reading at row 3, column 1
    Ca_dyIMB(:,t+1,1)=CaData(:,3); % mean contact angle
    Ca_dyIMB(:,t+1,2)=CaData(:,2); % number of data points (to calculate weighted averages)
    Ca_dyIMB(:,t+1,3)=CaData(:,4); % standard deviation
    
end

Kc_dyIMB = nan(numPores,numTimeSteps,3); % this will be the matrix with the values
%folder with the input files
prefix_kcName1 = './data/SurfaceData/Kc_files_curvatureData/perPore_geodesic/'; 
% last bit of the name of the files
postfix_kcName1 = 'KcOutputFile.txt';

for t = 0:51 %t here refers to the numbers on the file names
    n_strPaddedBefore = sprintf( '%03d', (t));
    KcData = dlmread([prefix_kcName1 n_strPaddedBefore  '_' postfix_kcName1], '\t', 2,0); %Sharon: '\t' means horizontal tab acts as the delimiter. Start reading at row 3, column 1
    Kc_dyIMB(:,t+1,1)=-KcData(:,3); % mean Kc
    Kc_dyIMB(:,t+1,2)=KcData(:,2); % number of data points (to calculate weighted averages)
    Kc_dyIMB(:,t+1,3)=KcData(:,4); % standard deviation
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplingDiff = 1; %step size for loop

for t = 2:samplingDiff:numTimeSteps
    disp(t)
    
    % identifying snap-off events
    possibleSnapOff = find((eThroatFillingClean(:,t) ~= 0) & (eThroatFillingClean(:,t-1) == 0) & (eThroatSurfaceFillingClean(:,t) ~= 0) & (eThroatSurfaceFillingClean(:,t-1) == 0) )%& (eThroatfillingFractionClean(:,t)<0.0001)); %identifies all throats that changed occupancy i.t.o. balls and surfaces %SE 12/07/2023 added & (eThroatfillingFractionClean(:,t)<0.3)
    possibleSnapOff = possibleSnapOff(indexHasBeenFilledSnapOff(possibleSnapOff) ~=1); %to remove jitters 
    StandAloneSnapOff = possibleSnapOff(ePoreFilling(linksClean(possibleSnapOff,1), t) == 0 & ePoreFilling(linksClean(possibleSnapOff,2), t)==0); %only keeps throats whose neighbouring pores are still oil
    nSnapOff = [nSnapOff, length(StandAloneSnapOff)]; %stores number of stand alone snap-offs per time step
    indexHasBeenFilledSnapOff(StandAloneSnapOff) = 1; 
    if StandAloneSnapOff;
        AllStandAloneSnapOff{t} = StandAloneSnapOff; %stores throat numbers of throats where stand alone snap-off occurred
    else; 
        AllStandAloneSnapOff{t} = NaN; %if no stand-alone snap off occurs in a timestep, a NaN is assigned
    end
    
    %Calc saturation and filled Volume
    listT = find(ePoreFilling(:,t)~= 1) ; %oil filled pores at time t
    fillingVol(t) =sum(poreVolumesAvizo.*ePorefillingFraction(:,t));
    saturation(t) = 1 - fillingVol(t)/totalPoreVolume;
           
    %find events
    event = find(ePoreFilling(:,t) ~= 0 &  ePoreFilling(:,t-1) == 0); %list of pores which have been filled between t and t-1
    event = event(indexHasBeenFilled(event) ~= 1); %to remove jitters
    indexHasBeenFilled(event) = 1; %keeping track of pores having been in an event or not
    
    %create non-sparse connection matrix for events
    Conn = C(event, event);
    Conn = full(Conn);

    %cluster pore events into connected cascade events
    if event
        [nEvents(t),eventSizesT,separateEventsT] = networkComponents(Conn);
        
        n_strPaddedBefore = sprintf( '%03d', t-KcFileOffset );
        
        kcMatrix = dlmread([prefix_kcName n_strPaddedBefore  '_' postfix_kcName], '\t', 2,0); 
        CaMatrix = dlmread([prefix_CaName n_strPaddedBefore  '_' postfix_CaName], '\t', 2,0); 
        
        %Determine sizes, volumes, rates, etc
        poreEventT = separateEventsT;
        nMultiPoreEvents(t) = sum(eventSizesT>1);
        eventSizes = [eventSizes eventSizesT];    
        
        for j = 1:length(separateEventsT)
            eventT = [eventT t]; % to track events over time
            
            poreEvent{j} = event(separateEventsT{j}); %event pores
                       
            eventPoreNeighbours = [];       
            eventThroatNeighbours = [];
            eventThroats = [];
                     
            eventCurvaturesIn = -kcMatrix(poreEvent{j}, curvatureColumnValue); % mean curvature
            eventWeightsIn = kcMatrix(poreEvent{j}, weightColumnValue);
            eventStDevIn = kcMatrix(poreEvent{j}, stdColumnValue);
            
            avEventCurvatureIn = [avEventCurvatureIn nansum(eventCurvaturesIn.*eventWeightsIn)/nansum(eventWeightsIn)]; %weighted average                              
            maxEventCurvatureIn = [maxEventCurvatureIn nanmax(eventCurvaturesIn)]; 
            avEventStdIn = [avEventStdIn nanmean(eventStDevIn)] ;
            
            
            for k=1:length(poreEvent{j})
                poreNeighbourhood = poreNeighbours{poreEvent{j}(k)};
                throatNeighbourhood = throatNeighbours{poreEvent{j}(k)};
                 
                %checks if the neighbourhood is not empty
                if ~isempty(poreNeighbourhood) && ~isempty(throatNeighbourhood)
                    eventPoreNeighbours = [eventPoreNeighbours poreNeighbourhood(~ismember(poreNeighbourhood, poreEvent{j}))]; %neighbours of the event, excluding pores of the event itself
                    eventPoreNeighbours = unique(eventPoreNeighbours); %to remove duplicates
                    
                    throatIsInternalToEvent = all(ismember(links(throatNeighbourhood,:),poreEvent{j}),2); %S.Ellman: throat is internal if both pores it is connected to are in poreEvent{j}
                    eventThroatNeighbours = [eventThroatNeighbours throatNeighbourhood(~throatIsInternalToEvent)]; %neighbours of the event, excluding throats of the event itself
                    
                    if ~isempty(throatNeighbourhood(throatIsInternalToEvent))
                        eventThroats = [eventThroats throatNeighbourhood(throatIsInternalToEvent)];
                    end
                    
                    if isempty(throatNeighbourhood(throatIsInternalToEvent))
                        eventThroats = [eventThroats nan];
                    end
                 
                end
                
                %if poreneighbourhood is empty
                if isempty(poreNeighbourhood) || isempty(throatNeighbourhood)               
                    eventPoreNeighbours = [eventPoreNeighbours nan]; % neighbors and event pores combined
                    eventThroatNeighbours = [eventThroatNeighbours nan]; %neighbours of the event, excluding pores of the event itself
                    eventThroats = [eventThroats nan];                
                    noNeighborEvent = [noNeighborEvent t];
                
                end
                       
                fillingEventPerPore((poreEvent{j}(k)),1) = eventNumber;
                fillingTPerPore((poreEvent{j}(k)),1) = t;
                
                % I_n pores
                if nansum(eThroatFilling(throatNeighbours{poreEvent{j}(k)}, t-1)) ~=0; % ensures pore is not internal, in which case it is impossible to say which event type filled it
                    PoreI_nType((poreEvent{j}(k)), 1) = nansum(~eThroatFilling(throatNeighbours{poreEvent{j}(k)}, t-1)); % counts the number of throat neighbours of pore k in event j that were oil filled in timestep t-1
                end
                
            end
            
            % combines things back together
            allEventPores = unique([poreEvent{j}' eventPoreNeighbours]); 
            allEventThroats = unique([eventThroats eventThroatNeighbours]);
            
            % keep track of all event pores
            inEventPoreNumbers = [inEventPoreNumbers; poreEvent{j}'];
            inEventThroatNumbers = [inEventThroatNumbers;eventThroats];
            outEventPoreNumbers = [outEventPoreNumbers;eventPoreNeighbours];
            outEventThroatNumbers = [outEventThroatNumbers;eventThroatNeighbours];
            allEventPoreNumbers = [allEventPoreNumbers;allEventPores];
            allEventThroatNumbers = [allEventThroatNumbers;allEventThroats];
            
            
            
            % for no neighbour cases
            if ~isnan(allEventPores)
                eventCurvaturesAll = -kcMatrix(allEventPores, curvatureColumnValue);  %S.Ellman 09/04/2021: -ve added
                eventWeightsAll = kcMatrix(allEventPores, weightColumnValue);
                eventStdAll = kcMatrix(allEventPores, stdColumnValue);  

            end
            
            if isnan(allEventPores)
                eventCurvaturesAll = nan;  
                eventWeightsAll = nan;
                eventStdAll = nan;  

            end
            
            avEventCurvatureAll = [avEventCurvatureAll nansum(eventCurvaturesAll.*eventWeightsAll)/nansum(eventWeightsAll)];
            avEventStdAll = [avEventStdAll nanmean(eventStdAll)];            
            maxEventCurvatureAll = [maxEventCurvatureAll nanmax(eventCurvaturesAll)];
            
            if isnan(eventThroatNeighbours) % to account for not connected pores
                 maxInvasionThroatRadius = [maxInvasionThroatRadius nan];
                 eventShapeFactor =  [eventShapeFactor nan];                       
                 eventAvgAspectRatio = [eventAvgAspectRatio nan];
                 eventMaxAspectRatio = [eventMaxAspectRatio nan];
                 eventLinkID = [eventLinkID nan];
                 minInvasionThroatRadius = [minInvasionThroatRadius nan];            
            

            
            end
            
            % Max pore radius of event
            [eventNodeRadius, eventNodeIndex] = nanmax(nodeRadius(poreEvent{j}));
            [minEventNodeRadius, eventNodeIndex1] = nanmin(nodeRadius(poreEvent{j}));%% added 03/05/2022
            eventNode = poreEvent{j}(eventNodeIndex);
            eventNode1 = poreEvent{j}(eventNodeIndex1);%% added 03/05/2022
            if (isempty(eventNodeRadius))
                eventNodeRadius = nan;
                eventNode = nan;
            end
            
            if (isempty(minEventNodeRadius))
                minEventNodeRadius = nan;
                eventNode1 = nan;
            end
            MaxEventPoreRadius = [MaxEventPoreRadius eventNodeRadius];
            MinEventPoreRadius = [MinEventPoreRadius minEventNodeRadius]; 
            eventNodeID = [eventNodeID eventNode];
            eventNodeID1 = [eventNodeID1 eventNode1]; 
            if(isnan(eventNode))
                evShapeFactor_node = nan;
                eventAvgAspectRatio_node = [eventAvgAspectRatio_node nan];
                eventMaxAspectRatio_node = [eventMaxAspectRatio_node nan]; 
            else
                evShapeFactor_node = nodeShapeFactor(eventNode);
                eventAvgAspectRatio_node = [eventAvgAspectRatio_node nanmean(linkRadius(cell2mat(throatNeighbours(eventNode, :)))./nodeRadius(eventNode))];
                eventMaxAspectRatio_node = [eventMaxAspectRatio_node nanmax(linkRadius(cell2mat(throatNeighbours(eventNode, :)))./nodeRadius(eventNode))];
            end
            eventShapeFactor_node =  [eventShapeFactor_node evShapeFactor_node];
            

            % Event type detection
            
            Event_I_n_type = [Event_I_n_type nansum(~eThroatFilling(eventThroatNeighbours, t-1))]; % counts the number of throat neighbours of pore k in event j that were oil filled in timestep t-1
 
            
                
            if ~isnan(eventThroatNeighbours)          
                eventThroatNeighbours_saturationTimeT = eThroatFilling(eventThroatNeighbours, t);
                eventThroatNeighbours_saturationTimeTmin1 = eThroatFilling(eventThroatNeighbours, t-1);
                
                [eventLinkRadius, eventLinkIndex] = nanmax(linkRadius(eventThroatNeighbours(~eventThroatNeighbours_saturationTimeTmin1)));
                eventLink = eventThroatNeighbours(eventLinkIndex);
                if(isempty(eventLinkRadius))
                    eventLinkRadius = nan;
                    eventLink = nan;
                end
                
                maxInvasionThroatRadius = [maxInvasionThroatRadius eventLinkRadius];

                eventLinkID = [eventLinkID eventLink];
                if(isnan(eventLink))
                    evShapeFactor = nan;
                    eventAvgAspectRatio = [eventAvgAspectRatio nan];
                    eventMaxAspectRatio = [eventMaxAspectRatio nan];                    
                    minInvasionThroatRadius = [minInvasionThroatRadius nan];

                else
                    evShapeFactor = linkShapeFactor(eventLink);                
                    eventAvgAspectRatio = [eventAvgAspectRatio nanmean(linkRadius(eventLink)./nodeRadius(links(eventLink,:)))];
                    eventMaxAspectRatio = [eventMaxAspectRatio nanmax(linkRadius(eventLink)./nodeRadius(links(eventLink,:)))];               
                    minInvasionThroatRadius = [minInvasionThroatRadius min(linkRadius(eventThroatNeighbours))];

                end
                eventShapeFactor =  [eventShapeFactor evShapeFactor];
                minInternalThroatRadiusE = 0;
                
                %Calculations of min internal throat radius
                if ~isnan(eventThroats)
                    if(min(linkRadius(eventThroats)))                
                        minInternalThroatRadiusE=min(linkRadius(eventThroats));                
                    end
                    minInternalThroatRadius = [minInternalThroatRadius minInternalThroatRadiusE];  
                end
            
                if isnan(eventThroats)
                    minInternalThroatRadius = [minInternalThroatRadius nan];  
                end
            
            end            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Contact angle calculations

            eventCaIn = CaMatrix(poreEvent{j}, CaColumnValue);
            eventCaStdIn = CaMatrix(poreEvent{j}, CaColumnValue);           
            eventCaWeightsIn = CaMatrix(poreEvent{j}, CaWeightColumnValue);         
            
            eventCaAll = CaMatrix(allEventPores, CaColumnValue);
            eventCaStdAll = CaMatrix(allEventPores, CaColumnValue);
            eventCaWeightsAll = CaMatrix(allEventPores, CaWeightColumnValue);

            avEventCaIn = [avEventCaIn nansum(eventCaIn.*eventCaWeightsIn)/nansum(eventCaWeightsIn)];
            avEventCaStdIn = [avEventCaStdIn nanmean(eventCaWeightsIn)];
                                 
            avEventCaAll = [avEventCaAll nansum(eventCaAll.*eventCaWeightsAll)/nansum(eventCaWeightsAll)];
            avEventCaStdAll = [avEventCaStdAll nanmean(eventCaWeightsAll)];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %keeping track of the event numbers
            eventNumber = eventNumber+1;
        end
  

        
    else
        nMultiPoreEvents(t) = 0;
    end
    
        

end
disp('Finished loop')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Calculating event properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mean radius of curvature
avEventRadiiIn = (1./avEventCurvatureIn) *1e-6 ;
maxEventRadiiIn = (1./maxEventCurvatureIn) *1e-6 ;
avEventRadiiAll = (1./(avEventCurvatureAll)) *1e-6 ;
maxEventRadiiAll = (1./maxEventCurvatureAll) *1e-6 ;

disp('Finished event properties')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating normalised filling volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normPoreVolFrac = (eventVolumesAll/totalPoreVolume);

%Calculates the cumulative filling volume 
for i = 1:length(cumulFillingVol)
    cumulFillingVol(i) =  sum(fillingVol(1:i));
end

% some statistics about events
nansum((normPoreVolFrac > 0.01))
 nansum((normPoreVolFrac > 0.01).* normPoreVolFrac);
% nansum((normPoreVolFrac));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the saturation value when each event occurred. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saturation_event_Imb = [];
for i = 1:length(eventT)
    saturation_event_Imb(i) = saturation(eventT(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real time per time step and real time per event.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real_time_imb = [145;258;372;485;598;711;824;938;1051;1164;1278;1391;1504;1731;1844;1957;2070;2184;2297;2410;2523;2863;3034;3147;3260;3373;3487;3600;3713;3940;4054;4167;4280;4393;4506;4620;4733;4846;4960;5073;5186;5300;5413;5526;5640;5753;6144;6485;6824;7050;7503;9581];
real_eventT_imb = [];
for i = 1:length(eventT)
    real_eventT_imb(i) = real_time_imb(eventT(i));
end
real_fillingTPerPore_imb = [];
for i = 1:length(fillingTPerPore)
    if isnan(fillingTPerPore(i))
        real_fillingTPerPore_imb(i)= NaN;
    else
        real_fillingTPerPore_imb(i) = real_time_imb(fillingTPerPore(i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display("Finished")
