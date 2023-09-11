%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to set up the data for the eventClustering script

% By Tom Bultreys, Arjen Mascini and Sharon Ellman
% Adapted for Schluter et al. data no pad
% Runs with statoil format of pn extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Some variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imageSize = 751*751*576;
numTimeSteps = 52;
numPores = 604;%%%%%%%%%%%%%%%%%%%%%%%%%%% Will change with new pores
numThroats = 2162;%%%%%%%%%%%%%%%%%%%%%%%%%%% Will change with new pores
interfacialTension = 36 *1e-3; % in N/m
totalPoreVolume = 0.000000047402; %m3, from image
% porosity = 0.22;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% files to load 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Will change with new pores
pnmBaseFileName = './data/network_data/schluter_ICL'; % location of the statoil formatted PNM files
load('./data/maxBall_data/pores_majorityGV.mat'); % majority in pore maximum inscribed spheres saved in .mat format
load('./data/maxBall_data/throats_majorityGV.mat'); % majority in throat maximum inscribed spheres
load('./data/volume_data/Mean_fillFrac.mat'); % mean filing fraction
load('./data/maxBall_data/ThroatMean.mat');
load('./data/maxBall_data/throatSurfaces_Majority.mat')
load('./data/maxBall_data/throatSurfaces_Mean.mat')

poreVolumesAvizo = table2array(readtable('./data/volume_data/poreVolumeAvizo.csv','ReadVariableNames',false)); 
pcGlobal =table2array(readtable('./data/global_data/pcGlobal.csv','ReadVariableNames',true)); %Pc global csv 
fullKc =  -table2array(readtable('./data/global_data/Kc_geod.csv','ReadVariableNames',false)); %clipped weighted average Kc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kc Ca Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clippingValue = '0.2';
prefix_kcName = './data/SurfaceData/Kc_files_curvatureData/perPore_geodesic/';  
postfix_kcName = 'KcOutputFile.txt';
KcFileOffset = 2; %because first two timesteps don't have datapoints, starts at 0 and at T-1 = 4 
curvatureColumnValue = 13; %column 3 for curvatures, 13 for weighted and cropped curvatures
weightColumnValue = 2;
stdColumnValue = 9;

prefix_CaName = './data/contactAngle_data/';
postfix_CaName = 'CaOutputFile.txt';
CaColumnValue = 3;
CaWeightColumnValue = 2;
CaStdColumnValue = 4;

CA_Global = table2array(readtable('./data/global_data/s00051.csv','ReadVariableNames',false)); 
CA_Global = CA_Global(:,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%`Folder for plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path2 = strcat('./Plots/plots_c',clippingValue,'_paper/'); 
if ~exist(path2, 'dir')
   mkdir(path2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%`Filling files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = 0.5; % oil = 1, brine+rock = 0

ePoreFilling = reshape(PoreMajority, [numPores, numTimeSteps]); % reshape the data file matrix to pores x timesteps 
ePoreFilling = ePoreFilling < threshold; % less than threshold. Reverses 1s and 0s. So brine =1 and oil = 0

eThroatFilling = reshape(ThroatMajority, [numThroats, numTimeSteps]); %same for throats
eThroatFilling = eThroatFilling < threshold; %Reverses 1s and 0s. So brine =1 and oil = 0
eThroatFilling_snapoff = (ThroatMean<0.3);

eThroatSurfaceFilling = reshape(throatSurfaces_Majority, [numThroats, numTimeSteps]); %same for throats
eThroatSurfaceFilling = eThroatSurfaceFilling < threshold; %Reverses 1s and 0s. So brine =1 and oil = 0
eThroatSurfaceFilling_snapoff = (throatSurfaces_Mean<0.3);

ePorefillingFraction = reshape(Mean, [numPores, numTimeSteps]);
eThroatfillingFraction = reshape(ThroatMean, [numThroats, numTimeSteps]);

%ePorefillingFraction = ePorefillingFraction(2:length(ePorefillingFraction),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%`node data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodeCoordNum = dlmread([pnmBaseFileName '_node1.dat'], '', [1 4 numPores 4] );
nodeCoords = dlmread([pnmBaseFileName '_node1.dat'], '', [1 1 numPores 3] );
nodeRadius = dlmread([pnmBaseFileName '_node2.dat'], '', [0 2 numPores-1 2] );
nodeVolume = dlmread([pnmBaseFileName '_node2.dat'], '', [0 1 numPores-1 1] );
nodeShapeFactor = dlmread([pnmBaseFileName '_node2.dat'], '', [0 3 numPores-1 3] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% link data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
links = dlmread([pnmBaseFileName '_link1.dat'], '', [1 1 numThroats 2] );
linkRadius = dlmread([pnmBaseFileName '_link1.dat'], '', [1 3 numThroats 3] );
linkVolume = dlmread([pnmBaseFileName '_link2.dat'], '', [0 6 numThroats-1 6] );
linkShapeFactor = dlmread([pnmBaseFileName '_link1.dat'], '', [1 4 numThroats 4] );
linkTotalLength = dlmread([pnmBaseFileName '_link1.dat'], '', [1 5 numThroats 5] );
linkAspectRatio = zeros(numThroats, 1);

%removes all inlet-outlet nodes and links
linksClean = links(links(:,1)>0 & links(:,2)>0, :);
eThroatFillingClean = eThroatFilling(links(:,1)>0 & links(:,2)>0, :);
eThroatFilling_snapoffClean = eThroatFilling_snapoff(links(:,1)>0 & links(:,2)>0, :);
linkRadiusClean = linkRadius(links(:,1)>0 & links(:,2)>0, :)

eThroatSurfaceFillingClean = eThroatSurfaceFilling(links(:,1)>0 & links(:,2)>0, :);
eThroatSurfaceFilling_snapoffClean = eThroatSurfaceFilling_snapoff(links(:,1)>0 & links(:,2)>0, :);

eThroatfillingFraction
eThroatfillingFractionClean = eThroatfillingFraction(links(:,1)>0 & links(:,2)>0, :);


% finds all connecting pores to throats
throatNeighbours = cell(numPores,1);
for thr=1:length(links)
    p1 = links(thr,1);
    p2 = links(thr,2);
    if(p1>0 && p2>0)
        throatNeighbours{p1} = [throatNeighbours{p1} thr];
        throatNeighbours{p2} = [throatNeighbours{p2} thr];
    end
    clear p1 p2
end
clear thr

%find all neighboring pores 
poreNeighbours = cell(numPores,1);
for thr=1:length(links)
    p1 = links(thr,1);
    p2 = links(thr,2);
    if(p1>0 && p2>0)
        poreNeighbours{p1} = [poreNeighbours{p1} p2];
        poreNeighbours{p2} = [poreNeighbours{p2} p1];
    end
    clear p1 p2
end
clear thr

% calculates the aspect ratio of each pore
for i=1:numThroats
    
    if(links(i,1)> 0 && links(i,2)>0)
        linkAspectRatio(i) = 0.5*(nodeRadius(links(i,1))+nodeRadius(links(i,2)) ) ./ linkRadius(i);
    else
        linkAspectRatio(i) = 0;
    end
end

%another definition of aspect ratio
nodeAspectRatio = zeros(numPores,1);
for i=1:numPores   
    if ~isempty(throatNeighbours{i})
        nodeAspectRatio(i) = nodeRadius(i) ./ max(linkRadius(throatNeighbours{i}));
    end
end


% check to see the difference in saturation and prints this
saturationTimeZero = (sum(ePoreFilling(:,1).*nodeVolume) + sum(eThroatFilling(:,1).*linkVolume)) ./ (sum(nodeVolume) + sum(linkVolume));
saturationTimeEnd = (sum(ePoreFilling(:,end).*nodeVolume) + sum(eThroatFilling(:,end).*linkVolume)) ./ (sum(nodeVolume) + sum(linkVolume));

disp('Set up complete')
disp(["Water Saturation (Sw) evolved from " saturationTimeZero " to " saturationTimeEnd])


