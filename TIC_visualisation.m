%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script to visualise EICs from TAPP raw data
%
% Authors: Peter Horvatovich, Vikram Mitra and Gooitzen Zwanenburg
% Date: July 17, 2015
%
% set mainDir variable to exact location of the ASCA main directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
% change the directory according to the exact coce location
mainDir = 'e:\data\publications\2012\vikram_fact_des\code\codeOkPerm\';
% mainDir = 'd:\users\peter\munka\publications\2012\vikram_fact_des\code\codeOk\';
cd(mainDir)
fileName='TicVisualisation\TIC_16FactDes.txt';%List of EICs
fileTmap='TicVisualisation\inputTIC_16FactDes_warp.txt';% list of .tmap files
imageOutputDir = 'TICFigures\';%output directory of EICs plots
dummyDir = exist(imageOutputDir,'dir');
if dummyDir ~= 7,
    mkdir(imageOutputDir)
end
clear dummyDir

load('peakList.mat')
load('factorLabel.mat')

f = importdata('FMatrix.txt'); % experiment design file (the same as the original fractional factorial design)
f([4 11 14],:) = [];
f = f([2,7,5,14,1,15,3,6,12,10,11,8,13,16,4,9],:);

formatSpec = '%s%f%[^\n\r]';
fileID = fopen(fileTmap,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', ' ', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
tmapFileName = dataArray{:, 1};
class = zeros(size(f,1),1);
clearvars filename formatSpec fileID dataArray ans;
legendFlag=0;
rtCenter=peakList(:,3);
rtHalfWindow=30;

fid = fopen([mainDir fileName]);
tline = fgetl(fid);
lineCounter=1;
TICCounter=1;
TIC=struct('Intensity',[],'MS',[],'rt',[],'param',[],'numberOfTIC',[],'TIC',[],'rtWarped',[]);
TIC.param=struct('m1',[],'m2',[],'t1',[],'t2',[],'Mass',[],'Delta',[], 'FileName', [], 'tmapFileName', []);
color=[1,0,0;0,0,1;0,1,0;1,1,0;0,1,1;1,0,1;0,0,0];
factSel = [2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6];
peaksNum = [1:10,1:10,1:10];
nonMZ=30;
nonClass=length(f);

while strcmp(tline,'Intensity'),
    temp=[];
    lineCounter=lineCounter+1;
    while not(strcmp(tline,'Mass')),
        dummy=[];
        dummy = sscanf(tline,'%e');
        temp = [temp;dummy'];
        tline = fgetl(fid);
        lineCounter=lineCounter+1;
    end
    TIC.Intensity{TICCounter}=temp;
    tline = fgetl(fid);
    TIC.MS{TICCounter} = sscanf(tline,'%e');
    tline = fgetl(fid);
    tline = fgetl(fid);
    TIC.rt{TICCounter} = sscanf(tline,'%e');
    dummy=[];
    tline = fgetl(fid);
    tline = fgetl(fid);
    dummy = sscanf(tline,'%e');
    TIC.param.m1{TICCounter} = dummy(1);
    TIC.param.m2{TICCounter} = dummy(2);
    TIC.param.t1{TICCounter} = dummy(3);
    TIC.param.t2{TICCounter} = dummy(4);
    dummy = [];
    tline = fgetl(fid);
    dummy = sscanf(tline,'%e');
    TIC.param.Mass{TICCounter} = dummy(1);
    TIC.param.Delta{TICCounter} = dummy(2);
    tline = fgetl(fid);
    tline = fgetl(fid);
    dummy = sscanf(tline,'%s');
    TIC.param.FileName{TICCounter} = dummy;
    tline = fgetl(fid);
    tline = fgetl(fid);
    TIC.TIC{TICCounter}=sum(TIC.Intensity{TICCounter});
    TICCounter=TICCounter+1;
end
TIC.numberOfTIC=TICCounter-1;
fclose(fid);

TICCounter=1;
for i=1:nonMZ,
    for ii=1:nonClass,
        delimiter = '\t';
        formatSpec = '%f%f%[^\n\r]';
        fileID = fopen(tmapFileName{ii},'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
        fclose(fileID);
        TIC.param.tmapFileName{i+nonMZ*(ii-1)} = tmapFileName{ii};
        originalTmap = dataArray{:, 1};
        warpedTmap = dataArray{:, 2};
        clearvars filename delimiter formatSpec fileID dataArray ans;
        dummyRt=TIC.rt{i+nonMZ*(ii-1)};
        TIC.rtWarped{i+nonMZ*(ii-1)}=interp1(originalTmap,warpedTmap,dummyRt,'linear');
        TICCounter=TICCounter+1;
    end
end

TICCounter=1;
warpedFlag=1;
legendFlag=0;
FontSize = 12;
for i=1:nonMZ,
    class = f(:,factSel(i),:);
    fHandle=figure;
    hold on
    set(gca,'FontSize',FontSize);
    title([num2str(TIC.param.Mass{i+nonMZ*(ii-1)}),'±',num2str(TIC.param.Delta{i+nonMZ*(ii-1)},'%1.3f'),' Da, ',num2str(rtCenter(i)),' min'])
    legendstring=cell(nonClass,1);
    for ii=1:nonClass,
        if warpedFlag,
            rtVec=TIC.rtWarped{i+nonMZ*(ii-1)};
        else
            rtVec=TIC.rt{i+nonMZ*(ii-1)};
        end
        if class(ii)==1,
            plotString='-';
            stringLegend=[factorLabel{factSel(i)} ' ' num2str(peaksNum(i)) '. peak'];
        else
            plotString='-';
            stringLegend=[factorLabel{factSel(i)} ' ' num2str(peaksNum(i)) '. peak'];
        end
        plot(rtVec,TIC.TIC{i+nonMZ*(ii-1)},plotString,'color',color(class(ii)+1,:), 'LineWidth', 1.25);
        legendstring{ii}=stringLegend;
        TICCounter=TICCounter+1;
    end
    axis tight
    dummyAxis = axis;
    axis([rtCenter(i)-rtHalfWindow, rtCenter(i)+rtHalfWindow, dummyAxis(3:4)])
    xlabel('Time (min)');
    ylabel('Ioncounts (cts)');
    axisPos = axis;
    X = axisPos(1)+(axisPos(2)-axisPos(1))*0.60;
    Y = axisPos(3)+(axisPos(4)-axisPos(3))*0.90;
    text(X,Y,[factorLabel{factSel(i)} ' ' num2str(peaksNum(i)) '. peak'],'FontSize',FontSize)
    if legendFlag,
        legend(legendstring);
    end
    saveas(fHandle,[imageOutputDir, num2str(TIC.param.Mass{i+nonMZ*(ii-1)}), '_Da_', num2str(TIC.param.Delta{i+nonMZ*(ii-1)},'%1.3f'),'_Da_',num2str(rtCenter(i)),'_min.emf'],'emf')
    saveas(fHandle,[imageOutputDir, num2str(TIC.param.Mass{i+nonMZ*(ii-1)}), '_Da_', num2str(TIC.param.Delta{i+nonMZ*(ii-1)},'%1.3f'),'_Da_',num2str(rtCenter(i)),'_min.fig'],'fig')
    close(fHandle)
end
