%% BIONEGG 2586 - Assignment - 1 Reaching Assignment
% Assignment done in group by Bobby Mohan and Sai Sravanthi Joshi
%% Loading the Data
clc;
close all;
load('neuraldata.mat')
% To get the standard target location
TP = [R.TrialParams];
allTarget1X = [TP.target1X];
allTarget1Y = [TP.target1Y];

% plot gives the target location
scatter(allTarget1X, allTarget1Y);
hold on;
HandX = [TP.touchPointX];
HandY = [TP.touchPointY];
scatter(HandX,HandY);
hold on;
% Ploting the initial eye position
scatter(20,20);

%% PLot showing the hand position of all the trials
figure();
plot([R.hhp],[R.vhp]);

% eliminating the outliers which is the trial 1739 by binary process
R_left=R(1:1738);
R_right=R(1740:end);
RN=cat(2,R_left,R_right);

%plot shows the animal hand position without outliers


figure();
plot([RN.hhp],[RN.vhp]);
     
    hold on;
    scatter(allTarget1X, allTarget1Y,'LineWidth',8);
    hold on;

plot([0,-34],[0,93],'color','red','LineWidth',3);
plot([0,-86],[0,50],'color','yellow','LineWidth',3);
plot([0,-98],[0,-17],'color','magenta','LineWidth',3);
plot([0,-64],[0,-76],'color','white','LineWidth',3);
plot([0,64],[0,-76],'color','black','LineWidth',3);
plot([0,98],[0,-17],'color','cyan','LineWidth',3);
plot([0,86],[0,50],'color','green','LineWidth',3);
plot([0,34],[0,93],'color','#50F','LineWidth',3);
% hold on;
% scatter([RN.mhep],[RN.mvep]);


%%
% Calculation of Perimovement period
for i = 1:1848
RN1(i).hhpN=[RN(i).hhp(:,(RN(i).timeGoCue):(RN(i).timeTargetAcquire))];
RN1(i).vhpN=[RN(i).vhp(:,(RN(i).timeGoCue):(RN(i).timeTargetAcquire))];
RN1(i).TimeN=[(RN(i).timeGoCue):(RN(i).timeTargetAcquire)];
RN1(i).hepN=[RN(i).hep(:,(RN(i).timeGoCue):(RN(i).timeTargetAcquire))];
RN1(i).vepN=[RN(i).vep(:,(RN(i).timeGoCue):(RN(i).timeTargetAcquire))];
end

% Sample trial of perimovement removal
figure();
subplot(2,1,1)
plot(RN(100).hhp,RN(100).vhp);
subplot(2,1,2)
plot(RN1(100).hhpN,RN1(100).vhpN);
hold on;
scatter(TP(100).target1X, TP(100).target1Y);

%% Calculating the velocity of the hand movement

for i=1:1848
    RN1(i).Distance=[sqrt((RN1(i).hhpN).^2 + (RN1(i).vhpN).^2)];
    RN1(i).Velocity=[RN1(i).Distance]./[RN1(i).TimeN]; % the units obtained for velocity here is m/s
    RN1(i).Maxvelocity=max(RN1(i).Velocity);  % the units obtained for max velocity here is m/s
    RN1(i).tenperVel= RN1(i).Maxvelocity/20;
    [RN1(i).val, RN1(i).idx]=max(abs((RN1(i).Velocity)-(RN1(i).tenperVel)));
    RN1(i).timeMovBegins=(RN1(i).idx)+(RN(i).timeGoCue);
    RN1(i).ReactionTime=(RN1(i).timeMovBegins)-(RN(i).timeGoCuePHOTO);
    RN1(i).TargetNumbers=allTarget1X(i);
    RN1(i).TargetNumbers=changem(RN1(i).TargetNumbers, [1 2 3 4 5 6 7 8],[-34 -86 -98 -64 64 98 86 34]);
    RN1(i).latency=(RN(i).timeGoCuePHOTO)-(RN(i).timeGoCue); % the latency value here is in ms
end

% Ploting all the trials with different colours
A = cat(RN1.TargetNumbers);
figure();
for i=1:1:8
    hold on;
    hhpT(i,:)=(RN1(A==i).hhpN);
    vhpT(i,:)=(RN1(A==i).vhpN);
%     plot((RN1(A==i).hhpN),(RN1(A==i).vhpN));
end


%% Statistical Analysis of Reaction Time based on each target
%BOX PLOT -> Reaction Time Vs Target-> Up - slow-> Down - Fast
figure();
boxplot([RN1.ReactionTime],[RN1.TargetNumbers]);
MeanReactionTime=mean([RN1.ReactionTime]);
StdReactionTime=std([RN1.ReactionTime]);

%% Calculating the monitor latency
figure();
boxplot([RN1.latency]);
MeanLatency=mean([RN1.latency]);
Refreshingrate=(1/MeanLatency)*1000; % 52.88 Hz is the refershing rate
%% Neural Analysis - Part II
% Eye movement 
Sam_Num=1;
figure();
subplot(2,2,1)
plot(RN(Sam_Num).hep,RN(Sam_Num).vep); % Raw data of eye movement
subplot(2,2,2)
plot(RN1(Sam_Num).hepN,RN1(Sam_Num).vepN); % Eye movement trimed to perimovement period
subplot(2,2,3)
plot(smooth(RN1(Sam_Num).hepN),smooth(RN1(Sam_Num).vepN)); % Smoothen Eye movement
subplot(2,2,4)
plot(RN1(Sam_Num).hhpN,RN1(Sam_Num).vhpN);
%% Neural Analysis on individual trial

figure();
subplot(2,1,1);
plot(RN1(Sam_Num).hhpN,RN1(Sam_Num).vhpN);
hold on;
scatter((RN1(Sam_Num).hhpN(1)),(RN1(Sam_Num).vhpN(1)));
hold on;
Diff=(RN(Sam_Num).timeGoCuePHOTO)-(RN(Sam_Num).timeGoCue);
scatter(((RN1(Sam_Num).hhpN(1+Diff))),((RN1(Sam_Num).vhpN(1+Diff))));
subplot(2,1,2);
plot(smooth(RN1(Sam_Num).hepN),smooth(RN1(Sam_Num).vepN));

%%
figure();
%Normalization of the hand position
Nor_hhp=normalize([RN1(1).hhpN],'range',[0 4]);
plot([RN1(1).TimeN],Nor_hhp);
hold on;
plot((RN(1).timeGoCue),Nor_hhp(1),'*');
hold on;
TimeIndex=(RN(1).timeGoCuePHOTO)-(RN(1).timeGoCue);
plot((RN(1).timeGoCuePHOTO),(Nor_hhp(TimeIndex)),'*');

% COde to plot Raster

for itrial=1:4
    spks=(RN(1).cells(itrial).spikeTimes)';
    xspikes=repmat(spks,3,1);
    yspikes=nan(size(xspikes));
    
    if (~isempty(yspikes))
        yspikes(1,:)=itrial-1;
        yspikes(2,:)=itrial;
    end
    hold on;
    plot(xspikes,yspikes,'color','k');
    
end
        
%% To plot PSTH curve
% psth(times, binsize, fs, ntrials, triallen, varargin);

% spiketimes=RN(1).cells;
% binsize=(((RN(1).timeCueOnset)-300):((RN(1).timeCueOnset)+600));
% fs=1000;
% ntrials=4;
% triallen=0;
% varargin=0;
% psth(spiketimes, binsize, fs, ntrials, triallen, varargin);

% TIMES - spike times (samples)
% BINSIZE - binwidth (ms)
% FS - sampling rate (hz)
% NTRIALS - number of trials
% TRIALLEN - length of a trial (samples)




 
