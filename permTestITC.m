%% permutation test for ITC maps
clear ITCall randItcDiff

%numFreq = 9;
numFreq = 6;
FS
timeReduction = [200:700]; %change according to the NaNs in the data matrix. Try to get rid of as many as possible. 
numTime = length(timeReduction);
itcTD_s = itcTD(:,:,timeReduction);
itcASD_s = itcASD(:,:,timeReduction);
%TDitc = rand(20,numFreq,numTime); % TD group itc maps

%ASDitc = rand(30,numFreq,numTime); % ASD group itc maps
realDiffItc = squeeze(mean(itcTD_s,1) - mean(itcASD_s,1)); % the real difference in ITC between TD and ASD groups
ITCall = cat(1, itcTD_s, itcASD_s);
%%

N = 10000; % number of permutations

% loop over permutations
nASD = 31; % number of subjects in the ASD group
nTD = 20; % number of subjects in the TD group
randItcDiff = zeros(N,numFreq,numTime); % matrix to hold the random differences results
for n = 1 : N
    %disp(n)
    % randomly assign subjects to TD or ASD groups
    randSubj = randperm(size(ITCall,1));
    TDitcRand = ITCall(randSubj(1:nTD),:,:);
    ASDitcRand = ITCall(randSubj(nASD+1:end),:,:);
    randItcDiff(n,:,:) = squeeze(mean(TDitcRand,1) - mean(ASDitcRand,1));
end

% plot example histogram for a single time and frequency point
figure
histogram(squeeze(randItcDiff(:,3,300)))
hold on
line([realDiffItc(3,300) realDiffItc(3,300)],[0 600],'Color','r')

% find the position of the real data compared to the random distribution
% for each time and frequency point
pVal = zeros(numFreq, numTime);
for freqInd = 1 : numFreq
    for timeInd = 1 : numTime
        pVal(freqInd, timeInd) = sum(squeeze(randItcDiff(:,freqInd,timeInd)) > realDiffItc(freqInd, timeInd))/N;
        %pVal(freqInd, timeInd) = (length(find(abs(squeeze(randItcDiff(:,freqInd,timeInd))) > abs(realDiffItc(freqInd, timeInd))))+1)/(N+1);
    end
end
pVal(pVal<0.0001)=nan;

sum(sum(pVal<0.0505))/numel(pVal)*100  % percentage of significant p from the whole map

% set a threshold for 10 consecutive 
pValSig = pVal<0.0505;
tooSmall= find(S<10 & S>0); %find the row in the matrix where there are too few (<10) significant p vals
pValSig(tooSmall,:) = 0;
pVal(tooSmall,find(pVal(5,:)<0.0505)) = 0.06;
%
%% Plot a map of p values - black & white

lines = PARAMS.lines;
Nint = 500
Timeoi = timeoi(timeReduction);
[x,y] = meshgrid(Timeoi,freqoi); % low-res grid
[x2,y2] = meshgrid(Timeoi(1):1/Nint/5:Timeoi(end),freqoi(1):.01:freqoi(end));  %high-res grid
dataInterp_p = interp2(x,y,pVal, x2,y2, 'linear'); %interpolate up 


figure;
h =  imagesc(Timeoi(1):1/Nint/5:Timeoi(end),freqoi(1):.01:freqoi(end),dataInterp_p)
pSig = dataInterp_p < 0.0505;
    hold on; 
    contour(Timeoi(1):1/Nint/5:Timeoi(end),freqoi(1):.01:freqoi(end),pSig,'r', 'LineWidth',1.5);
    colormap gray;
pSig2 = dataInterp_p < 0.03;
hold on; 
contour(Timeoi(1):1/Nint/5:Timeoi(end),freqoi(1):.01:freqoi(end),pSig2,'g', 'LineWidth',1.5);
colorbar;
    caxis([0 1.12]);
    set(gca,'YDir','normal')
    for k = 2:length(lines)
        line([lines(k), lines(k)],[y2(1,1),y2(end,1)]...
        ,'Color','w','LineWidth',2,'LineStyle','--');
    end
    ylabel('Frequency (Hz)');
    xlabel('Time (Sec.)');
    title('P values for TD/ASD ITPC differences')
    set(gca,'fontsize', 14);
    
%% correcting for multiple comparisons - only the significant p

[row,col,v] = find(dataInterp_p(dataInterp_p < 0.05));

    
    
    
    
    
    
    %% averaging across the cues
clear ITCall_short itcTD_short itcASD_short
timeStartInds = [70,237,403]; %indices 
TIMEWIN= 0.1;
FREQWIN = freqoi; 

for i = 1:3 
    ITCall_short{i} = ITCall(:,:,round(timeStartInds(i)-TIMEWIN*FS):round(timeStartInds(i)+TIMEWIN*FS));
    itcTD_short{i} = itcTD_s(:,:,round(timeStartInds(i)-TIMEWIN*FS):round(timeStartInds(i)+TIMEWIN*FS));
    itcASD_short{i} = itcASD_s(:,:,round(timeStartInds(i)-TIMEWIN*FS):round(timeStartInds(i)+TIMEWIN*FS));

end
itcTD_short_mean = mean(cat(4,itcTD_short{1},itcTD_short{2},itcTD_short{3}),4);
itcASD_short_mean = mean(cat(4,itcASD_short{1},itcASD_short{2},itcASD_short{3}),4);
realDiffItc_s = squeeze(mean(itcTD_short_mean,1) - mean(itcASD_short_mean,1)); % the real difference in ITC between TD and ASD groups

ITCall_short_mean = mean(cat(4,ITCall_short{1},ITCall_short{2},ITCall_short{3}),4);


% run permutation test on the small window

clear randItcDiff
N = 10000; % number of permutations
numTime = size(ITCall_short_mean,3)
% loop over permutations
nASD = 31; % number of subjects in the ASD group
nTD = 20; % number of subjects in the TD group
randItcDiff = zeros(N,numFreq,numTime); % matrix to hold the random differences results
for n = 1 : N
    % randomly assign subjects to TD or ASD groups
    randSubj = randperm(size(ITCall_short_mean,1));
    TDitcRand = ITCall_short_mean(randSubj(1:nTD),:,:);
    ASDitcRand = ITCall_short_mean(randSubj(nASD+1:end),:,:);
    randItcDiff(n,:,:) = squeeze(mean(TDitcRand,1) - mean(ASDitcRand,1));
end

% find the position of the real data compared to the random distribution
% for each time and frequency point
pVal = zeros(numFreq, numTime);
for freqInd = 1 : numFreq;
    for timeInd = 1 : numTime
        pVal(freqInd, timeInd) = sum(squeeze(randItcDiff(:,freqInd,timeInd)) > realDiffItc_s(freqInd, timeInd))/N;
        %pVal(freqInd, timeInd) = (length(find(abs(squeeze(randItcDiff(:,freqInd,timeInd))) > abs(realDiffItc(freqInd, timeInd))))+1)/(N+1);
    end
end
%pVal(pVal<0.0001)=nan;

sum(sum(pVal<0.05))./numel(pVal)*100

%% trim the pVal matrix 
pVal = pVal(2:4,:);
%pValSig = pVal(pVal>0.05)=1;
[q,r] = find(pVal<0.05);
%% correction for multiple comparison

 [H, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pVal,0.05,'pdep','yes');   
    
 reshaped = reshape(pVal,1,[])

 [FDR,Q,aPrioriProb]  = mafdr(reshaped,'Showplot',true)   
    
 Psort = sort(reshaped);
 Psort(2,:) = FDR;
 Psort(3,:) = Q;
  %  clear wind TOI_win
%  TIMEWIN = 0.15; % in sec, time window before and after each cue, to consider for stat test
%  FREQWIN = freqoi; %[2:4];
%  % limit to known windows
%  
%  timeStartInds = [70,237,403]; %indices 
%  % isolate time windows around each cue. 

%  for i = 1:length(timeStartInds)
%     TOI_win{i} = pVal(FREQWIN,round(timeStartInds(i)-TIMEWIN*FS):round(timeStartInds(i)+TIMEWIN*FS));
%     %pSigWin{i} = TOI_win{i}<0.05;
%     %[H{i}, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(TOI_win{i},0.05,'pdep','yes');
%  end
%  
 %reshaped = reshape(TOI_win{2},1,[]);
 [FDR,Q,aPrioriProb]  = mafdr(reshaped,'Showplot',true)

 %average all areas around the cues
%meanWind = mean(cat(3,TOI_win{1},TOI_win{2},TOI_win{3}),3);

%% plot interpolated mean of p values around the cue times
%Timeoi_m = round(timeStartInds(i)-TIMEWIN*FS):round(timeStartInds(i)+TIMEWIN*FS);
Timeoi_m = 0:1/FS:numTime/FS;

[x,y] = meshgrid(Timeoi_m,FREQWIN); % low-res grid
[x2,y2] = meshgrid(Timeoi_m(1):1/Nint/5:Timeoi_m(end),FREQWIN(1):.01:FREQWIN(end));  %high-res grid
dataInterp_p = interp2(x,y,pVal, x2,y2, 'linear'); %interpolate up 

figure;
h =  imagesc(Timeoi_m(1):1/Nint/5:Timeoi_m(end),pVal(1):.01:FREQWIN(end),dataInterp_p)
colormap gray;
set(gca,'YDir','normal')
[H, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(meanWind,0.05,'pdep','yes');
 
%  for k = 2:4
%      wind{k} = pVal(:,round(WOI(k)-TIMEWIN*FS-timeReduction(1)):...
%          round(WOI(k)+TIMEWIN/2*FS)-timeReduction(1)); %location of pVal in samples
%  
%      pSigWin = wind{k}<0.05;
%      [h{k}, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pSigWin,0.05,'pdep','yes');
%  end
%  
%  meanWind = mean(cat(3,wind{1},wind{2},wind{3}))
%  
% [FDR,Q,aPrioriProb]  = mafdr(reshaped,'Showplot',true)
%  

% Timeoi_win = TOI_win(1,:);
% 
% [x,y] = meshgrid(Timeoi_win,freqoi); % low-res grid
% [x2,y2] = meshgrid(Timeoi_win(1):1/Nint/5:Timeoi_win(end),freqoi(1):.01:freqoi(end));  %high-res grid
% 
% dataInterp_p_win = interp2(x,y,h{2}, x2,y2, 'linear'); 
%  