%clear itcAll dataTemp spectAll spectrumEst itcGroup DATA_I sumItc DATA_reref data
% analysed for the supp figure 1 in VAMP. 

% Shlomit Beker 2020-2023

%% calculate TFR (and ITC) for each channel 
clear tfrAll dataTemp spectrumEst tfr Group data DATA_r DATA_I sumtfr

Group = 4;
% tic
gwidthWlt = 3;
freq = 1:45; %1:8 2:7;

widthWlt = linspace(3,10,45);
CHAN = 'PO3'; %used in the itpc time-freq plot
REF = 'all';
chns = find(strcmp(ERPb{1}{1}.label,CHAN)); % enter the ID of channels to analyze
LENGTH = size(ERP{1}{1}.avg,2);
TIMEBEFORE_STIM = 1;
t = [1:LENGTH]/256-TIMEBEFORE_STIM;
refCHAN = find(strcmp(ERPb{1}{1}.label,REF));

% calculate TFR
% try one to initialize size of timeoi and freqoi
% data is a fieldtrip structure with single trial data
% inputs: 
% dataTemp - single trial data matrix (chans X time)
% t - time vector in seconds

dataTemp = DATA{1}{1}{1}(chns,:);
[~,freqoi,timeoi] = ft_specest_wavelet(dataTemp, t, 'freqoi', freq,'width', widthWlt, 'gwidth',gwidthWlt);
% DATA_r = rerefData(DATA,refCHAN); 
DATA_r = DATA;
      clear DATA_I
      for k = 1:length(DATA_r{Group})
          for l = 1:length(DATA_r{Group}{k})
              DATA_I{k}(:,:,l) = DATA_r{Group}{k}{l}; % DATA_I includes only participants in the analysed datasets. Number of cells = number of participants. 
          end
      end
 %numTrl = size(DATA_I{Group},3);
 numTrl = size(DATA_I,3);
       
     %for CHAN = chns; %1:length(ERPb{1}{1}.label);
        for k =  1:length(DATA_I)  % loop on subjects
            % initialize matrix for spectrogram
            %for chnI = 1 : length(chns)
            % this matrix holds the instantaneous phase for each trial, frequency
            % and time point
            numTrl = size(DATA_I{k},3);
            spectAll = zeros(numTrl,length(freqoi),length(timeoi));
            for trlI = 1 : numTrl
                % select the correct trial and channel
                %dataTemp = data.trial{trlI}(chnI,:);
                dataTemp = DATA_I{k}(chns,:,trlI);
                % calculate time frequency analysis using wavelets
                [spectrumEst,freqoi,timeoi] = ft_specest_wavelet(dataTemp, t, 'freqoi', freqoi, 'width', widthWlt, 'gwidth',gwidthWlt);
                spectAll(trlI,:,:) = spectrumEst(1,:,:);
                
            end
            spectAll_abs = abs(spectAll);
            spectAll_abs_mean = squeeze(mean(spectAll_abs,1)); %2d matrix of one individual subject
            
            %itc{k} = it_calcITC(spectAll); comment out if only tfr
            tfr{k}(1,:,:) = spectAll_abs_mean;
            %tfr_sub{k} = spectAll_abs; %calculate per subject
        end
        
        sumtfr = zeros(length(freqoi),length(timeoi));
        nans = 0;
        for k = 1:length(tfr) %for averaging across group
             sqtfr = squeeze(tfr{k});
%              if sum(sum(isnan(sqtfr))) == size(sqtfr,1)*size(sqtfr,2)
%                  sqtfr = [];
%                  nans = nans+1;
%                  k
%                  continue
%              end
             sumtfr = sumtfr+sqtfr;
        end
        
%%
    if Group == 1
            tfrAllTd = sumtfr./length(DATA_r{Group});
            data = tfrAllTd;
            tfrTD = cell2mat(tfr');
        elseif Group == 4
            tfrAllAsd = sumtfr./length(DATA_r{Group});
            data = tfrAllAsd;
            tfrASD = cell2mat(tfr'); 
        end


                    tfrAllAsd = sumtfr./length(DATA_r{Group})-nans;
                    data = tfrAllAsd;
%%
data = tfrAllTd - tfrAllAsd;

%% Baselining the TFR by subtracting the TFR from the pre-stim period, at the average level (as in Hu and Iannetti, 2014). 
% done per one group each time on "data". 
clear data_bl_TD data_bl_ASD

blMethod = 1;
timeBefore = [128:256]; %time before the first Epoch; DO IT PER EACH FREQUENCY!

switch blMethod
    case 1      %  baseline normalization (and log to get db)
            normWindow = nanmean(tfrAllTd(:,timeBefore),2);
            normWindow = nanmean(tfrAllAsd(:,timeBefore),2);

        for i = 1:length(timeoi)
             %data_bl_TD(:,i) = tfrAllTd(:,i)./normWindow;
             data_bl_ASD(:,i) = tfrAllAsd(:,i)./normWindow;
        end
    
   % data_bl_groupsDiff = data_bl_TD-data_bl_ASD;
    
    
    
    case 2                    % Baseline on the group level
    
    
    for i = 1:length(freq)
        %data_bl_TD(i,:) = tfrAllTd(i,:)-nanmean(tfrAllTd(i,timeBefore));
        data_bl_ASD(i,:) = tfrAllAsd(i,:)-nanmean(tfrAllAsd(i,timeBefore));
    end
    
    data_bl_groupsDiff = data_bl_TD-data_bl_ASD;
    
    case 3                % Baseline on the individual level
        
        for k = 1:size(tfrTD,1)
            for m = 1:length(freq)
                tfrTD_bl(k,m,:) = tfrTD(k,m,:)-nanmean(tfrTD(k,m,timeBefore));
            end
        end
        
        
        for k = 1:size(tfrTD,1)
            for m = 1:length(freq)
                tfrASD_bl(k,m,:) = tfrASD(k,m,:)-nanmean(tfrASD(k,m,timeBefore));
            end
        end
end
    

%% Plot
%Data = data;
Data = db(data_bl_ASD);

N = 500;
timeReduction = [1:length(timeoi)]; %change according to the NaNs in the data matrix. Try to get rid of as many as possible. 
%timeReduction = [100:length(timeoi)-100];
Timeoi = timeoi(timeReduction);
%clear Data

%Data = data_bl_groupsDiff(:,timeReduction); %data_bl is baselined data (see VAMP_TFR);
  
    [n, m] = size(Data);
    [x,y] = meshgrid(Timeoi,freqoi); % low-res grid
    [x2,y2] = meshgrid(Timeoi(1):1/N/5:Timeoi(end),freqoi(1):.01:freqoi(end));  %high-res grid
    dataInterp = interp2(x,y,Data, x2,y2, 'linear'); %interpolate up


    figure;
    %subplot(4,4,[1 2 2 8])

    f = surf(x2,y2,dataInterp);
    f.EdgeColor = 'none';
    f.FaceColor = 'interp';
    f.FaceLighting = 'gouraud';
    set(gca,'ydir','normal')
    ylabel('Frequency (Hz)')
    xlabel('Time (Sec.)');
    colorbar;
    colormap jet;
    ax = gca;
    
    caxis(ax.CLim)
    %caxis([-20 20]);
    view(0,90)
    axis tight