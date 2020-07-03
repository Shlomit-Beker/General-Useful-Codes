% this function rereferenced data to a different channel. 
% Inputs: preRefData: old data to reref
% refChan: the channel to ref to. 
% shlomit Beker 2018


function refDat =  rerefData(preRefData,refChan)

for i = 1:length(preRefData) 
    for j = 1:length(preRefData{i})
      for k = 1:length(preRefData{i}{j})  
        A{i}{j}{k} = ft_preproc_rereference(preRefData{i}{j}{k}, refChan);
        %tempDat{i}{j}{k} = preRefData{i}{j}{k};
        tempDat{i}{j}{k} = A{i}{j}{k};
      end
    end  
end


refDat = tempDat; 
