% Calculating dprime for binary detection test results.
% Input vars: 4*N matrix of  absolute numbers. N is number of observations/participatns.
% 4 rows of: hit, miss, false alarms, correct rejections, in this order.
% Output vars = dprime, calculated as z[hit rate]-z[false alarm rate] for
% each observation/participant. Z being the inverse of the standard
% normal cumulative distribution function (cdf)
% TPR (aka sensitivity or hit rate) = Hit/(Hit+miss)
% FPR (aka false alarm rate or fall-out) = (false alarm)/(false alarm+correct Rejection)
%
% Shlomit Beker May 2020

function [dprime,TPR,FPR] = dPrimeCalc(fullMat);

    
    Hit = fullMat(1,:);
    Miss = fullMat(2,:);
    
    TPR = Hit./(Hit+Miss); %True positive rate
    TPR(find(TPR==1))=0.99;
    
  if size(fullMat,1) == 4
    
    falseAlarm = fullMat(3,:);
    corRej = fullMat(4,:);
    
    FPR = falseAlarm./(falseAlarm+corRej); %False positive rate
    
    FPR(find(FPR==0))=0.01;
         
    for i = 1:length(TPR)
        
        dprime(i) = norminv(TPR(i))-norminv(FPR(i)); % calculating Z-scores of the difference
        
    end
  else
      
      for i = 1:length(TPR)
        
        dprime(i,1) = TPR(i); % calculating Z-scores of the difference
        
      end
      
end

end