function [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(array)


Mean            = mean(array);
SD              = std(array);
SEM             = Mean/sqrt(length(array));
Min             = min(array);
Med             = median(array);
Max             = max(array);
cv              = SD/Mean;
FanoFactor      = SD^2/Mean;


if length(array)>2
    
    %cv2
    for i = 1:length(array)-1
        cv2_data(i) = 2*abs(array(i)-array(i+1))/(array(i)+array(i+1));
    end
    
    cv2             = mean(cv2_data);
    
else
    cv2         = NaN;
    
    
end