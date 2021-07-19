function cv2 = cv2(array)

%only for arrays with one raw
if size(array,1) > 1
    error('input array dimension must be 1xN')
end
    
if size(array,2)>2
    for i = 1:length(array)-1
        cv2(i) = 2*abs(array(1,i)-array(1,i+1))/(array(1,i)+array(1,i+1));
    end
else
    cv2 = NaN;    
end