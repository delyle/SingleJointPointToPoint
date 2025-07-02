function lTarget = lengthShifting(t,tShift,lShift)
% Creates a vector of lengths lTarget of size equal to t at a constant
% value lSwitch(i) from t = tSwitch(i) to tSwitch(i+1)

lTarget = ones(size(t)); % initialize target muscle length vector
for i = 1:length(tShift)-1
    lTarget(tShift(i) <= t & t <= tShift(i+1)) = lShift(i);
end
