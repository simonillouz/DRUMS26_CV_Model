function E = Elastance(t,EM, Em, TC, TR)
%Input the following variables:
%t = time from 0 to T
%T = heart rate
%TC = fraction of time to go from min to max elasticity
%TR = fraction of time to go from max to min elasticity
%Em  = minimum elasticity
%EM  = maximum elasticity

% if else statement for left heart
if t<=TC
   E = Em + ((EM-Em)/2)*(1-cos(pi*t/TC)); 
elseif t <=TR+TC    
   E = Em + ((EM-Em)/2)*(cos((t-TC)*pi/TR) + 1);     
else
   E = Em;
end