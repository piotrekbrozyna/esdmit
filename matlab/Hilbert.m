%% initialize varaible
clear all;
input = csvread('../100/Input.csv');
samplingFrequency = 360;
deartiveKernel = [ -0.125; -0.250; 0; 0.250; 0.125 ];

lengthOfIntegrationWindow = floor(samplingFrequency*0.15);
integartionKernel = ones(lengthOfIntegrationWindow, 1)/lengthOfIntegrationWindow;
gradientKernel = [-1,1];
nThRadius = ceil(0.1*samplingFrequency);
skip = false;
pos = 0;
posTemp = 0;
fidualMarkTemp = 0;
output = [];
%% derative signal
derativedSignal = conv(input, deartiveKernel);
%% hilbert
hilbertSignal = abs(hilbert(derativedSignal(2:end-2)));

k = 5*60*samplingFrequency;
size = length(input);

for i=1:k:size
   beginElem = i;
   endElem = i + k;
   if(endElem > size)
       endElem = size;
   end
   
   tempInput = input(beginElem:endElem);
   tempHilbertSignal = hilbertSignal(beginElem:endElem);

   maxValue = max(tempHilbertSignal(1:5*samplingFrequency));
   threshold = 0.8*maxValue;
   
   tempOutput = 1;
   
   v = floor(0.7*samplingFrequency);
   
   RR = ones(8,1)*v;
   tmpECG = zeros(5,1);
   
   for j=1:1:length(tempInput)
       dt = j - tempOutput(end);
       mR = mean(RR);

       if (tempInput(j) > threshold)
            if (dt <= 0.2*samplingFrequency | dt <= 0.55*mR)
                if(tempInput(j) > tempInput(tempOutput(end)))    
                    tempOutput(end) = j;
                    tmpECG(mod(length(tempOutput),5)+1) = tempInput(tempOutput(end));
                end

                if (length(tempOutput) > 1)
                	RR(mod(length(tempOutput),8)+1) = tempOutput(end) - tempOutput(end-1);
                end
            else
                tempOutput = [tempOutput j];
                tmpECG(mod(length(tempOutput),5)+1) = tempInput(tempOutput(end));
                RR(mod(length(tempOutput),8)+1) = tempOutput(end) - tempOutput(end-1);        
            end
       end
       if(length(tempOutput) > 4)    
            threshold = 0.55/4*(sum(tmpECG)-max(tmpECG'));
       end
        
   end
   
   tempOutput = tempOutput(2:end);
   
   tempOutput = tempOutput +i-1;
   output = [output tempOutput];
   
    
end





