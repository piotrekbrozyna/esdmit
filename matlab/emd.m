function [output, executionTime] = EMD(input)
% clear all;
% input = csvread('../104/Input.csv');
tic
samplingFrequency = 360;
input = transpose(input(:));
inputLength = length(input);

maxParts = ceil(inputLength/1000);
imf = [];
for j =1:3
    %tempImf = [];
%     for i=0:1:maxParts-1
%         beginPosition = i*1000+1;
%         endPosition = (i+1)*1000;
%         if(endPosition>inputLength)
%             endPosition = inputLength;
%         end
        %x1=input(beginPosition:endPosition);
        x1 = input;
        sd = Inf;
        while (sd > 0.3) 
            s1 = getspline(x1);
            s2 = -getspline(-x1);
            x2 = x1-(s1+s2)/2;
            sd = sum((x1-x2).^2)/sum(x1.^2);
            x1 = x2;
        end
        %tempImf = [tempImf x1];

    %end
    input = input - x1;
    imf{j} = x1;
end

imf = imf{1}+imf{2}+imf{3};

deartiveKernel = [ -0.125; -0.250; 0; 0.250; 0.125 ];
output = [];
%% derative signal
derativedSignal = conv(imf, deartiveKernel);
%% hilbert
hilbertSignal = abs(hilbert(derativedSignal(3:end-2)));

hilbertSignal = sqrt(hilbertSignal.^2+derativedSignal(3:end-2).^2);

k = 2*60*samplingFrequency;
size = length(imf);
for i=1:k:size
   beginElem = i;
   endElem = i + k;
   if(endElem > size)
       endElem = size;
   end

   tempInput = imf(beginElem:endElem);
   tempHilbertSignal = hilbertSignal(beginElem:endElem);

   maxValueEnd = 5*samplingFrequency;
   if maxValueEnd > length(tempHilbertSignal)
       maxValueEnd = length(tempHilbertSignal);
   end
   maxValue = max(tempHilbertSignal(1:maxValueEnd));
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

executionTime=toc;
end


% FUNCTIONS

function s = getspline(x)

N = length(x);
[peaks, pos] = findpeaks(x);
s = spline([1 pos N] ,[x(1) peaks x(N)],1:N);
end