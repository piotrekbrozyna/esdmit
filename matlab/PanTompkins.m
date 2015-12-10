%% initialize varaible
%clear all;
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
%% square signal
squaredSignal = derativedSignal.^2;
%% integrate
integartedSignal = conv(squaredSignal, integartionKernel);
dl = length(integartedSignal)-length(input);
integratedSignal=integartedSignal((floor(dl/2)+1):(end-ceil(dl/2)),1);
for i=1:1:length(integratedSignal)
     if (integratedSignal(i)<=0.0001)
         integratedSignal(i)=0;
     end
end
%% coutn gradient
signalGradient = conv(integratedSignal, gradientKernel);
 
%% findFidualMark
fidualMark = [];
for i=3:1:length(signalGradient)-1
    if(signalGradient(i-1)>0 && signalGradient(i)<=0)
        if(~isempty(fidualMark))
            if (i-fidualMark(end) >= ceil(0.2*samplingFrequency))
                fidualMark = [fidualMark i]; 
            elseif(integratedSignal(i)>=integratedSignal(fidualMark(end)))
                fidualMark(end) = i;
            end
        else
            fidualMark = i; 
        end
    end
end
%% initialize var to threshold
spki = max(integratedSignal(1:2*samplingFrequency,1))/3;
npki = sum(integratedSignal(1:2*samplingFrequency,1))/(2*samplingFrequency*2);
thi1 = spki;
thi2 = npki;
spkf = max(input(1:2*samplingFrequency,1))/3;
npkf = sum(input(1:2*samplingFrequency,1))/(2*samplingFrequency*2);
thf1 = spkf;
thf2 = npkf;
%% threshold
for i=1:1:length(fidualMark)
    startMax = fidualMark(i) - nThRadius;
    endMax = fidualMark(i)+nThRadius;
    if(startMax < 1)
        startMax = 1;
    end
    if(endMax > length(input))
        endMax = length(input);
    end
    [val, pos] = max(input(startMax:endMax));
    pos = pos + startMax-1;
    if(length(output) > 8)
       for j=1:1:8
          RR(j) = output(end-8+j)-output(end-9+j); 
       end
    end
    if(exist('RR', 'var') && length(RR) == 8)
       mRR =  mean(RR);
       if(RR(end) <= 0.92*mRR || RR(end) >=1.16*mRR)
          thf1 = 0.5*thf1;
          thf2 = 0.5*thf2;
       end
    end
    
    if(exist('mRR', 'var') && mRR > 0)
       if(fidualMark(i) - output(end) >= 1.66*mRR)
           [val,fidualMarkTmp] = max(integratedSignal((output(end)+floor(0.2*samplingFrequency)):(fidualMark(i)-floor(0.2*samplingFrequency))));
           fidualMarkTmp = fidualMarkTmp +output(end)+floor(0.2*samplingFrequency);
            startMax = fidualMark(i) - nThRadius;
            endMax = fidualMark(i)+nThRadius;
            if(startMax < 1)
                startMax = 1;
            end
            if(endMax > length(input))
                endMax = length(input);
            end
            [val, posTmp] = max(input(startMax:endMax));
            posTmp = posTmp + startMax;
            if(integratedSignal(posTmp) > thi2)
               if(input(posTmp)>thf2)
                  output = [output posTmp];
                  spkf = 0.25*input(posTmp) + 0.75*spkf;
               end
               spki = 0.25*integratedSignal(posTmp)+0.75*spki;
            end
       end
    end
    
    if(integratedSignal(fidualMark(i)) >= thi1)
       if((length(output) > 2) && (fidualMark(i) - output(end) <= ceil(0.360*samplingFrequency)))
          slope1 = [-1 -2 0 2 1]*integratedSignal((fidualMark(i)-4):fidualMark(i),1);
          slope2 = [-1 -2 0 2 1]*integratedSignal((output(end)-4):output(end),1);
          if(abs(slope1) <= 0.5*abs(slope2))
              skip = true;
              npki = 0.125*integratedSignal(fidualMark(i))+0.875*npki;
              npkf = 0.125*input(pos)+0.875*npkf;
          else
              skip=false;
          end
          
       end
       
       if(skip == false)
          if(input(pos)>=thf1)
             
             output = [output pos];
             spkf = 0.125*input(pos)+0.875*spkf;
          end
          spki = 0.125*integratedSignal(fidualMark(i))+0.875*spki;
       end
    elseif((thi2 <= integratedSignal(fidualMark(i)) && integratedSignal(fidualMark(i)) < thi1) || integratedSignal(fidualMark(i)) < thi2)
        npki = 0.125*integratedSignal(fidualMark(i)) + 0.875*npki;
		npkf = 0.125*input(pos) + 0.875*npkf;
    end
    thi1 = npki + 0.25*(spki-npki);
    thi2 = 0.5*thi1;
    thf1 = npkf + 0.25*(spkf-npkf);
    thf2 = 0.5*thf1;
    
    skip = false;
    
end

