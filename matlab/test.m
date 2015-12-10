clear all;

A = dir('..');
result = [];
for k=1:1:length(A)
   clearvars -except A k result
   pathInput = strcat('../', A(k).name, '/PanTompkinsInput.csv');
   pathOutput = strcat('../', A(k).name, '/PanTompkinsOutput.csv');
   if(A(k).isdir && exist(pathInput, 'file') && exist(pathOutput, 'file'))
        A(k).name
        input = csvread(pathInput);
        samplingFrequencis=360;
        PanTompkins;
        test1 = csvread(pathOutput);
        temp = output-1;
%         result(k).lengthDiff = length(test1)-length(temp);
%         result(k).inRef = setdiff(test1,temp);
%         result(k).inTest = setdiff(temp,test1);
%         %;
%         %
%         %setdiff(temp,test1)
%         result = abs(temp-test1);
%         if((length(test1)-length(temp))==0)
%             max(abs(temp-test1))     
%         else
%             pathInput
%             'dupa'
             setdiff(test1,temp)
             setdiff(temp,test1)
%         end
   end
end
