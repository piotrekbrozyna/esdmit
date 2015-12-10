function test_rpeaks(alg, test_no)

switch (alg)
    case 'pan-tompkins'
        ecg_filename = sprintf('../%d/PanTompkinsInput.csv', test_no);
        ref_rpeaks_filename = sprintf('../%d/PanTompkinsOutput.csv', test_no);
        res_rpeaks_filename = sprintf('../%d/PanTompkinsResultsPython.csv', test_no);
    case 'hilbert'
        ecg_filename = sprintf('../%d/HilbertInput.csv', test_no);
        ref_rpeaks_filename = sprintf('../%d/HilbertOutput.csv', test_no);
        res_rpeaks_filename = sprintf('../%d/HilbertResultsPython.csv', test_no);
    case 'emd'
        ecg_filename = sprintf('../%d/PanTompkinsInput.csv', test_no);
        ref_rpeaks_filename = sprintf('../%d/PanTompkinsOutput.csv', test_no);
        res_rpeaks_filename = sprintf('../%d/EMDResultsPython.csv', test_no);
end

ecg = csvread(ecg_filename);
ref_rpeaks = csvread(ref_rpeaks_filename);
res_rpeaks = csvread(res_rpeaks_filename);

ref_rpeaks = ref_rpeaks + ones(size(ref_rpeaks));
res_rpeaks = res_rpeaks + ones(size(res_rpeaks));

ref_rpeaks_val = ecg(ref_rpeaks);
res_rpeaks_val = ecg(res_rpeaks);

plot(ecg, 'r');
hold on;
plot(ref_rpeaks, ref_rpeaks_val, 'bo');
plot(res_rpeaks, res_rpeaks_val, 'm+');

title(sprintf('MIT-BIH Arythmia Database - signal no. %d', test_no));
xlabel('sample');
ylabel('voltage');
legend('ECG', 'reference', 'python');

end

