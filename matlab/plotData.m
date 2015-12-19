function plotData( ecg, ref_rpeaks, res_rpeaks, test_no)

ref_rpeaks_val = ecg(ref_rpeaks);
res_rpeaks_val = ecg(res_rpeaks);
figure;
plot(ecg, 'r');
hold on;
plot(ref_rpeaks, ref_rpeaks_val, 'bo');
plot(res_rpeaks, res_rpeaks_val, 'm+');

title(sprintf('MIT-BIH Arythmia Database - signal no. %d', test_no));
xlabel('sample');
ylabel('voltage');
legend('ECG', 'reference', 'matlab');

end

