function test_rpeaks(alg, lang, test_no)

switch (alg)
    case 'Pan-Tompkins'
        ref_rpeaks_filename = sprintf('../%d/PanTompkinsOutput.csv', test_no);
        res_rpeaks_filename = sprintf('../%d/PanTompkinsResults%s.csv', test_no, lang);
    case 'Hilbert'
        ref_rpeaks_filename = sprintf('../%d/HilbertOutput.csv', test_no);
        res_rpeaks_filename = sprintf('../%d/HilbertResults%s.csv', test_no, lang);
    case 'EMD'
        ref_rpeaks_filename = sprintf('../%d/EMDOutput.csv', test_no);
        res_rpeaks_filename = sprintf('../%d/EMDResults%s.csv', test_no, lang);
end

ecg_filename = sprintf('../%d/Input.csv', test_no);

ecg = csvread(ecg_filename);
f = 360.0;
T = 1 / f;
t_vec = 0:(size(ecg) - 1);
t_vec = t_vec * T;

ref_rpeaks = csvread(ref_rpeaks_filename)';
ref_rpeaks = ref_rpeaks(1:(end-1));
res_rpeaks = csvread(res_rpeaks_filename);

size(setdiff(ref_rpeaks, res_rpeaks), 1)
% setdiff(ref_rpeaks, res_rpeaks)
size(setdiff(res_rpeaks, ref_rpeaks), 1)
% setdiff(res_rpeaks, ref_rpeaks)
size(ref_rpeaks, 1)
size(unique(res_rpeaks), 1)

% ref_rpeaks = ref_rpeaks + ones(size(ref_rpeaks));
% res_rpeaks = res_rpeaks + ones(size(res_rpeaks));

ref_rpeaks_val = ecg(ref_rpeaks + ones(size(ref_rpeaks)));
res_rpeaks_val = ecg(res_rpeaks + ones(size(res_rpeaks)));

figure('units','normalized','outerposition',[0 0 1 1])
plot(t_vec, ecg, 'r');
hold on;
grid on;
plot(ref_rpeaks * T, ref_rpeaks_val, 'go', 'LineWidth', 2);
plot(res_rpeaks * T, res_rpeaks_val, 'b+', 'MarkerSize', 10, 'LineWidth', 1);

title(sprintf('MIT-BIH Arythmia Database - sygna³ nr %d - algorytm %s', test_no, alg));
xlabel('czas [s]');
ylabel('napiêcie [mV]');
legend('EKG', 'C++', 'Python');

% window = 25000;
% ref_rpeaks_indices = find(ref_rpeaks < window);
% ref_rpeaks_window = ref_rpeaks(ref_rpeaks_indices);
% ref_rpeaks_window_val = ecg(ref_rpeaks_window);
% res_rpeaks_indices = find(res_rpeaks < window);
% res_rpeaks_window = res_rpeaks(res_rpeaks_indices);
% res_rpeaks_window_val = ecg(res_rpeaks_window);

% figure
% plot(t_vec(1:window), ecg(1:window), 'r');
% hold on;
% grid on;
% plot(ref_rpeaks_window * T, ref_rpeaks_window_val, 'go');
% plot(res_rpeaks_window * T, res_rpeaks_window_val, 'b+');
% 
% title(sprintf('MIT-BIH Arythmia Database - sygna³ nr %d - algorytm %s', test_no, alg));
% xlabel('czas [s]');
% ylabel('napiêcie [mV]');
% legend('EKG', 'C++', lang);

end

