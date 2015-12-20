const DIFF_KERNEL = [-0.125; -0.25; 0; 0.25; 0.125]
const SIGNAL_LOW_LIMIT = 0.0001
const RR_SIZE = 8
const RR_LOW_LIMIT_FACTOR = 0.92
const RR_HIGH_LIMIT_FACTOR = 1.16
const RR_MISSED_LIMIT_FACTOR = 1.66

const SIGNAL_LOW_LIMIT = 0.0001

# 360 ms
const T_WAVE_DETECTION_RR_HIGH_LIMIT = 0.36

# thresholding factors
const SECOND_THRESHOLD_SPKF_PEAKF_FACTOR = 0.25
const SECOND_THRESHOLD_SPKF_SPKF_FACTOR = 0.75
const SECOND_THRESHOLD_SPKI_PEAKI_FACTOR = 0.25
const SECOND_THRESHOLD_SPKI_SPKI_FACTOR = 0.75

const SPKI_UPDATE_FACTOR_1 = 0.125
const SPKI_UPDATE_FACTOR_2 = 0.875
const NPKI_UPDATE_FACTOR_1 = 0.125
const NPKI_UPDATE_FACTOR_2 = 0.875
const SPKF_UPDATE_FACTOR_1 = 0.125
const SPKF_UPDATE_FACTOR_2 = 0.875
const NPKF_UPDATE_FACTOR_1 = 0.125
const NPKF_UPDATE_FACTOR_2 = 0.875
const THRESHOLDI_UPDATE_FACTOR_1 = 0.25
const THRESHOLDI_UPDATE_FACTOR_2 = 0.5
const THRESHOLDF_UPDATE_FACTOR_1 = 0.25
const THRESHOLDF_UPDATE_FACTOR_2 = 0.5

function processPanTompkins(ecg_signal, sampling_frequency=360.0, t_int_window=0.15, t_radius=0.1)

	ecg_signal_size = size(ecg_signal)
	n_radius = convert(Int, floor(sampling_frequency * t_radius))
	fiducial_marks = Int[]
	r_peaks = Int[]
	rr_mean = 0.0
	rr = Array{Int}(RR_SIZE)

	differentiated_signal = conv(ecg_signal, DIFF_KERNEL)[3:end-2]
	squared_signal = differentiated_signal .^ 2

	n = convert(Int, floor(t_int_window * sampling_frequency))
	int_kernel_val = 1.0 / n
	int_kernel = fill(1.0, n)
	size_diff = n - 1.0
	int_start = convert(Int, floor(size_diff / 2.0))
	int_end = convert(Int, ceil(size_diff / 2.0))
	integrated_signal = conv(squared_signal, int_kernel)[int_start:end-int_end]

	for i = 1:size(integrated_signal, 1)
		integrated_signal[i] = (integrated_signal[i] < SIGNAL_LOW_LIMIT) ? 0.0 : integrated_signal[i]
	end

	# findFiducialMarks
	gradient_signal = [integrated_signal[1], diff(integrated_signal)]
	for i = 3:size(gradient_signal, 1)
		if gradient_signal[i-1] > 0 && gradient_signal[i] <= 0
			if size(fiducial_marks, 1) > 0
				if (i - fiducial_marks[end]) >= ceil(0.2 * sampling_frequency)
					append!(fiducial_marks, [i])
				elseif integrated_signal[i] >= integrated_signal[fiducial_marks[end]]
					fiducial_marks[end] = i
				end
			else
				append!(fiducial_marks, [i])
			end
		end
	end

	# initializeThresholds
	end_index = convert(Int, 2 * sampling_frequency)
	spki = maximum(integrated_signal[1:end_index]) / 3.0
	npki = sum(integrated_signal[1:end_index]) / (end_index * 2.0)
	thresholdi1 = spki
	thresholdi2 = npki
	spkf = maximum(ecg_signal[1:end_index]) / 3.0
	npkf = sum(ecg_signal[1:end_index]) / (end_index * 2.0)
	thresholdf1 = spkf
	thresholdf2 = npkf

	#threshold
	skip = false
	for mark in fiducial_marks

		pos = find_maximum(ecg_signal, mark, n_radius)

		if size(r_peaks, 1) >= RR_SIZE
			rr = diff(r_peaks[(RR_SIZE-1):end]) # calculateAvarageRR
			rr_mean = sum(rr) / RR_SIZE

			if is_arythmia(rr, rr_mean)
				thresholdf1 *= 0.5
				thresholdf2 *= 0.5
			end

			if possible_peak_omission(mark, r_peaks, rr_mean)
				# searchBack
				window_samples = convert(Int, floor(0.2 * sampling_frequency))

				if (mark - window_samples) < (r_peaks[end] + window_samples)
					mark_temp = indmax(integrated_signal[(mark - window_samples):(r_peaks[end] + window_samples)])
				else
					mark_temp = indmax(integrated_signal[(r_peaks[end] + window_samples):(mark - window_samples)])
				end

				pos_temp = find_maximum(ecg_signal, mark_temp, n_radius)

				if integrated_signal[mark_temp] > thresholdi2
					if ecg_signal[pos_temp] > thresholdf2
						append!(r_peaks, [pos_temp])
						spkf = SECOND_THRESHOLD_SPKF_PEAKF_FACTOR * ecg_signal[pos_temp] + SECOND_THRESHOLD_SPKF_SPKF_FACTOR * spkf
					end
					spki = SECOND_THRESHOLD_SPKI_PEAKI_FACTOR * integrated_signal[mark_temp] + SECOND_THRESHOLD_SPKI_SPKI_FACTOR * spki
				end
			end
		end

		if integrated_signal[mark] >= thresholdi1
			if size(r_peaks, 1) > 2 && possible_t_wave(r_peaks, mark, sampling_frequency)
				if compare_peaks_slopes()
					npki = NPKI_UPDATE_FACTOR_1 * integrated_signal[mark] + NPKI_UPDATE_FACTOR_2 * npki
        			npkf = NPKF_UPDATE_FACTOR_1 * ecg_signal[pos] + NPKF_UPDATE_FACTOR_2 * npkf
        			skip = true
        		else
        			skip = false
        		end
			end

			if !skip
				if ecg_signal[pos] > thresholdf1
					append!(r_peaks, [pos])
					spkf = SPKF_UPDATE_FACTOR_1 * ecg_signal[pos] + SPKF_UPDATE_FACTOR_2 * spkf
				end
				spki = SPKI_UPDATE_FACTOR_1 * integrated_signal[mark] + SPKI_UPDATE_FACTOR_2 * spki
			elseif (thresholdi2 <= integrated_signal[mark] && integrated_signal[mark] < thresholdi1) || (integrated_signal[mark] < thresholdi2)
				spkf = SPKF_UPDATE_FACTOR_1 * ecg_signal[pos] + SPKF_UPDATE_FACTOR_2 * spkf
        		spki = SPKI_UPDATE_FACTOR_1 * integrated_signal[mark] + SPKI_UPDATE_FACTOR_2 * spki
        	end
		end

		thresholdi1 = npki + THRESHOLDI_UPDATE_FACTOR_1 * (spki - npki)
        thresholdi2 = THRESHOLDI_UPDATE_FACTOR_2 * thresholdi1
        thresholdf1 = npkf + THRESHOLDF_UPDATE_FACTOR_1 * (spkf - npkf)
        thresholdf2 = THRESHOLDF_UPDATE_FACTOR_2 * thresholdf1
        skip = false
	end

	return r_peaks
end

function find_maximum(signal, index, radius)
	signal_size = size(signal, 1)
	start_index = (index - radius) < 1 ? 1 : (index - radius)
	if index + radius > signal_size
		return start_index + indmax(signal[start_index:end])
	else
		return start_index + indmax(signal[start_index:(index + radius)])
	end
end

function is_arythmia(rr, rr_mean)
	return (rr[end] <= RR_LOW_LIMIT_FACTOR * rr_mean) || (rr[end] >= RR_HIGH_LIMIT_FACTOR * rr_mean)
end

function possible_peak_omission(mark, r_peaks, rr_mean)
	return (mark - r_peaks[end]) >= (RR_MISSED_LIMIT_FACTOR * rr_mean)
end

function possible_t_wave(r_peaks, mark, sampling_frequency)
	return (mark - r_peaks[end]) <= ceil(T_WAVE_DETECTION_RR_HIGH_LIMIT * sampling_frequency)
end

ref_dir = dirname(pwd())
ref_num_str = "100"
ref_in_file_dir = joinpath(ref_dir, ref_num_str, "Input.csv")
ref_out_file_dir = joinpath(ref_dir, ref_num_str, "PanTompkinsResultsJulia.csv")

ref_in_vector = readcsv(ref_in_file_dir, Float64)
# convert 2D to 1D array and delete last element (trailing comma in CSV file...)
ref_in_vector = float(reshape(ref_in_vector, size(ref_in_vector, 2))[1:end-1])

result_vector = processPanTompkins(ref_in_vector)

println(result_vector)

writecsv(ref_out_file_dir, result_vector)