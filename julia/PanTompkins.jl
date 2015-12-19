function processPanTompkins(ecgSignal, samplingFrequency=360.0, tIntegrationWindow=0.15, tThRadius=0.1)
	const DIFF_KERNEL = [-0.125, -0.25, 0, 0.25, 0.125]
	const SIGNAL_LOW_LIMIT = 0.0001
	ecg_signal_size = size(ecgSignal)
	nThRadius = convert(Int, floor(samplingFrequency * tThRadius))
	fiducial_marks = Int[]
	rPeaks = Int[]
	rr = 0.0
	rrMean = 0.0

	differentiated_signal = conv(ecgSignal, DIFF_KERNEL)[3:end-2]
	squared_signal = differentiated_signal .^ 2

	n = convert(Int, floor(tIntegrationWindow * samplingFrequency))
	int_kernel_val = 1.0 / n
	int_kernel = fill(1.0, n)
	size_diff = n - 1.0
	int_start = convert(Int, floor(size_diff / 2.0))
	int_end = convert(Int, ceil(size_diff / 2.0))
	integrated_signal = conv(squared_signal, int_kernel)[int_start:end-int_end]
	println(size(integrated_signal))

	for i = 1:size(integrated_signal, 1)
		integrated_signal[i] = (integrated_signal[i] < SIGNAL_LOW_LIMIT) ? 0.0 : integrated_signal[i]
	end

	return rPeaks
end

refDir = dirname(pwd())
refNumStr = "100"
refInFileDir = joinpath(refDir, refNumStr, "Input.csv")
resOutFileDir = joinpath(refDir, refNumStr, "PanTompkinsResultsJulia.csv")

refInVector = readcsv(refInFileDir, Float64)
# convert 2D to 1D array and delete last element (trailing comma in CSV file...)
refInVector = float(reshape(refInVector, size(refInVector, 2))[1:end-1])

resultVector = processPanTompkins(refInVector)

writecsv(resOutFileDir, resultVector)