import DSP.hilbert
DIFF_KERNEL = [-0.125, -0.250, 0, 0.250, 0.125]
RR_SIZE = 8
function processHilbert(ecg_signal, sampling_frequency=360.0, tWindowSize=300.0)
  ecg_signal_size = size(ecg_signal,1)
  println(ecg_signal_size)
  nWindowSize = convert(Int, floor(tWindowSize * sampling_frequency))
	r_peaks = Int[]
  differentiated_signal = conv(ecg_signal, DIFF_KERNEL)[3:end-2]
  analyticalSignal = hilbert(differentiated_signal)
  envelope = abs(analyticalSignal)
  startIndex = 1
  while startIndex < ecg_signal_size
    endIndex = startIndex + nWindowSize
    endIndex = (ecg_signal_size < endIndex) ? ecg_signal_size : endIndex
    output = doThresholding(ecg_signal[startIndex:endIndex], envelope[startIndex:endIndex], sampling_frequency)+startIndex-1
    append!(r_peaks, output)
    startIndex += nWindowSize
  end
  return r_peaks
end

function doThresholding(ecgInterval, envelopeInterval, sampling_frequency)
  endIndex = convert(Int, floor(5 * sampling_frequency))
  signalSize = size(ecgInterval,1)
  endIndex = (signalSize<endIndex) ? signalSize : endIndex

  max = maximum(envelopeInterval[1:endIndex])
  threshold = 0.8 * max
  output = Int[1]
  rr = ones(RR_SIZE)*floor(0.7 *sampling_frequency)
  tmpEcg = zeros(5)

  for i=1:size(ecgInterval,1)
    data = ecgInterval[i]

    dt = i - output[end]
    rrMean = mean(rr)

    if data > threshold
      if dt <= 0.2 * sampling_frequency || dt <= 0.55 * rrMean
        if data > ecgInterval[output[end]]
          output[end] = i
          tmpEcg[(size(output,1) % 5+1)] = ecgInterval[output[end]]
        end
        if size(output,1) > 1
            rr[(size(output,1) % 8+1)] = output[end] - output[end-1]
        end
      else
        push!(output, i)
        tmpEcg[(size(output,1) % 5+1)] = ecgInterval[output[end]]
        rr[(size(output,1) % 8+1)] = output[end] - output[end-1]
      end
    end
    if (size(output,1) > 4)
      threshold = ((sum(tmpEcg) - indmax(tmpEcg)[1])/4)*0.55
    end
  end
  shift!(output)

  return output
end

ref_dir = "E:\\medyczny\\esdmit"
ref_num_str = "100"
ref_in_file_dir = joinpath(ref_dir, ref_num_str, "Input.csv")
ref_out_file_dir = joinpath(ref_dir, ref_num_str, "HilbertResultsJulia.csv")

ref_in_vector = readcsv(ref_in_file_dir, Float64)
# convert 2D to 1D array and delete last element (trailing comma in CSV file...)
ref_in_vector = float(reshape(ref_in_vector, size(ref_in_vector, 2))[1:end-1])
tic()
result_vector = processHilbert(ref_in_vector)
toc()

println(result_vector)

writecsv(ref_out_file_dir, result_vector)
