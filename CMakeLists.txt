set(R_PEAKS_source_files
#  test.cpp
#  another_file.cpp
r_peaks_hard_coded_data.cpp
R_PEAKS_PanTompkins.cpp
R_PEAKS_Hilbert.cpp
R_PEAKS_EMD.cpp
R_PEAKS_Tools.cpp
)
add_library(R_PEAKS ${R_PEAKS_source_files})

#qt5_use_modules(R_PEAKS Widgets Gui Core)
# double check if the module names are correct
