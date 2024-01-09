using DelimitedFiles
using Plots
plotlyjs()

cd("./Data/OSA_data/Anritsu/")
files=readdir(".")

for i in files[]
    file1=readdlm(i,',',skipstart=29)
    a=plot!(file1[:,1],file1[:,2],ylabel="Intensity (mw)",xlabel="wavelength (nm)",title=i[1:end-4])
    gui(a)
    #savefig(i[1:end-4])
end

savefig("superposition_peaks_and_chirp")
