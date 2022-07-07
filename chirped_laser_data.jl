struct Experiment
    OSA_file::String
    tfbg_file::String
    system_risetime_file::String
    λ0::Float64
    ΔI::Float64
end

exp_27_281_16=Experiment(
    "./Data/OSA_data/TFBG1_4deg_20220705.csv",
    "./Data/20220706_varied_measurements/27C 2kHz 16 281mA_ALL.csv",
    "./Data/20220706_varied_measurements/27C 2kHz 16 281mA POWER.csv",
    1547.439,
    12.1,
)
