module ChirpedLaserData
export experiment_list

struct Experiment
    OSA_file::String
    tfbg_file::String
    system_risetime_file::String
    λ0::Float64
    λ_high::Float64
    ΔI::Float64
    amplitude::Float64
    offset::Float64
    temp::Float64
    i_low::Float64
    i_high::Float64
end


exp_372_60=Experiment(
    "./Data/OSA_data/TFBG1_4deg_20220705.csv",
    "./Data/20220706_varied_measurements/372C 2kHz 13_313mA_ALL.csv",
    "./Data/20220706_varied_measurements/372C 2kHz 13 313mA POWER.csv",
    1548.252,
    1548.451,
    18.3,
    60,
    30,
    37.2,
    13,
    31.3,
)

exp_273_60=Experiment(
    "./Data/OSA_data/TFBG1_4deg_20220705.csv",
    "./Data/20220706_varied_measurements/273C 2kHz 129 312mA_ALL.csv",
    "./Data/20220706_varied_measurements/273C 2kHz 129 312mA POWER.csv",
    1547.447,
    1547.639,
    18.3,
    60,
    30,
    27.3,
    12.9,
    31.2
)
# exp_362_60=Experiment(
#     "./Data/OSA_data/TFBG1_4deg_20220705.csv",
#     "./Data/20220706_varied_measurements/372C 2kHz 13_313mA_ALL.csv",
#     "./Data/20220706_varied_measurements/372C 2kHz 13 313mA POWER.csv",
#     1548.252,
#     1548.451,
#     18.3,
#     60,
#     30,
#     37.2,
#     13,
#     31.3,
# )
exp_27_40=Experiment(
    "./Data/OSA_data/TFBG1_4deg_20220705.csv",
    "./Data/20220706_varied_measurements/27C 2kHz 16 281mA_ALL.csv",
    "./Data/20220706_varied_measurements/27C 2kHz 16 281mA POWER.csv",
    1547.439,
    1547.56,
    12.1,
    40,
    30,
    27,
    16,
    28.1
)

exp_255_60_25C=Experiment(
    "./Data/friday080722/tfbg_spectrum_20220708.csv",
    "./Data/friday080722/255C TEC 25C TFBG_ALL.csv",
    "./Data/friday080722/255C TEC 25C TFBG POWER_ALL.csv",
    1547.279,
    1547.474,
    18.3,
    60,
    30,
    25.5,
    13.7,
    32
)

exp_272_60_25C=Experiment(
    "./Data/friday080722/tfbg_spectrum_20220708.csv",
    "./Data/friday080722/272C TEC 25C TFBG_ALL.csv",
    "./Data/friday080722/272C TEC 25C TFBG POWER_ALL.csv",
    1547.419,
    1547.614,
    18.3,
    60,
    30,
    27.2,
    13.7,
    32
)
experiment_list=[exp_372_60,exp_273_60,exp_27_40,exp_255_60_25C,exp_272_60_25C]

end
