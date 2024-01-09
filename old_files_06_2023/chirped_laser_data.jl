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

exp_255_60_25C_2=Experiment(
    "./Data/Monday11072022/25C spectrum monday.csv",
    "./Data/Monday11072022/255C TEC 25C TFBG_ALL.csv",
    "./Data/Monday11072022/255C TEC 25C TFBG POWER_ALL.csv",
    1547.265,
    1547.451,
    18.3,
    60,
    30,
    25.5,
    13.7,
    32
)

exp_272_60_25C_2=Experiment(
    "./Data/Monday11072022/25C spectrum monday.csv",
    "./Data/Monday11072022/272C TEC 25C TFBG_ALL.csv",
    "./Data/Monday11072022/272C TEC 25C TFBG POWER_ALL.csv",
    1547.419,
    1547.614,
    18.3,
    60,
    30,
    27.2,
    13.7,
    32
)

exp_257_60_25C=Experiment(
    "./Data/Monday11072022/retaped_grating/retaped_spectrum monday.csv",
    "./Data/Monday11072022/retaped_grating/257C TEC ROOM TFBG_ALL.csv",
    "./Data/Monday11072022/retaped_grating/257C TEC ROOM TFBG POWER_ALL.csv",
    1547.267,
    1547.438,
    18.3,
    60,
    30,
    25.7,
    13.7,
    32
)

exp_272_60_25C_3=Experiment(
    "./Data/Monday11072022/retaped_grating/retaped_spectrum monday.csv",
    "./Data/Monday11072022/retaped_grating/272C TEC ROOM TFBG_ALL.csv",
    "./Data/Monday11072022/retaped_grating/272 TEC ROOM TFBG POWER_ALL.csv",
    1547.381,#previously .381
    1547,
    18.3,
    60,
    30,
    27.2,
    13.7,
    32
)

exp_wc_272c_18c = Experiment(
    "./Data/Tuesday07122022/Spectrums/18C spectrum water.csv",
    "./Data/Tuesday07122022/272C TEC 18C WATER3_ALL.csv",
    "./Data/Tuesday07122022/272C TEC POWER_ALL.csv",
    1547.381,#from anritsu measpreviously .381
    1547.568,
    18.3,
    60,
    30,
    27.2,
    13.7,
    32,
)

#Friday measurement, with power reference at grating minimum #12,333mv
exp_272c_18c = Experiment(
    "./Data/Friday07152022/18C spectrum Friday15072022.csv",
    "./Data/Friday07152022/272C TEC 18C TFBG 20220715_ALL.csv",
    "./Data/Tuesday07122022/272C TEC POWER_ALL.csv",#should be the same power didn't measure on Friday
    1547.381,#from anritsu measpreviously .381
    1547.568,
    18.3,
    60,
    30,
    27.2,
    13.7,
    32,
)

#12.83-12.87. minimumpower
exp_20220718_01=Experiment(
    "./Data/20220718_DWDD/20220718_tfbg_air_2passhann.csv",
    "./Data/20220718_DWDD/tfbg-air-2718c-22-9ma_000_ALL.csv",
    "./Data/20220718_DWDD/tfbg-air-2718c-22-9ma-power_000_ALL.csv",#should be the same power didn't measure on Friday
    1547.381,#from anritsu measpreviously .381
    1547.568,
    18.3,
    60,
    30,
    27.2,
    13.7,
    32,
)

#12.0?mv checkreadme
exp_20220721=Experiment(
    "./Data/20220721_varied_T/TFBG_20220720_h2_polarized_adjusted_to_1547_18c.csv",
    "./Data/20220721_varied_T/S_data/18c repeat_000_ALL.csv",
    "./Data/20220721_varied_T/power reference_000_ALL.csv",
    1547.1,#?
    1547.272,
    0,
    45,#high
    2,#low
    23.54,
    0,
    0,
)

#37.25
exp_20220721_100hz=Experiment(
    "./Data/20220721_varied_T/TFBG_20220720_h2_polarized_adjusted_to_1547_18c.csv",
    "./Data/20220721_varied_T/100hz/18c 100hz_ALL.csv",
    "./Data/20220721_varied_T/100hz/power reference 100hz_ALL.csv",
    1547.1,#?
    1547.272,
    0,
    45,#high
    2,#low
    23.54,
    0,
    0,
)
experiment_list=[exp_20220718_01,exp_20220721, exp_20220721_100hz]#[exp_372_60,exp_273_60,exp_27_40,exp_255_60_25C,exp_272_60_25C,exp_255_60_25C_2,exp_272_60_25C_2,exp_257_60_25C,exp_272_60_25C_3]

end
