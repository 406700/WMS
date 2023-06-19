module ChirpedLaser_S_Data
export s_experiment_list


struct S_Experiment
        path::String
        #chirp_param::Dict
end


exp_20220722=S_Experiment(
"Data/20220721_varied_T/S_data/"
)
## series of experiments with polarization controll, varied T and varied chirp.
exp_20220722_45mv_2mv=S_Experiment(
"Data/20220722/45mv_2mv/"
)

exp_20220722_45mv_10mv=S_Experiment(
"Data/20220722/45mv_10mv/"
)

exp_20220722_45mv_35mv=S_Experiment(
"Data/20220722/45mv_35mv/"
)

s_experiment_list=[exp_20220722,exp_20220722_45mv_2mv,exp_20220722_45mv_10mv,exp_20220722_45mv_35mv]

end
