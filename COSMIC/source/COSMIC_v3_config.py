"""
* COSMIC - CONFIG
* Version 3
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""

import numpy as np
import amuse.units.units as u
import datetime

def params_dic():
    params_list = ["N_stars", "timestep", "output_times", "R0_king", "W0_king", "metallicity_stellar", "M_min_salpeter", "M_max_salpeter", "a_salpeter", "binary_fraction", "mean_period", "std_period", "X_coefficient_OBA", "X_coefficient_GKM", "out_directory", "in_directory", "config", "filename", "workers_stellar", "workers_gravity", "format_type"]
    params = {i:0 for i in params_list}
    return params, params_list

def verbatim(params):
    print("\033[94m"+"Info: Simulation time of "+"\033[0m"+"{:.2f} Myr".format(np.max(params["output_times"].value_in(u.Myr))))
    print("\033[94m"+"Info: Simulation of "+"\033[0m"+"{:d} stars".format(params["N_stars"]))
    print("\033[94m"+"Info: Simulation will be saved with the prefix "+"\033[0m"+"{}".format(params["filename"]))
    print("\033[94m"+"Info: Gravity: "+"\033[0m"+"{:d} workers".format(params["workers_gravity"])+"\033[94m"+", Stellar: "+"\033[0m"+"{:d} workers".format(params["workers_stellar"]))
    return 0

def default():
    params, params_list = params_dic()

    date = datetime.datetime.now()
    date = "{:04d}-{:02d}-{:02d}_{:02d}h{:02d}".format(date.year, date.month, date.day, date.hour, date.minute)

    params["N_stars"] = np.int32(10) # Number of particle
    params["timestep"] = np.float32(1)|u.Myr # Time increment
    params["output_times"] =  np.float32([1,2,5,10,20,50,100,200,500,1000,2000,5000])|u.Myr # Output intants
    params["R0_king"] = np.float32(1)|u.pc # Cluster virial radius
    params["W0_king"] = np.float32(9) # Cluster King's parameter
    params["metallicity_stellar"] = np.float32(0.02) # Metallicity
    params["M_min_salpeter"] = np.float32(0.1)|u.MSun # Minimum mass for Salpeter's law
    params["M_max_salpeter"] = np.float32(125)|u.MSun # Maximum mass for Salpeter's law
    params["a_salpeter"] = np.float32(-2.35) # Salpeter's law coefficient
    params["binary_fraction"] = np.float32(0.13) # Fraction of binary stars
    params["mean_period"] = np.float32(4.54) # Mean orbital period of binary systems (log10(T/days))
    params["std_period"] = np.float32(2.4) # Orbital period deviation of binary systems
    params["X_coefficient_OBA"] = np.float32(1.4e-7) # X-ray luminosity coefficient for O, B and A type stars
    params["X_coefficient_GKM"] = np.float32(1.4e27) # X-ray luminosity coefficient for G, K and M type stars
    params["out_directory"] = "./output/"
    params["in_directory"] = "./input/"
    params["config"] = "default.cfg"
    params["filename"] = "COSMIC_{}".format(date)
    params["workers_stellar"] = np.int32(1)
    params["workers_gravity"] = np.int32(1)
    params["format_type"] = "csv"
    return params, params_list

def change_params(params, cfg_param, cfg_values):
    if cfg_param == "":
        return params
    elif cfg_param == "N_stars": params["N_stars"] = np.int32(cfg_values[0])
    elif cfg_param == "timestep": params["timestep"] = np.float32(cfg_values[0])|u.Myr
    elif cfg_param == "output_times": params["output_times"] = np.float32(cfg_values)|u.Myr
    elif cfg_param == "R0_king": params["R0_king"] = np.float32(cfg_values[0])|u.pc
    elif cfg_param == "W0_king": params["W0_king"] = np.float32(cfg_values[0])
    elif cfg_param == "metallicity_stellar": params["metallicity_stellar"] = np.float32(cfg_values[0])
    elif cfg_param == "M_min_salpeter": params["M_min_salpeter"] = np.float32(cfg_values[0])|u.MSun
    elif cfg_param == "M_max_salpeter": params["M_max_salpeter"] = np.float32(cfg_values[0])|u.MSun
    elif cfg_param == "a_salpeter": params["a_salpeter"] = np.float32(cfg_values[0])
    elif cfg_param == "binary_fraction": params["binary_fraction"] = np.float32(cfg_values[0])
    elif cfg_param == "mean_period": params["mean_period"] = np.float32(cfg_values[0])
    elif cfg_param == "std_period": params["std_period"] = np.float32(cfg_values[0])
    elif cfg_param == "X_coefficient_OBA": params["X_coefficient_OBA"] = np.float32(cfg_values[0])
    elif cfg_param == "X_coefficient_GKM": params["X_coefficient_GKM"] = np.float32(cfg_values[0])
    elif cfg_param == "out_directory": params["out_directory"] = cfg_values[0].replace("\n", "")
    elif cfg_param == "in_directory": params["in_directory"] = cfg_values[0].replace("\n", "")
    elif cfg_param == "config": params["config"] = cfg_values[0].replace("\n", "")
    elif cfg_param == "filename": params["filename"] = cfg_values[0].replace("\n", "")
    elif cfg_param == "workers_stellar": params["workers_stellar"] = np.int32(cfg_values[0])
    elif cfg_param == "workers_gravity": params["workers_gravity"] = np.int32(cfg_values[0])
    elif cfg_param == "format_type": params["format_type"] = cfg_values[0].replace("\n", "")
    else:
        print("\033[93m"+"Warning: cannot understand the parameter in config: "+"\033[0m" + cfg_param)
    return params
            

def from_cfg_file(params):
    cfg_path = params["in_directory"] + params["config"]
    with open(cfg_path, "r") as cfg_file:
        for cfg_line in cfg_file.readlines():
            cfg_split = cfg_line.split(" ")
            if len(cfg_split) == 1:
                print("\033[93m"+"Warning: cannot understand the line in config: "+"\033[0m" + cfg_line)
                continue
            cfg_param = cfg_split[0]
            cfg_values = cfg_split[1:]
            params = change_params(params, cfg_param, cfg_values) 
    verbatim(params)
    return params

