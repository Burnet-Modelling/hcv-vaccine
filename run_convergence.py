# -*- coding: utf-8 -*-
"""
Created on Fri May 22 13:53:37 2026

@author: farah.houdroge

Script to run convergence analysis on sample size (number of simulations).
"""
from hcv.generate_databooks import generate_databook
from hcv import utils as ut
from hcv import utils_plotting as ut_plt
import sys

dir = ut.get_project_root()

country = "VNM"
rand_seed = 250711
n_samples = 1000
results_folder = dir / "results"
cal_folder = results_folder / "calibration" / "y-factors"
savedir_scens = results_folder / "convergence" / country / f"{n_samples}_samples"
savedir_scens.mkdir(parents=True, exist_ok=True)
savedir_pp = results_folder / "convergence" / "results" / country / f"{n_samples}_samples"
savedir_pp.mkdir(parents=True, exist_ok=True)

ut.run_scenario_sampling(country, cal_folder, rand_seed=rand_seed, n_samples=n_samples, savedir=savedir_scens,scenarios_to_run=["counter_0", "scenario_0"])
ut.econ_eval(country, savedir_scens, savedir_pp, rand_seed=rand_seed, n_samples=n_samples,scenarios_to_run=["counter_0", "scenario_0"])

#%%

    # ref_countries = pd.read_excel(
    #     str(rootdir) + "/data/flat_datasheet.xlsx",
    #     sheet_name="Cost - YLL and productivity",
    # ).iloc[:, np.r_[1, 4]]
    # who_reg = list(pd.unique(ref_countries.WHO_reg)) + ["global", "top10"]
    # loop_folder = pathlib.Path(str(scens_folder) + "/agg_outputs/")
    # countries = [
    #     country.stem[:3] for country in loop_folder.iterdir() if country.is_file()
    # ]

    # scenarios = pd.read_excel(
    #     str(rootdir) + "/data/progbook_inputs.xlsx", sheet_name="scenarios"
    # )
    # scen_name = list(pd.unique(scenarios.scenario_name))

    # # Produce lists of countries within each region, this can be adjusted pretty easily now
    # ref_countries = ref_countries[ref_countries["ISO3"].isin(countries)].reset_index(
    #     drop=True
    # )
    # top_10_burden = [
    #     "PAK",
    #     "IND",
    #     "CHN",
    #     "RUS",
    #     "USA",
    #     "IDN",
    #     "NGA",
    #     "UKR",
    #     "UZB",
    #     "BGD",
    # ]

    # AFR, AMR, EMR, EUR, SEAR, WPR, TOP = [], [], [], [], [], [], []
    # for i in range(len(ref_countries)):
    #     if ref_countries.iloc[i, 1] == "AFR":
    #         AFR.append(ref_countries.iloc[i, 0])
    #     if ref_countries.iloc[i, 1] == "AMR":
    #         AMR.append(ref_countries.iloc[i, 0])
    #     if ref_countries.iloc[i, 1] == "EMR":
    #         EMR.append(ref_countries.iloc[i, 0])
    #     if ref_countries.iloc[i, 1] == "EUR":
    #         EUR.append(ref_countries.iloc[i, 0])
    #     if ref_countries.iloc[i, 1] == "SEAR":
    #         SEAR.append(ref_countries.iloc[i, 0])
    #     if ref_countries.iloc[i, 1] == "WPR":
    #         WPR.append(ref_countries.iloc[i, 0])
    #     if ref_countries.iloc[i, 0] in top_10_burden:
    #         TOP.append(ref_countries.iloc[i, 0])

    # # Produce empty dataframes to collate results (fillna means can do iterative sums)
    # epis_outs = ["HCV incidence", "HCC incidence", "HCV mortality", "DALYs"]

    # epi_agg = {}
    # for reg in who_reg:
    #     epi_agg[reg] = {}
    #     for scen in scen_name:
    #         epi_agg[reg][scen] = {}
    #         for out in epis_outs:
    #             epi_agg[reg][scen][out] = pd.DataFrame(
    #                 columns=["year", "central"] + [f"run_{i}" for i in range(n_samples)]
    #             )
    #             epi_agg[reg][scen][out].year = np.arange(2000.5, 2051.1, 1)
    #             epi_agg[reg][scen][out] = epi_agg[reg][scen][out].fillna(0)

    # # Economic aggregates
    # econ_outs = [
    #     "PWID Ab Tests",
    #     "PWID RNA Tests",
    #     "Prisoner Ab Tests",
    #     "Prisoner RNA Tests",
    #     "Gen Ab Tests",
    #     "Gen RNA Tests",
    #     "Total Ab Tests",
    #     "Total RNA Tests",
    #     "PWID Ab Positive",
    #     "PWID RNA Positive",
    #     "Prisoner Ab Positive",
    #     "Prisoner RNA Positive",
    #     "Gen Ab Positive",
    #     "Gen RNA Positive",
    #     "Total Ab Positive",
    #     "Total RNA Positive",
    #     "Total treatment",
    #     "PWID vaccinations",
    #     "Prisoner vaccinations",
    #     "Gen vaccinations",
    #     "Total vaccinations",
    #     "Diagnosis Cost",
    #     "Disease Management Cost",
    #     "Treatment Cost",
    #     "Morbidity Cost",
    #     "Mortality Cost",
    #     "Direct Cost",
    #     "Productivity Cost",
    #     "Total Cost",
    #     "Vaccine 5USD Cost",
    #     "Vaccine 10USD Cost",
    #     "Vaccine 15USD Cost",
    #     "Vaccine 20USD Cost",
    # ]

    # econ_agg = {}
    # for reg in who_reg:
    #     econ_agg[reg] = {}
    #     for scen in scen_name:
    #         econ_agg[reg][scen] = {}
    #         for out in econ_outs:
    #             econ_agg[reg][scen][out] = pd.DataFrame(
    #                 columns=["year", "central"] + [f"run_{i}" for i in range(n_samples)]
    #             )
    #             econ_agg[reg][scen][out].year = np.arange(2000.5, 2051.1, 1)
    #             econ_agg[reg][scen][out] = econ_agg[reg][scen][out].fillna(0)

    # for country in countries:
    #     country_data = sc.load(scens_folder / f"agg_outputs/{country}_econ_eval.pkl")

    #     if country in AFR:
    #         adds = ["AFR", "global"]
    #     if country in AMR:
    #         adds = ["AMR", "global"]
    #     if country in EMR:
    #         adds = ["EMR", "global"]
    #     if country in EUR:
    #         adds = ["EUR", "global"]
    #     if country in SEAR:
    #         adds = ["SEAR", "global"]
    #     if country in WPR:
    #         adds = ["WPR", "global"]
    #     if country in TOP:
    #         adds += ["top10"]

    #     for add in adds:
    #         for scen in scen_name:
    #             for i in range(n_samples + 1):
    #                 # Epidemiological Outs
    #                 epi_agg[add][scen]["HCV incidence"].iloc[:, i + 1] = (
    #                     epi_agg[add][scen]["HCV incidence"].iloc[:, i + 1]
    #                     + country_data["agg_data"][scen]["HCV incidence"].iloc[:, i + 1]
    #                 )
    #                 epi_agg[add][scen]["HCC incidence"].iloc[:, i + 1] = (
    #                     epi_agg[add][scen]["HCC incidence"].iloc[:, i + 1]
    #                     + country_data["agg_data"][scen]["HCC incidence"].iloc[:, i + 1]
    #                 )
    #                 epi_agg[add][scen]["HCV mortality"].iloc[:, i + 1] = (
    #                     epi_agg[add][scen]["HCV mortality"].iloc[:, i + 1]
    #                     + country_data["agg_data"][scen]["HCV mortality"].iloc[:, i + 1]
    #                 )
    #                 epi_agg[add][scen]["DALYs"].iloc[:, i + 1] = (
    #                     epi_agg[add][scen]["DALYs"].iloc[:, i + 1]
    #                     + country_data["econ_ests"][scen]["DALYs"].iloc[:, i + 1]
    #                 )

    #                 # Tests completed
    #                 econ_agg[add][scen]["PWID Ab Tests"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["PWID Ab Tests"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["PWID Ab Tests"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["PWID RNA Tests"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["PWID RNA Tests"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["PWID RNA Tests"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Prisoner Ab Tests"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Prisoner Ab Tests"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Prisoner Ab Tests"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Prisoner RNA Tests"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Prisoner RNA Tests"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Prisoner RNA Tests"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Gen Ab Tests"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Gen Ab Tests"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Gen Ab Tests"].iloc[:, i + 1]
    #                 )
    #                 econ_agg[add][scen]["Gen RNA Tests"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Gen RNA Tests"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Gen RNA Tests"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Total Ab Tests"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["PWID Ab Tests"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Prisoner Ab Tests"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Gen Ab Tests"].iloc[:, i + 1]
    #                 )
    #                 econ_agg[add][scen]["Total RNA Tests"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["PWID RNA Tests"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Prisoner RNA Tests"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Gen RNA Tests"].iloc[:, i + 1]
    #                 )

    #                 # Positive Tests
    #                 econ_agg[add][scen]["PWID Ab Positive"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["PWID Ab Positive"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["PWID Ab Positive"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["PWID RNA Positive"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["PWID RNA Positive"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["PWID RNA Positive"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Prisoner Ab Positive"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Prisoner Ab Positive"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Prisoner Ab Positive"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Prisoner RNA Positive"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Prisoner RNA Positive"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Prisoner RNA Positive"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Gen Ab Positive"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Gen Ab Positive"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Gen Ab Positive"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Gen RNA Positive"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Gen RNA Positive"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Gen RNA Positive"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Total Ab Positive"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["PWID Ab Positive"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Prisoner Ab Positive"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Gen Ab Positive"].iloc[:, i + 1]
    #                 )
    #                 econ_agg[add][scen]["Total RNA Positive"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["PWID RNA Positive"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Prisoner RNA Positive"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Gen RNA Positive"].iloc[:, i + 1]
    #                 )

    #                 # Vaccinations and Treatments Administered
    #                 econ_agg[add][scen]["Total treatment"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Total treatment"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Total treatment"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["PWID vaccinations"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["PWID vaccinations"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["PWID vaccinations"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Prisoner vaccinations"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Prisoner vaccinations"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Prisoner vaccinations"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Gen vaccinations"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Gen vaccinations"].iloc[:, i + 1]
    #                     + country_data["util_data"][scen]["Gen vaccinations"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Total vaccinations"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["PWID vaccinations"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Prisoner vaccinations"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Gen vaccinations"].iloc[:, i + 1]
    #                 )

    #                 # Costs
    #                 econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1]
    #                     + country_data["econ_ests"][scen]["Diagnosis Cost"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1]
    #                     + country_data["econ_ests"][scen][
    #                         "Disease Management Cost"
    #                     ].iloc[:, i + 1]
    #                 )
    #                 econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1]
    #                     + country_data["econ_ests"][scen]["Treatment Cost"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Morbidity Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Morbidity Cost"].iloc[:, i + 1]
    #                     + country_data["econ_ests"][scen][
    #                         "Morbidity Productivity"
    #                     ].iloc[:, i + 1]
    #                 )
    #                 econ_agg[add][scen]["Mortality Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Mortality Cost"].iloc[:, i + 1]
    #                     + country_data["econ_ests"][scen][
    #                         "Mortality Productivity"
    #                     ].iloc[:, i + 1]
    #                 )
    #                 econ_agg[add][scen]["Direct Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1]
    #                 )
    #                 econ_agg[add][scen]["Productivity Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Morbidity Cost"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Mortality Cost"].iloc[:, i + 1]
    #                 )
    #                 econ_agg[add][scen]["Total Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1]
    #                     + econ_agg[add][scen]["Productivity Cost"].iloc[:, i + 1]
    #                 )

    #                 # Vaccination Cost
    #                 econ_agg[add][scen]["Vaccine 5USD Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Vaccine 5USD Cost"].iloc[:, i + 1]
    #                     + country_data["econ_ests"][scen]["Vaccine_5USD"].iloc[:, i + 1]
    #                 )
    #                 econ_agg[add][scen]["Vaccine 10USD Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Vaccine 10USD Cost"].iloc[:, i + 1]
    #                     + country_data["econ_ests"][scen]["Vaccine_10USD"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Vaccine 15USD Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Vaccine 15USD Cost"].iloc[:, i + 1]
    #                     + country_data["econ_ests"][scen]["Vaccine_15USD"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )
    #                 econ_agg[add][scen]["Vaccine 20USD Cost"].iloc[:, i + 1] = (
    #                     econ_agg[add][scen]["Vaccine 20USD Cost"].iloc[:, i + 1]
    #                     + country_data["econ_ests"][scen]["Vaccine_20USD"].iloc[
    #                         :, i + 1
    #                     ]
    #                 )

    # filename = sc.makefilepath(
    #     filename="epi_agg.pkl", folder=scens_folder, makedirs=True
    # )
    # sc.save(filename=filename, obj=epi_agg)

    # filename = sc.makefilepath(
    #     filename="econ_agg.pkl", folder=scens_folder, makedirs=True
    # )
    # sc.save(filename=filename, obj=econ_agg)