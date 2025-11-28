import numpy as np
from hcv import atomica as at
import sciris as sc
import pathlib
import shutil
import os
import pandas as pd
from os.path import isfile, join
from .parameters import iso_to_country
import re
import math
from .atomica.plotting import Series
from scipy.stats import pearsonr

#%% Root directory
def get_project_root():
    """Get the root directory of the project.
    
    This function attempts to locate the root directory of a Python project by searching for
    the presence of specific project files ('pyproject.toml' or 'setup.py') in the parent
    directories of the current file's location. If the current file's location cannot be
    determined, it defaults to the current working directory.
    
    Returns:
        pathlib.Path: The path to the project root directory.
    
    Raises:
        FileNotFoundError: If neither 'pyproject.toml' nor 'setup.py' is found in any
        parent directory.
    """
    
    # Start checking from the current file's directory
    # Note: Use try/except to handle when __file__ is undefined (e.g., REPL)
    try:
        current_dir = pathlib.Path(__file__).resolve().parent
    except NameError:
        # Fallback for REPL/Interactive environment
        current_dir = pathlib.Path.cwd().resolve() 

    # Look for a unique marker file (the "anchor")
    for parent in current_dir.parents:
        if (parent / 'pyproject.toml').exists():
            return parent
        if (parent / 'setup.py').exists():
             return parent
    
    # If no anchor is found after traversing, raise an error
    raise FileNotFoundError("Project root anchor file not found!")

rootdir = get_project_root()

# %% Project & Calibration
def project(
    country: str,
    load_calibration: bool = True,
    cal_folder: str = None,
    cal_version: str = None,
    load_programs: bool = False,
) -> at.Project:
    """Atomica project of HCV vaccine for a specified country.
    
    Args:
        country (str): The name of the country for which the project is being created.
        load_calibration (bool, optional): Flag to indicate whether to load calibration data. Defaults to True.
        cal_folder (str, optional): The folder path where calibration files are stored. Defaults to None.
        cal_version (str, optional): The version of the calibration file to load. Defaults to None.
        load_programs (bool, optional): Flag to indicate whether to load program data. Defaults to False.
    
    Returns:
        at.Project: Atomica Project class containing the country-specific vaccination framework and data.
    """

    # Load project framework
    fw_path = rootdir / "framework" / "hcv_vaccine_framework.xlsx"
    F = at.ProjectFramework(fw_path)
    db_dir = rootdir / "databooks"

    db_files = [f for f in os.listdir(db_dir) if isfile(join(db_dir, f))]
    db = [i for i in db_files if country in i][0]
    db_path = join(db_dir, db)
    D = at.ProjectData.from_spreadsheet(db_path, framework=F)
    
    P = at.Project(framework=F, databook=D, do_run=False)
    P.settings.update_time_vector(start=2000, end=2051, dt=0.25)

    if load_calibration:
        if cal_version is None:
            if (cal_folder / f"{country}_calibration_v2.xlsx").is_file():
                fname = cal_folder / f"{country}_calibration_v2.xlsx"
            else:
                fname = cal_folder / f"{country}_calibration.xlsx"
        else:
            if cal_version == "v1":
                fname = cal_folder / f"{country}_calibration.xlsx"
            else:
                fname = cal_folder / f"{country}_calibration_{cal_version}.xlsx"
        parset = P.make_parset()
        parset = parset.load_calibration(fname)
    else:
        parset = P.make_parset()

    if load_programs:
        _ = P.load_progbook(rootdir / "progbooks" / f"progbook_{country}.xlsx")

    return P


def country_scope(file=None):
    """Retrieve a list of countries from an Excel file.
    
    This function reads an Excel file containing country data and returns a list of countries. 
    If no file is provided, it defaults to reading from a predefined file path.
    
    Args:
        file (str, optional): The path to the Excel file. If None, the function uses 
                              the default file path "/data/country_scope.xlsx".
    
    Returns:
        list: A list of countries extracted from the specified Excel file.
    """
    if file is None:
        df = pd.read_excel("/data/country_scope.xlsx")
    else:
        df = pd.read_excel(file)

    return list(df.Countries)


def return_fw_db(country):
    """Returns the framework and databook for a specified country.
    
    Args:
        country (str): The name of the country for which the databook is to be retrieved.
    
    Returns:
        tuple: A tuple containing:
            - framework (ProjectFramework): The project framework associated with the HCV vaccine.
            - databook (ProjectData): The databook loaded from the corresponding spreadsheet for the specified country.
    
    Raises:
        IndexError: If no databook file is found for the specified country.
        FileNotFoundError: If the framework or databook files cannot be found.
    """
    framework = at.ProjectFramework(
        rootdir / "framework" / "hcv_vaccine_framework.xlsx"
    )
    db_path = rootdir / "databooks"
    db_files = [f for f in os.listdir(db_path) if isfile(join(db_path, f))]
    db = [i for i in db_files if country in i][0]
    databook = at.ProjectData.from_spreadsheet(join(db_path, db), framework=framework)

    return framework, databook


def run_calibration(country, savedir, yaml_folder=None):
    """Run the calibration process for a specified country.
    
    This function initializes a project for the given country, sets up the necessary YAML configuration files based on the country's characteristics, and performs calibration using the specified parameters. The results are saved to an Excel file in the designated directory.
    
    Args:
        country (str): The country code for which the calibration is to be run.
        savedir (Path): The directory where the calibration results will be saved.
        yaml_folder (Path, optional): The folder containing the YAML configuration files. If not provided, defaults to a predefined directory.
    
    Raises:
        ValueError: If the country code is not recognized or supported.
    
    Returns:
        None
    """
    P = project(country, load_calibration=False)
    cal = P.make_parset()

    if yaml_folder is None:
        yaml_folder = rootdir / "yaml"

    # YAML file for pop calibration
    pop_yaml = yaml_folder / "model_config_HepC_pop.yaml"

    # YAML file for epi calibration
    if country in [
        "ARM",
        "BDI",
        "EGY",
        "GAB",
        "GEO",
        "KAZ",
        "KGZ",
        "LVA",
        "MNG",
        "PAK",
        "TJK",
        "UKR",
    ]:
        epi_yaml = (
            yaml_folder / "model_config_HepC_epi_gpop.yaml"
        )  # Countries with gpop transmission
    else:
        epi_yaml = (
            yaml_folder / "model_config_HepC_epi.yaml"
        )  # countries with no gpop transmission

    print(f"\n\nCurrent country: {country} --------------------------------")

    # calibrate pop transfers first
    cal = P.calibrate(
        cal,
        yaml=pop_yaml,
        savedir=savedir / country,
        save_intermediate=False,
        log_output=True,
    )

    # calibrate epi second
    cal = P.calibrate(
        cal,
        yaml=epi_yaml,
        savedir=savedir / country,
        save_intermediate=False,
        log_output=True,
    )
    cal.save_calibration(savedir / f"{country}_calibration.xlsx")


def calibration_outputs(cal_folder, scens_folder):
    """Calibrates and outputs epidemiological data for hepatitis C across various countries.
    
    Args:
        cal_folder (str or Path): The directory where the calibration output Excel file will be saved.
        scens_folder (str or Path): The directory containing scenario output files for economic evaluation.
    
    Returns:
        None: The function saves the calibration outputs to an Excel file and does not return any value.
    
    The function performs the following steps:
    1. Loads incidence data from specified Excel files.
    2. Iterates through WHO regions and their corresponding countries.
    3. For each country, it checks if the country is in the relevant datasets.
    4. Extracts and calculates various epidemiological metrics, including incidence and prevalence estimates.
    5. Compiles the results into dataframes and saves them to an Excel file with multiple sheets.
    """
    
    inci_data_ghr = pd.read_excel(
        str(rootdir) + "/data/GHR24_data.xlsx",
        sheet_name="Table2.11 Countries burden",
    )

    regions_sheet = pd.read_excel(
        str(rootdir) + "/data/flat_datasheet.xlsx",
        sheet_name="Info-Country Region Allocation",
    )
    regions = regions_sheet["WHO region"].unique()
    calib_out = []
    inci_print = []
    inci_plot = []

    inci_polaris = pd.read_excel(
        str(rootdir)
        + "/data/Burnet Polaris Datapull 250528.xlsm",
        sheet_name="2022 estimates",
    )
    countries = country_scope()
    for region in regions:
        print(region)
        countries_df = regions_sheet[regions_sheet["WHO region"] == region]
        for country in countries_df.ISO3:
            if (country not in countries) or (
                country not in inci_data_ghr["ISO3"].values
            ):
                continue
            print(country)
            result = sc.load(scens_folder / "agg_outputs" / f"{country}_econ_eval.pkl")[
                "agg_data"
            ]["counter_0"]
            t = result["HCC incidence"]["year"]
            idx_15 = list(t).index(2015 + 0.5)
            idx_16 = list(t).index(2016 + 0.5)
            idx_22 = list(t).index(2022 + 0.5)
            country_label = countries_df.loc[
                countries_df.ISO3 == country, "Country"
            ].values[0]
            data_ghr = inci_data_ghr[inci_data_ghr.ISO3 == country]
            data_inci_ghr = data_ghr["Hepatitis C incidence"].values[0]
            data_plhcv_ghr = data_ghr[
                "Total hepatitis C infections (all ages) in 2022"
            ].values[0]

            data_polaris = inci_polaris[inci_polaris.ISO3 == country]
            data_inci_polaris = data_polaris["Incidence (chronic)"].values[0]
            data_plhcv_polaris = data_polaris["Prevalent cases"].values[0]
            # multiply polaris inci by spontaneous clearance for total
            spont_clear = sc.load(
                scens_folder
                / "central"
                / f"{country}"
                / "counter_0_central_extracted.pkl"
            )[("spontaneous_clearance", "Total")][idx_22]
            data_inci_polaris = data_inci_polaris * (1 + spont_clear)
            idx_25 = list(t).index(2025 + 0.5)
            idx_26 = list(t).index(2026 + 0.5)
            idx_50 = list(t).index(2050 + 0.5)

            inci_pe = result["HCV incidence"]["central"]
            model_inci = result["HCV incidence"].drop(["year", "central"], axis=1)
            plhcv_pe = result["HCV prevalence"]["central"]
            model_plhcv = result["HCV prevalence"].drop(["year", "central"], axis=1)
            inci_py_pe = result["HCV incidence per 100py"]["central"]
            model_inci_py = result["HCV incidence per 100py"].drop(
                ["year", "central"], axis=1
            )
            # year 2022
            inci_pe_22 = inci_pe[idx_22]
            inci_lb_22 = np.percentile(model_inci.iloc[idx_22], 2.5)
            inci_ub_22 = np.percentile(model_inci.iloc[idx_22], 97.5)
            plhcv_pe_22 = plhcv_pe[idx_22]
            plhcv_lb_22 = np.percentile(model_plhcv.iloc[idx_22], 2.5)
            plhcv_ub_22 = np.percentile(model_plhcv.iloc[idx_22], 97.5)

            inci_print.append(
                [
                    country_label,
                    country,
                    f"{inci_pe_22:,.0f} ({inci_lb_22:,.0f} to {inci_ub_22:,.0f})",
                    f"{data_inci_ghr:,.0f}",
                    f"{data_inci_polaris:,.0f}",
                    f"{plhcv_pe_22:,.0f} ({plhcv_lb_22:,.0f} to {plhcv_ub_22:,.0f})",
                    f"{data_plhcv_ghr:,.0f}",
                    f"{data_plhcv_polaris:,.0f}",
                ]
            )

            inci_plot.append(
                [
                    country_label,
                    inci_pe_22,
                    inci_lb_22,
                    inci_ub_22,
                    data_inci_ghr,
                    data_inci_polaris,
                    plhcv_pe_22,
                    plhcv_lb_22,
                    plhcv_ub_22,
                    data_plhcv_ghr,
                    data_plhcv_polaris,
                ]
            )

            plhcv_2015_pe = plhcv_pe[idx_15]
            plhcv_2015_lb = np.percentile(model_plhcv.iloc[idx_15], 2.5)
            plhcv_2015_ub = np.percentile(model_plhcv.iloc[idx_15], 97.5)
            plhcv_2025_pe = plhcv_pe[idx_25]
            plhcv_2025_lb = np.percentile(model_plhcv.iloc[idx_25], 2.5)
            plhcv_2025_ub = np.percentile(model_plhcv.iloc[idx_25], 97.5)

            inci_py_2015_pe = inci_py_pe[idx_15]
            inci_py_2015_lb = np.percentile(model_inci_py.iloc[idx_15], 2.5)
            inci_py_2015_ub = np.percentile(model_inci_py.iloc[idx_15], 97.5)
            inci_py_2025_pe = inci_py_pe[idx_25]
            inci_py_2025_lb = np.percentile(model_inci_py.iloc[idx_25], 2.5)
            inci_py_2025_ub = np.percentile(model_inci_py.iloc[idx_25], 97.5)

            inci_1625_pe = np.sum(inci_pe.iloc[idx_16 : idx_25 + 1])
            inci_1625_lb = np.percentile(
                np.sum(model_inci.iloc[idx_16 : idx_25 + 1]), 2.5
            )
            inci_1625_ub = np.percentile(
                np.sum(model_inci.iloc[idx_16 : idx_25 + 1]), 97.5
            )
            inci_2650_pe = np.sum(inci_pe.iloc[idx_26 : idx_50 + 1])
            inci_2650_lb = np.percentile(
                np.sum(model_inci.iloc[idx_26 : idx_50 + 1]), 2.5
            )
            inci_2650_ub = np.percentile(
                np.sum(model_inci.iloc[idx_26 : idx_50 + 1]), 97.5
            )

            deaths_pe = result["HCV mortality"]["central"]
            deaths = result["HCV mortality"].drop(["year", "central"], axis=1)
            deaths_1625_pe = np.sum(deaths_pe.iloc[idx_16 : idx_25 + 1])
            deaths_1625_lb = np.percentile(
                np.sum(deaths.iloc[idx_16 : idx_25 + 1]), 2.5
            )
            deaths_1625_ub = np.percentile(
                np.sum(deaths.iloc[idx_16 : idx_25 + 1]), 97.5
            )
            deaths_2650_pe = np.sum(deaths_pe.iloc[idx_26 : idx_50 + 1])
            deaths_2650_lb = np.percentile(
                np.sum(deaths.iloc[idx_26 : idx_50 + 1]), 2.5
            )
            deaths_2650_ub = np.percentile(
                np.sum(deaths.iloc[idx_26 : idx_50 + 1]), 97.5
            )

            calib_out.append(
                [
                    country_label,
                    country,
                    f"{plhcv_2015_pe:,.0f} ({plhcv_2015_lb:,.0f} to {plhcv_2015_ub:,.0f})",
                    f"{plhcv_2025_pe:,.0f} ({plhcv_2025_lb:,.0f} to {plhcv_2025_ub:,.0f})",
                    f"{inci_py_2015_pe*1000:.1f} ({inci_py_2015_lb*1000:.1f} to {inci_py_2015_ub*1000:.1f})",
                    f"{inci_py_2025_pe*1000:.1f} ({inci_py_2025_lb*1000:.1f} to {inci_py_2025_ub*1000:.1f})",
                    f"{inci_1625_pe:,.0f} ({inci_1625_lb:,.0f} to {inci_1625_ub:,.0f})",
                    f"{inci_2650_pe:,.0f} ({inci_2650_lb:,.0f} to {inci_2650_ub:,.0f})",
                    f"{deaths_1625_pe:,.0f} ({deaths_1625_lb:,.0f} to {deaths_1625_ub:,.0f})",
                    f"{deaths_2650_pe:,.0f} ({deaths_2650_lb:,.0f} to {deaths_2650_ub:,.0f})",
                ]
            )

    df_calib = pd.DataFrame(
        calib_out,
        columns=[
            "Country",
            "ISO3",
            "PLHCV 2015",
            "PLHCV 2025",
            "inci 100,000py 2015",
            "inci 100,000py 2025",
            "New infections 2016-2025",
            "New infections 2026-2050",
            "HCV-related deaths 2016-2025",
            "HCV-related deaths 2026-2050",
        ],
    )
    df_inci = pd.DataFrame(
        inci_print,
        columns=[
            "Country",
            "ISO3",
            "Model incidence 2022",
            "GHR24 incidence 2022",
            "Polaris incidence 2022",
            "Model PLHCV 2022",
            "GHR24 PLHCV 2022",
            "Polaris PLHCV 2022",
        ],
    )

    df_plot = pd.DataFrame(
        inci_plot,
        columns=[
            "country",
            "inci_model_pe",
            "inci_model_lb",
            "inci_model_ub",
            "inci_ghr",
            "inci_polaris",
            "plhcv_model_pe",
            "plhcv_model_lb",
            "plhcv_model_ub",
            "plhcv_ghr",
            "plhcv_polaris",
        ],
    )

    savepath = cal_folder / "calibration_validation.xlsx"
    with pd.ExcelWriter(savepath, engine="xlsxwriter") as writer:
        df_calib.to_excel(writer, sheet_name="Countries outcomes", index=False)
        df_inci.to_excel(writer, sheet_name="Countries inci validation", index=False)
        df_plot.to_excel(writer, sheet_name="Plot data", index=False)
    print(f"Saved calibration outputs: {savepath}")


# %% Vaccine scenarios
def define_scenarios(P, progset):
    """Defines scenarios based on program settings and coverage data.
    
    Args:
        P (object): Atomica project.
        progset (object): An object containing a set of programs to be analyzed.
    
    Returns:
        tuple: A tuple containing:
            - list: A list of program instructions for each scenario.
            - list: A list of unique scenario names.
    
    Raises:
        FileNotFoundError: If the specified Excel file does not exist.
        ValueError: If there are issues with the data in the Excel file.
    """
    # Control programs with coverage since we don't care about costs (calculated outside)
    t_vec = list(np.arange(P.settings.sim_start, P.settings.sim_end))

    data_progbook = str(rootdir) + "/data/progbook_inputs.xlsx"
    df_scenarios = pd.read_excel(pd.ExcelFile(data_progbook), sheet_name="scenarios")
    scenarios = list(df_scenarios.scenario_name.unique())

    progset_instructions = []
    for scen_name in scenarios:
        df = df_scenarios[df_scenarios.scenario_name == scen_name]
        coverage = {
            prog_all: at.TimeSeries(t_vec, [0] * len(t_vec))
            for prog_all in progset.programs
        }
        for _, row in df.iterrows():
            prog = row["programs"]
            start_year = row["start_year"]
            end_year = row["end_year"]
            if P.progsets[0].programs[prog].is_one_off:
                cov_vals = [
                    0 if (t < start_year or t > end_year) else 1 / P.settings.sim_dt
                    for t in t_vec
                ]
            else:
                cov_vals = [0 if (t < start_year or t > end_year) else 1 for t in t_vec]
            coverage[prog] = at.TimeSeries(t_vec, cov_vals)
        instructions = at.ProgramInstructions(start_year=2026, coverage=coverage)
        progset_instructions.append(instructions)

    return progset_instructions, scenarios


def run_scenario_sampling(country, cal_folder, rand_seed, n_samples, savedir):
    """Run scenario sampling for a given country and save the results.
    
    This function initializes a project for the specified country, loads calibration data, and runs simulations based on defined scenarios. It generates outputs for various population groups and aggregates results, saving them to specified directories.
    
    Args:
        country (str): The name of the country for which the scenario sampling is to be run.
        cal_folder (str): The folder path where calibration files are located.
        rand_seed (int): The random seed for reproducibility of the sampling process.
        n_samples (int): The number of samples to generate during the simulation.
        savedir (str): The directory where the output files will be saved.
    
    Returns:
        None: The function saves output files to the specified directory and does not return any value.
    """
    P = project(
        country, load_calibration=True, cal_folder=cal_folder, load_programs=False
    )
    pops = P.data.pops.keys()  # define population

    # Function to define outputs to save
    def mapping_function(x):
        # Outputs to be saved
        outputs_by_pop = [
            "alive",
            "prevalence",
            "total_hcv",
            "inci_all_m",
            "notifications_m",
            "inci_m",
            "ab_tests_m",
            "pcr_tests_m",
        ]
        outputs_tot = [
            "total_hcv",
            "deaths_hcv_total",
            "tx_m",
            "inci_m",
            "notifications_m",
            "ab_tests_m",
            "pcr_tests_m",
            "f4_all",
            "d_cirrhosis_all",
            "inci_all_m",
            "hcc_inci",
            "acute_all",
            "undiag_all",
            "abpos_chronic",
            "pcr",
            "treated",
            "dc_udx",
            "hcc_udx",
            "f0f2_dx",
            "f3_dx",
            "f4_dx",
            "dc_dx",
            "hcc_dx",
            "d_cirrhosis_all",
            "cancer_all",
        ]
        outputs_tot_weighted = [
            "prop_f0f2",
            "prop_f3",
            "prop_f4",
            "prop_dc",
            "prop_hcc",
            "prevalence",
            "inci_100py",
        ]

        outputs_all = at.PlotData(
            x, outputs=outputs_by_pop, pops=pops, t_bins=1
        )  # outputs for each population group
        extra = at.PlotData(
            x, outputs=outputs_tot, pops="total", pop_aggregation="sum", t_bins=1
        )  # total outputs / summed
        outputs_all.series.extend(extra.series)
        extra = at.PlotData(
            x,
            outputs=outputs_tot_weighted,
            pops="total",
            pop_aggregation="weighted",
            t_bins=1,
        )  # total outputs / weighed
        outputs_all.series.extend(extra.series)

        # Aging pars
        transfer_pars = [
            "age_10-17_males_to_18-64_males",
            "age_10-17_females_to_18-64_females",
        ]
        transfer_pops = ["10-17_males", "10-17_females"]
        for par, pop in zip(transfer_pars, transfer_pops):
            extra = at.PlotData(x, outputs=par, pops=pop, t_bins=1)
            outputs_all.series.extend(extra.series)

        # Specific pop aggregations for econ fcts
        pars_tot = [
            "inci_all_m",
            "ab_all_m",
            "pcr_all_m",
            "f0f3_utx",
            "f4_utx",
            "f0f3_c",
            "f4_c",
            "inci_100py",
        ]
        pops_tot = [
            {
                "gen_pop": [
                    "0-9_males",
                    "0-9_females",
                    "10-17_males",
                    "10-17_females",
                    "18-64_males",
                    "18-64_females",
                    "65+_males",
                    "65+_females",
                ]
            },
            {"pwid": ["PWID_males", "PWID_females"]},
            {"prisoners": ["Prisoners_males", "Prisoners_females"]},
            {
                "working_age": [
                    "18-64_males",
                    "18-64_females",
                    "PWID_males",
                    "PWID_females",
                ]
            },
        ]
        for par in pars_tot:
            for pop in pops_tot:
                extra = at.PlotData(
                    x, outputs=par, pops=pop, pop_aggregation="sum", t_bins=1
                )
                outputs_all.series.extend(extra.series)

        pars_tot = ["rna_test_efficiency", "ab_test_efficiency", "prevalence"]
        extra = at.PlotData(
            x, outputs=pars_tot, pops=pops_tot, pop_aggregation="weighted", t_bins=1
        )
        outputs_all.series.extend(extra.series)

        pops_tot = [
            {"children": ["0-9_males", "0-9_females"]},
            {"teens": ["10-17_males", "10-17_females"]},
            {"adults": ["18-64_males", "18-64_females", "65+_males", "65+_females"]},
            {"pwid": ["PWID_males", "PWID_females"]},
            {"prisoners": ["Prisoners_males", "Prisoners_females"]},
        ]
        par = "vaccine_uptake:flow"
        extra = at.PlotData(
            x, outputs=par, pops=pops_tot, pop_aggregation="sum", t_bins=1
        )
        outputs_all.series.extend(extra.series)

        # Weighted spontaneous clearance rate
        tvec = (
            at.PlotData(x, outputs="acute", pops=pops, pop_aggregation="sum", t_bins=1)
            .series[0]
            .tvec
        )
        numerator = np.zeros_like(tvec)
        denominator = np.zeros_like(tvec)
        for pop in pops:
            acute = at.PlotData(x, outputs="acute", pops=pop, t_bins=1).series[0].vals
            acute_vax = (
                at.PlotData(x, outputs="acute_vax_reinfected", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            clear = (
                at.PlotData(x, outputs="prop_clear", pops=pop, t_bins=1).series[0].vals
            )
            clear_vax = (
                at.PlotData(x, outputs="vaccine_prop_clear", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            acute_vax_i = (
                at.PlotData(x, outputs="acute_vax_initial", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            clear_vax_i = (
                at.PlotData(x, outputs="vaccine_prop_clear_initial", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            numerator += (
                acute * clear + acute_vax * clear_vax + acute_vax_i * clear_vax_i
            )
            denominator += acute + acute_vax + acute_vax_i
        with np.errstate(divide="ignore", invalid="ignore"):
            weighted_clearance = np.where(denominator > 0, numerator / denominator, 0)
        extra = Series(
            tvec=tvec,
            vals=weighted_clearance,
            pop="Total",
            output="spontaneous_clearance",
            units="probability",
        )
        outputs_all.series.append(extra)

        return outputs_all

    # Function to save outputs
    def write_outputs(result, i):
        if isinstance(result, list):
            assert len(result) == 1, "result is a list with more than 1 item"
            result = result[0]
        outputs_all = mapping_function(result)
        name = result.name
        extracted = dict()
        for serie in outputs_all.series:
            tvec = serie.tvec
            pop = serie.pop
            par = serie.output
            vals = serie.vals
            extracted[(par, pop)] = vals
        extracted[("_name")] = name
        extracted[("_tvec")] = tvec

        filename = sc.makefilepath(
            filename=f"{name}_{i}_extracted.pkl", folder=savedir, makedirs=True
        )
        sc.save(filename=filename, obj=extracted)
        print(f"Saved outputs: {filename}")

    # Function to generate progbook
    savedir_pb = pathlib.Path(savedir).parents[1]
    savedir_pb = savedir_pb / "progbooks"
    savedir_pb.mkdir(parents=True, exist_ok=True)

    def gen_pb(result, i):
        from hcv.generate_progbooks import generate_progbook
        if isinstance(result, list):
            assert len(result) == 1, "result is a list with more than 1 item"
            result = result[0]
        generate_progbook(
            country, result=result, savedir=savedir_pb / f"progbook_{country}_{i}.xlsx"
        )

    parset = P.make_parset()
    parset.load_calibration(cal_folder / f"{country}_calibration.xlsx")
    # Run with no sampling first
    result_central = P.run_sim(parset=parset)
    gen_pb(result_central, "central")
    pset = P.load_progbook(savedir_pb / f"progbook_{country}_central.xlsx")
    progset_instructions, scenarios = define_scenarios(P, pset)
    # Save central results to separate central folder as well
    savedir2 = savedir.parents[0] / "central" / f"{country}"
    savedir2.mkdir(parents=True, exist_ok=True)
    for p_i, scen in zip(progset_instructions, scenarios):
        result_central = P.run_sim(
            parset=parset, progset=pset, progset_instructions=p_i, result_name=scen
        )
        write_outputs(result_central, "central")
        filename = f"{scen}_central_extracted.pkl"
        shutil.copy(savedir / filename, savedir2 / filename)

    # Run with sampling
    progset = []
    P.run_sampled_sims(
        parset=parset, n_samples=n_samples, rand_seed=rand_seed, custom_fct=gen_pb
    )
    for i in np.arange(n_samples):
        pset = P.load_progbook(savedir_pb / f"progbook_{country}_{i}.xlsx")
        progset.append(pset)
    for p_i, scen in zip(progset_instructions, scenarios):
        P.run_sampled_sims(
            parset=parset,
            progset=progset,
            progset_instructions=p_i,
            result_names=[scen],
            n_samples=n_samples,
            rand_seed=rand_seed,
            custom_fct=write_outputs,
        )


def calc_averted_region(results_folder, file_name="epi_agg.pkl", n_samples=100):
    """Calculates averted health outcomes for different scenarios and regions based on epidemiological data.
    
    Args:
        results_folder (Path): The folder containing the results data.
        file_name (str, optional): The name of the file containing the results data. Defaults to "epi_agg.pkl".
        n_samples (int, optional): The number of samples to use for calculating averted outcomes. Defaults to 100.
    
    Returns:
        pd.DataFrame: A DataFrame containing the averted outcomes, including point estimates, lower and upper bounds, 
        for each region and scenario.
    
    Raises:
        FileNotFoundError: If the specified results file does not exist.
        ValueError: If the data in the results file is not in the expected format.
    """
    out = sc.load(results_folder / file_name)
    data_progbook = str(rootdir) + "/data/progbook_inputs.xlsx"
    df_scenarios = pd.read_excel(pd.ExcelFile(data_progbook), sheet_name="scenarios")
    df = df_scenarios[df_scenarios.counterfactual.str.len() > 0]
    scenario_counter = dict(zip(df["scenario_name"], df["counterfactual"]))
    alpha = 0.05

    epi_outcomes = ["HCV incidence", "HCC incidence", "HCV mortality", "DALYs"]
    epi_averted = [
        "HCV infections averted",
        "HCC cases averted",
        "HCV deaths averted",
        "DALYs averted",
    ]
    years = np.arange(2026, 2051)

    averted_dict = {
        region: {
            scen: pd.DataFrame(
                index=epi_averted,
                columns=["point estimate", "lower bound", "upper bound"],
            )
            for scen in scenario_counter
        }
        for region in out
    }

    for region in out:
        for scen, counter in scenario_counter.items():
            for par, par_avert in zip(epi_outcomes, epi_averted):
                scen_df = out[region][scen][par].set_index("year")
                counter_df = out[region][counter][par].set_index("year")

                years_shifted = years + 0.5
                runs = [f"run_{i}" for i in range(n_samples)]

                central_averted = (
                    counter_df.loc[years_shifted, "central"]
                    - scen_df.loc[years_shifted, "central"]
                )
                point_estimate = central_averted.sum()

                run_diff = (
                    counter_df.loc[years_shifted, runs]
                    - scen_df.loc[years_shifted, runs]
                )
                total_averted = run_diff.sum(axis=0)
                lower_bound = np.percentile(total_averted, alpha / 2 * 100)
                upper_bound = np.percentile(total_averted, (1 - alpha / 2) * 100)

                averted_dict[region][scen].loc[par_avert] = [
                    point_estimate,
                    lower_bound,
                    upper_bound,
                ]

    data = []
    scenarios = df_scenarios.scenario_name.unique()
    scenario_name = dict()
    for scen in scenarios:
        scenario_name[scen] = df_scenarios[
            df_scenarios.scenario_name == scen
        ].description.values[0]
    for region, scen_data in averted_dict.items():
        for scenario, df in scen_data.items():
            for outcome in df.index:
                row = df.loc[outcome]
                data.append(
                    {
                        "Region": region,
                        "Scenario": scenario_name[scenario].replace(" + ", "\n"),
                        "Outcome": outcome,
                        "Point Estimate": row["point estimate"],
                        "Lower Bound": row["lower bound"],
                        "Upper Bound": row["upper bound"],
                        "Error": (row["upper bound"] - row["lower bound"]) / 2,
                    }
                )
    return pd.DataFrame(data)


def calc_outcomes_region(scens_folder, n_samples=100, regions=None):
    """Calculates epidemiological and economic outcomes for specified regions based on scenario data.
    
    Args:
        scens_folder (Path): The folder containing scenario data files.
        n_samples (int, optional): The number of samples to use for calculations. Defaults to 100.
        regions (list, optional): A list of regions to calculate outcomes for. If None, outcomes for all regions will be calculated. Defaults to None.
    
    Returns:
        dict: A dictionary containing the calculated outcomes for each region, including central estimates and confidence intervals for various epidemiological parameters.
    """
    out = sc.load(scens_folder / "epi_agg.pkl")
    out_econ = sc.load(scens_folder / "econ_agg.pkl")
    data_progbook = str(rootdir) + "/data/progbook_inputs.xlsx"
    df_scenarios = pd.read_excel(pd.ExcelFile(data_progbook), sheet_name="scenarios")
    df = df_scenarios[df_scenarios.counterfactual.str.len() > 0]
    scenario_counter = dict(zip(df["scenario_name"], df["counterfactual"]))
    scenario_description = dict(zip(df["scenario_name"], df["description"]))
    alpha = 0.05

    if regions is None:
        regions = list(out.keys())

    out = {region: out[region] for region in regions}
    out_econ = {
        region: {
            scen: {"Total vaccinations": out_econ[region][scen]["Total vaccinations"]}
            for scen in scenario_counter
        }
        for region in regions
    }

    for region in regions:
        for scen in scenario_counter:
            out[region][scen].update(out_econ[region][scen])

    epi_outcomes = [
        "HCV incidence",
        "HCC incidence",
        "HCV mortality",
        "Total vaccinations",
    ]
    epi_outcomes_print = [
        "HCV incidence",
        "HCC incidence",
        "HCV deaths",
        "Vaccines administered",
    ]
    plot_data = {
        region: {par2: {} for par2 in epi_outcomes_print} for region in regions
    }
    years = np.arange(2020, 2051)

    for region in out:
        for scen, counter in scenario_counter.items():
            for par, par2 in zip(epi_outcomes, epi_outcomes_print):
                scen_df = out[region][scen][par].set_index("year")
                if par == "Total vaccinations":
                    counter_df = scen_df.copy()
                    counter_df.loc[:, :] = 0
                else:
                    counter_df = out[region][counter][par].set_index("year")

                years_shifted = years + 0.5
                runs = [f"run_{i}" for i in range(n_samples)]

                counter_central = list(counter_df.loc[years_shifted, "central"])
                counter_lb = np.percentile(
                    counter_df.loc[years_shifted, runs], alpha / 2 * 100, axis=1
                )
                counter_ub = np.percentile(
                    counter_df.loc[years_shifted, runs], (1 - alpha / 2) * 100, axis=1
                )

                scen_central = list(scen_df.loc[years_shifted, "central"])
                scen_lb = list(
                    np.percentile(
                        scen_df.loc[years_shifted, runs], alpha / 2 * 100, axis=1
                    )
                )
                scen_ub = list(
                    np.percentile(
                        scen_df.loc[years_shifted, runs], (1 - alpha / 2) * 100, axis=1
                    )
                )

                # if 'Vaccines' in par2:
                plot_data[region][par2][scenario_description[scen]] = {
                    "years": years,
                    "central": scen_central,
                    "lb": scen_lb,
                    "ub": scen_ub,
                    "cf_central": counter_central,
                }

    return plot_data


# %% Misc calculations
def calculate_pop_transfers(res):
    """Calculate population transfers based on model parameters.
    
    This function processes the model results to compute the transfer flows between different population groups based on specified parameters. It identifies relevant parameters, extracts transfer information, and aggregates the results into a structured dictionary.
    
    Args:
        res: An object containing the model results, which includes population data and parameters.
    
    Returns:
        dict: A dictionary where each key is a population group and the value is another dictionary containing:
            - 'tvec': A time vector associated with the transfers.
            - 'transfers': A list of computed transfer values for the population group.
    """

    # population transfers (summer over for each pop group)
    all_pops = [
        "PWID_males",
        "PWID_females",
        "0-9_males",
        "0-9_females",
        "10-17_males",
        "10-17_females",
        "18-64_males",
        "18-64_females",
        "65+_males",
        "65+_females",
        "Prisoners_males",
        "Prisoners_females",
    ]
    transfer_codes = ["idu_status_", "inc_", "age_"]
    transfer_flows_by_pop = {pop: None for pop in all_pops}

    disease_stages = ["f0", "f1", "f2", "f3", "f4", "dc", "hcc"]
    cascade_stages = [
        "ud_ltc",
        "ab_ud",
        "ab_d",
        "ab_ltfu",
        "pcr_d",
        "pcr_ltfu",
        "t1",
        "ft1",
        "ft1_ltfu",
        "t2",
        "ft2",
    ]
    list_stages = ["{}_{}".format(d, c) for d in disease_stages for c in cascade_stages]

    for pop in res.model.pops:
        for par in pop.pars:
            par_name = par.name
            if not any(code in par_name for code in transfer_codes):
                continue
            if "_to_" not in par_name:
                continue
            src, dst = par_name.rsplit("_to_", 1)
            for pop in all_pops:
                if pop in src:
                    src = pop
                if pop in dst:
                    dst = pop
            # src = src.split('_')[-1] # get pop name after last underscore (example idu_status_)

            if src not in all_pops and dst not in all_pops:
                continue

            # Get flow from the source pop only
            tvec = (
                at.PlotData(
                    res, outputs=f"{list_stages[0]}::{par_name}", pops=src, t_bins=1
                )
                .series[0]
                .tvec
            )
            vals = np.zeros(len(tvec))
            for comp in list_stages:
                vals += np.array(
                    at.PlotData(res, outputs=f"{comp}::{par_name}", pops=src, t_bins=1)
                    .series[0]
                    .vals
                )

            for pop_name, sign in [(src, -1), (dst, +1)]:
                if pop_name not in all_pops:
                    continue
                if transfer_flows_by_pop[pop_name] is None:
                    transfer_flows_by_pop[pop_name] = [0.0] * len(vals)
                transfer_flows_by_pop[pop_name] = [
                    a + sign * b for a, b in zip(transfer_flows_by_pop[pop_name], vals)
                ]
    transfers_dict = dict(map(lambda key: (key, dict()), all_pops))
    for pop in all_pops:
        transfers_dict[pop] = {"tvec": tvec, "transfers": transfer_flows_by_pop[pop]}

    return transfers_dict


def format_number_with_sf_and_suffix(number, sf=3):
    """Formats a number with a specified significant figure and a suffix based on its magnitude.
    
    Args:
        number (float): The number to be formatted.
        sf (int, optional): The number of significant figures to display. Defaults to 3.
    
    Returns:
        str: The formatted number as a string, including a suffix if applicable.
    
    Examples:
        >>> format_number_with_sf_and_suffix(1234567)
        '1.23M'
        >>> format_number_with_sf_and_suffix(0)
        '0'
        >>> format_number_with_sf_and_suffix(9876543210, sf=2)
        '9.88B
    """
    if number == 0:
        return "0"
    abs_num = abs(number)

    e = math.floor(math.log10(abs_num) / 3)
    suffixes = ["", "k", "M", "B", "T", "P", "E"]

    if e < 1:
        return f"{number:.{sf}g}"

    scale = 10 ** (e * 3)
    scaled_num = number / scale
    suffix = suffixes[e]
    formatted_num_str = f"{scaled_num:.{sf}g}"

    return formatted_num_str + suffix


def format_structured_string(input_str):
    """Formats a structured string by converting numeric patterns into a standardized format.
    
    This function takes an input string, cleans it by replacing non-breaking spaces, and searches for numeric patterns. It converts these patterns into a float and formats them with significant figures and suffixes. If the conversion fails, it retains the original match.
    
    Args:
        input_str (str): The input string containing numeric patterns to be formatted.
    
    Returns:
        str: The formatted string with numeric patterns converted and cleaned.
    """

    # 1. Precise Regex to capture the full number, including decimal part
    # Group 1: Integer part with optional commas (e.g., "10,812")
    # Group 2: Optional decimal part (e.g., ".51")
    pattern = r"(\d{1,3}(?:,\d{3})*)(?:\.(\d+))?"

    # Clean the input of any non-breaking spaces before processing
    clean_input = input_str.replace("\xa0", " ")

    # Helper to clean and convert the captured number string
    def clean_and_convert(integer_part, decimal_part):
        # 1. Construct the full number string: "10,812" + "." + "51" (if exists)
        num_str_dirty = integer_part
        if decimal_part:
            num_str_dirty += "." + decimal_part

        # 2. Remove all commas (thousands separator)
        cleaned = num_str_dirty.replace(",", "")

        # 3. Convert to float
        try:
            return float(cleaned)
        except ValueError:
            print(f"Warning: Could not convert '{num_str_dirty}' after cleaning.")
            return None

    def replacer(match):
        # match.group(1) is the integer part (e.g., "10,812")
        # match.group(2) is the decimal part (e.g., "51")

        integer_part = match.group(1)
        decimal_part = match.group(
            2
        )  # Will be None if no decimal part was in the input

        number = clean_and_convert(integer_part, decimal_part)

        return (
            format_number_with_sf_and_suffix(number)
            if number is not None
            else match.group(0)
        )  # Return original on failure

    # --- Step 2 & 3: Replace all captured numbers in the string ---
    temp_output = re.sub(pattern, replacer, clean_input)

    # --- Step 4: Re-insert the space before the parenthesis if missing ---
    # Pattern: [digit or k/M/B] + '('
    # Replacement: [digit or k/M/B] + ' ('
    final_output = re.sub(r"([0-9kMBTPE])\(", r"\1 (", temp_output)

    return final_output


def bcr_correlation():
    """Calculates the correlation between various predictors and BCR scenarios for different countries.
    
    This function reads data from CSV and Excel files, processes the data to extract relevant information for each country, and computes the Pearson correlation coefficients between GDP, prevalence of HCV RNA among PWID, percentage ever tested, and percentage ever treated against BCR scenarios. The results are saved to an Excel file.
    
    Args:
        None
    
    Returns:
        None
    
    Raises:
        FileNotFoundError: If the specified CSV or Excel files do not exist.
        KeyError: If expected columns are not found in the dataframes.
    
    Example:
        To use this function, simply call bcr_correlation() after setting the appropriate rootdir variable.
    """

    # Load BCR values (scenarios x countries)
    bcr_country = pd.read_csv(
        str(rootdir) + "results/scenarios/bcr_map.csv", index_col=0
    )
    scenarios = list(bcr_country.columns)
    countries = list(bcr_country.index)

    # Load predictors
    gdp_all = pd.read_excel(
        str(rootdir) + "/data/flat_datasheet.xlsx",
        sheet_name="Cost - YLL and productivity",
    )
    prev_all = pd.read_excel(
        str(rootdir) + "/data/flat_datasheet.xlsx",
        sheet_name="Country-Prison Prevalence",
    )
    cc_all = pd.read_excel(
        str(rootdir) + "/data/flat_datasheet.xlsx",
        sheet_name="Country-PWID Care Cascade",
    )

    data_list = []

    # Build merged dataset for correlation
    for country in countries:
        gdp = gdp_all.loc[gdp_all.ISO3 == country, "GDP"].values[0]
        prev = (
            prev_all.loc[prev_all.ISO3 == country, "HCV RNA prevalence (PWID)"]
            .dropna()
            .values[-1]
        )
        test = (
            cc_all.loc[cc_all.ISO3 == country, "Ever tested (proportion)"]
            .dropna()
            .values[-1]
        )
        treat = (
            cc_all.loc[cc_all.ISO3 == country, "Ever treated (proportion)"]
            .dropna()
            .values[-1]
        )

        row = {
            "country": country,
            "GDP": gdp,
            "Prevalence (PWID RNA)": prev,
            "% Ever Tested": test,
            "% Ever Treated": treat,
        }

        # Add BCR values for each scenario
        for scen in scenarios:
            row[f"BCR {scen}"] = bcr_country.loc[country, scen]

        data_list.append(row)

    df = pd.DataFrame(data_list)  # remove countries missing any key data

    # Generate correlation results
    predictors = ["GDP", "Prevalence (PWID RNA)", "% Ever Tested", "% Ever Treated"]
    corr_rows = []
    print_rows = []

    for scen in scenarios:
        scen_col = f"BCR {scen}"
        for var in predictors:
            r, p = pearsonr(df[var], df[scen_col])
            corr_rows.append(
                {
                    "Scenario": scen,
                    "Variable": var,
                    "Pearson r": round(r, 3),
                    "p-value": round(p, 4),
                }
            )
            # scen_row.append(f'r={round(r,3)} (p={round(p,3)})')

    corr_df = pd.DataFrame(corr_rows)
    for var in predictors:
        corr_var = corr_df[corr_df.Variable == var]
        scen_row = []
        for scen in scenarios:
            corr = corr_var[corr_var.Scenario == scen]
            r = corr["Pearson r"].values[0]
            p = corr["p-value"].values[0]
            scen_row.append(f"{r} ({p})")
        print_rows.append(scen_row)
    print_df = pd.DataFrame(print_rows, index=predictors, columns=scenarios)

    # Write to Excel
    out_path = str(rootdir) + "results/scenarios/bcr_correlations.xlsx"
    print_df.to_excel(out_path, index=True)

    print(f"Correlation results saved to:{out_path}")


# %% Economic functions
def econ_eval(country, savedir_scens, results_folder, rand_seed, n_samples):
    """Evaluate economic scenarios for a given country based on epidemiological data.
    
    Args:
        country (str): The ISO3 code of the country for which the evaluation is performed.
        savedir_scens (Path): The directory where scenario data is saved.
        results_folder (Path): The folder where results will be stored.
        rand_seed (int): The random seed for reproducibility of results.
        n_samples (int): The number of samples to generate for the evaluation.
    
    Returns:
        None: The function saves the evaluation results to a specified file.
        
    Details:
        This function loads scenario data from an Excel file, processes various epidemiological outputs, 
        and computes economic estimates related to healthcare costs, productivity losses, and vaccination 
        impacts over a specified time frame. The results are saved in a structured format for further analysis.
    Author: @ChrisSeaman-Burnet
    """

    # Import scenario names and definitions for loop and mapping
    scenarios = pd.read_excel(
        str(rootdir) + "/data/progbook_inputs.xlsx", sheet_name="scenarios"
    )
    scen_name = list(pd.unique(scenarios.scenario_name))
    sim_loop = ["central"] + list(np.arange(0, n_samples, 1))

    # Aggregate Summary Data (deaths will be used for counterfactual YLL/YPLLs)
    agg_data = {}
    for scen in scen_name:
        agg_data[scen] = {}

    colnames = ["year", "central"] + [f"run_{i}" for i in range(n_samples)]
    agg_outs = [
        "HCV incidence",
        "HCC incidence",
        "HCV incidence per 100py",
        "HCV mortality",
        "HCV prevalence",
        "PWID prevalence",
        "Prisoner prevalence",
    ]
    for scen in scen_name:
        for out in agg_outs:
            agg_data[scen][out] = pd.DataFrame(columns=colnames)
            agg_data[scen][out]["year"] = np.arange(2000.5, 2051.5, 1)

    for scen in scen_name:
        for i, sim in enumerate(sim_loop):
            data = sc.load(savedir_scens / f"{scen}_{sim}_extracted.pkl")
            # Total Population
            agg_data[scen]["HCV incidence"].iloc[:, i + 1] = data["inci_all_m", "Total"]
            agg_data[scen]["HCC incidence"].iloc[:, i + 1] = data["hcc_inci", "Total"]
            agg_data[scen]["HCV incidence per 100py"].iloc[:, i + 1] = data[
                "inci_100py", "Total"
            ]
            agg_data[scen]["HCV mortality"].iloc[:, i + 1] = data[
                "deaths_hcv_total", "Total"
            ]
            agg_data[scen]["HCV prevalence"].iloc[:, i + 1] = data["total_hcv", "Total"]
            # Key Populations
            agg_data[scen]["PWID prevalence"].iloc[:, i + 1] = (
                data["total_hcv", "Prisoners_males"]
                + data["total_hcv", "Prisoners_females"]
            )
            agg_data[scen]["Prisoner prevalence"].iloc[:, i + 1] = (
                data["total_hcv", "PWID_males"] + data["total_hcv", "PWID_females"]
            )

    # Diagnostic Testing, Positive Tests, Treatment and Vaccination Coverage (most can be used for costs)
    util_data = {}
    for scen in scen_name:
        util_data[scen] = {}

    util_outs = [
        "PWID Ab Tests",
        "Prisoner Ab Tests",
        "Gen Ab Tests",
        "Total Ab Tests",
        "PWID Ab Positive",
        "Prisoner Ab Positive",
        "Gen Ab Positive",
        "Total Ab Positive",
        "PWID RNA Tests",
        "Prisoner RNA Tests",
        "Gen RNA Tests",
        "Total RNA Tests",
        "PWID RNA Positive",
        "Prisoner RNA Positive",
        "Gen RNA Positive",
        "Total RNA Positive",
        "Total treatment",
        "PWID vaccinations",
        "Prisoner vaccinations",
        "Gen vaccinations",
        "Total vaccinations",
    ]

    # Extract 'spontaneous_clearance' parameter from central runs
    spont_clear = {}
    for scen in scen_name:
        data = sc.load(
            str(rootdir)
            + f"results/scenarios/central/{country}/{scen}_central_extracted.pkl"
        )
        spont_clear[scen] = pd.DataFrame(columns=["year", "spontaneous_clearance"])
        spont_clear[scen].year = np.arange(2000.5, 2051.5, 1)
        spont_clear[scen].spontaneous_clearance = data["spontaneous_clearance", "Total"]

    for scen in scen_name:
        for out in util_outs:
            util_data[scen][out] = pd.DataFrame(columns=colnames)
            util_data[scen][out]["year"] = np.arange(2000.5, 2051.5, 1)

    for scen in scen_name:
        for i, sim in enumerate(sim_loop):
            data = sc.load(savedir_scens / f"{scen}_{sim}_extracted.pkl")

            # Antibody Tests (total)
            util_data[scen]["PWID Ab Tests"].iloc[:, i + 1] = (
                data["ab_tests_m", "PWID_males"] + data["ab_tests_m", "PWID_females"]
            )
            util_data[scen]["Prisoner Ab Tests"].iloc[:, i + 1] = (
                data["ab_tests_m", "Prisoners_males"]
                + data["ab_tests_m", "Prisoners_females"]
            )
            util_data[scen]["Total Ab Tests"].iloc[:, i + 1] = data[
                "ab_tests_m", "Total"
            ]
            util_data[scen]["Gen Ab Tests"].iloc[:, i + 1] = util_data[scen][
                "Total Ab Tests"
            ].iloc[:, i + 1] - (
                util_data[scen]["PWID Ab Tests"].iloc[:, i + 1]
                + util_data[scen]["Prisoner Ab Tests"].iloc[:, i + 1]
            )  # acts as a double check

            # Antibody Tests (positive)
            util_data[scen]["PWID Ab Positive"].iloc[:, i + 1] = data[
                "ab_all_m", "pwid"
            ]
            util_data[scen]["Prisoner Ab Positive"].iloc[:, i + 1] = data[
                "ab_all_m", "prisoners"
            ]
            util_data[scen]["Gen Ab Positive"].iloc[:, i + 1] = data[
                "ab_all_m", "gen_pop"
            ]
            util_data[scen]["Total Ab Positive"].iloc[:, i + 1] = (
                util_data[scen]["PWID Ab Positive"].iloc[:, i + 1]
                + util_data[scen]["Prisoner Ab Positive"].iloc[:, i + 1]
                + util_data[scen]["Gen Ab Positive"].iloc[:, i + 1]
            )  # acts as a double check

            # RNA tests (total)
            util_data[scen]["PWID RNA Tests"].iloc[:, i + 1] = data[
                "ab_all_m", "pwid"
            ] * (1 + spont_clear[scen].iloc[:, 1])
            util_data[scen]["Prisoner RNA Tests"].iloc[:, i + 1] = data[
                "ab_all_m", "prisoners"
            ] * (1 + spont_clear[scen].iloc[:, 1])
            util_data[scen]["Gen RNA Tests"].iloc[:, i + 1] = data[
                "ab_all_m", "gen_pop"
            ] * (1 + spont_clear[scen].iloc[:, 1])
            util_data[scen]["Total RNA Tests"].iloc[:, i + 1] = (
                util_data[scen]["PWID RNA Tests"].iloc[:, i + 1]
                + util_data[scen]["Prisoner RNA Tests"].iloc[:, i + 1]
                + util_data[scen]["Gen RNA Tests"].iloc[:, i + 1]
            )  # acts as a double check

            # RNA tests (positive)
            util_data[scen]["PWID RNA Positive"].iloc[:, i + 1] = data[
                "pcr_all_m", "pwid"
            ]
            util_data[scen]["Prisoner RNA Positive"].iloc[:, i + 1] = data[
                "pcr_all_m", "prisoners"
            ]
            util_data[scen]["Gen RNA Positive"].iloc[:, i + 1] = data[
                "pcr_all_m", "gen_pop"
            ]
            util_data[scen]["Total RNA Positive"].iloc[:, i + 1] = (
                util_data[scen]["PWID RNA Positive"].iloc[:, i + 1]
                + util_data[scen]["Prisoner RNA Positive"].iloc[:, i + 1]
                + util_data[scen]["Gen RNA Positive"].iloc[:, i + 1]
            )  # acts as a a double check

            # Treatment initiations
            util_data[scen]["Total treatment"].iloc[:, i + 1] = data["tx_m", "Total"]

            # Vaccinations administered
            util_data[scen]["PWID vaccinations"].iloc[:, i + 1] = data[
                "vaccine_uptake:flow", "pwid"
            ]
            util_data[scen]["Prisoner vaccinations"].iloc[:, i + 1] = data[
                "vaccine_uptake:flow", "prisoners"
            ]
            util_data[scen]["Gen vaccinations"].iloc[:, i + 1] = (
                data["vaccine_uptake:flow", "children"]
                + data["vaccine_uptake:flow", "teens"] * 1 / 8
                + data["vaccine_uptake:flow", "adults"]
            )
            util_data[scen]["Total vaccinations"].iloc[:, i + 1] = (
                util_data[scen]["PWID vaccinations"].iloc[:, i + 1]
                + util_data[scen]["Prisoner vaccinations"].iloc[:, i + 1]
                + util_data[scen]["Gen vaccinations"].iloc[:, i + 1]
            )

    # Disease Stages (can be used for productivity, DALYs and management costs)
    dis_data = {}
    for scen in scen_name:
        dis_data[scen] = {}

    dis_outs = [
        "F0-F3 infected",
        "F4+ infected",
        "F0-F3 cured",
        "F4+ cured",
        "F0-F2 diagnosed",
        "F3 diagnosed",
        "F4 diagnosed",
        "DC diagnosed",
        "HCC diagnosed",
        "DC undiagnosed",
        "HCC undiagnosed",
        "Total F4",
        "Total Decompensated",
        "Total Liver Cancer",
    ]

    for scen in scen_name:
        for out in dis_outs:
            dis_data[scen][out] = pd.DataFrame(columns=colnames)
            dis_data[scen][out]["year"] = np.arange(2000.5, 2051.5, 1)

    for scen in scen_name:
        for i, sim in enumerate(sim_loop):
            data = sc.load(savedir_scens / f"{scen}_{sim}_extracted.pkl")
            # Productivity
            dis_data[scen]["F0-F3 infected"].iloc[:, i + 1] = data[
                "f0f3_utx", "working_age"
            ]
            dis_data[scen]["F4+ infected"].iloc[:, i + 1] = data[
                "f4_utx", "working_age"
            ]
            dis_data[scen]["F0-F3 cured"].iloc[:, i + 1] = data["f0f3_c", "working_age"]
            dis_data[scen]["F4+ cured"].iloc[:, i + 1] = data["f4_c", "working_age"]

            # Disease Management Costs
            dis_data[scen]["F0-F2 diagnosed"].iloc[:, i + 1] = data["f0f2_dx", "Total"]
            dis_data[scen]["F3 diagnosed"].iloc[:, i + 1] = data["f3_dx", "Total"]
            dis_data[scen]["F4 diagnosed"].iloc[:, i + 1] = data["f4_dx", "Total"]
            dis_data[scen]["DC diagnosed"].iloc[:, i + 1] = data["dc_dx", "Total"]
            dis_data[scen]["HCC diagnosed"].iloc[:, i + 1] = data["hcc_dx", "Total"]
            dis_data[scen]["DC undiagnosed"].iloc[:, i + 1] = data["dc_udx", "Total"]
            dis_data[scen]["HCC undiagnosed"].iloc[:, i + 1] = data["hcc_udx", "Total"]

            # DALYs
            dis_data[scen]["Total F4"].iloc[:, i + 1] = data["f4_all", "Total"]
            dis_data[scen]["Total Decompensated"].iloc[:, i + 1] = data[
                "d_cirrhosis_all", "Total"
            ]
            dis_data[scen]["Total Liver Cancer"].iloc[:, i + 1] = data[
                "cancer_all", "Total"
            ]

    # Call data for economic analysis
    global_pars = pd.read_excel(
        str(rootdir) + "/data/flat_datasheet.xlsx",
        sheet_name="Cost - Global Inputs",
    )
    dir_costs = pd.read_excel(
        str(rootdir) + "/data/flat_datasheet.xlsx",
        sheet_name="Cost - Direct Costs",
    )
    dir_costs = dir_costs[dir_costs["ISO3"] == country]
    prod_costs = pd.read_excel(
        str(rootdir) + "/data/flat_datasheet.xlsx",
        sheet_name="Cost - YLL and productivity",
    )
    prod_costs = prod_costs[prod_costs["ISO3"] == country]

    # Not Sampled Economic Pars
    start_year = global_pars.iloc[10, 1]
    c_disc, h_disc = global_pars.iloc[9, 1], global_pars.iloc[8, 1]
    vax_weight = dir_costs.vax_weight.values[0]
    abs_f, abs_fr, pres_f, pres_fr, abs_c, abs_cr, pres_c, pres_cr = global_pars[
        "value"
    ].iloc[0:8]
    waste, logistics, overheads, prog_costs = global_pars["value"].iloc[17:21]
    f0f3_care, f4dx_care, dchcc_care = (
        dir_costs.f0f3_dx_care.values[0],
        dir_costs.f4_dx_care.values[0],
        dir_costs.dchcc_care.values[0],
    )
    employed, gdp, hr_salary = (
        prod_costs.ER_total.values[0],
        prod_costs.GDP.values[0],
        prod_costs.hour_wage.values[0],
    )

    # Sampled Economic Pars (central value, n_sims sampled in array)
    np.random.seed(rand_seed)
    f0f2_cost = [dir_costs.f0f2.values[0]] + list(
        np.random.triangular(
            dir_costs.f0f2_lb.values[0],
            dir_costs.f0f2.values[0],
            dir_costs.f0f2_ub.values[0],
            n_samples,
        )
    )  # triangular
    f3_cost = [dir_costs.f3.values[0]] + list(
        np.random.triangular(
            dir_costs.f3_lb.values[0],
            dir_costs.f3.values[0],
            dir_costs.f3_ub.values[0],
            n_samples,
        )
    )  # triangular
    f4_cost = [dir_costs.f4.values[0]] + list(
        np.random.triangular(
            dir_costs.f4_lb.values[0],
            dir_costs.f4.values[0],
            dir_costs.f4_ub.values[0],
            n_samples,
        )
    )  # triangular
    dc_cost = [dir_costs.DC.values[0]] + list(
        np.random.triangular(
            dir_costs.DC_lb.values[0],
            dir_costs.DC.values[0],
            dir_costs.DC_ub.values[0],
            n_samples,
        )
    )  # triangular
    hcc_cost = [dir_costs.HCC.values[0]] + list(
        np.random.triangular(
            dir_costs.HCC_lb.values[0],
            dir_costs.HCC.values[0],
            dir_costs.HCC_ub.values[0],
            n_samples,
        )
    )  # triangular
    ab_test = [global_pars["value"].iloc[14]] + list(
        np.random.uniform(
            global_pars["lower"].iloc[14], global_pars["upper"].iloc[14], n_samples
        )
    )  # uniform
    rna_test = [global_pars["value"].iloc[15]] + list(
        np.random.uniform(
            global_pars["lower"].iloc[15], global_pars["upper"].iloc[15], n_samples
        )
    )  # uniform
    daa = [global_pars["value"].iloc[16]] + list(
        np.random.uniform(
            global_pars["lower"].iloc[16], global_pars["upper"].iloc[16], n_samples
        )
    )  # uniform
    f4_weight = [global_pars["value"].iloc[21]] + list(
        np.random.triangular(
            global_pars["lower"].iloc[21],
            global_pars["value"].iloc[21],
            global_pars["upper"].iloc[21],
            n_samples,
        )
    )  # triangular
    dc_weight = [global_pars["value"].iloc[11]] + list(
        np.random.triangular(
            global_pars["lower"].iloc[11],
            global_pars["value"].iloc[11],
            global_pars["upper"].iloc[11],
            n_samples,
        )
    )  # triangular
    hcc_weight = [global_pars["value"].iloc[12]] + list(
        np.random.triangular(
            global_pars["lower"].iloc[12],
            global_pars["value"].iloc[12],
            global_pars["upper"].iloc[12],
            n_samples,
        )
    )  # triangular

    # Discount Arrays
    cost_disc, health_disc = np.zeros((len(np.arange(2000, 2051, 1)), 2)), np.zeros(
        (len(np.arange(2000, 2051, 1)), 2)
    )
    cost_disc[:, 0], health_disc[:, 0] = np.arange(2000.5, 2051.5, 1), np.arange(
        2000.5, 2051.5, 1
    )
    for i in range(len(cost_disc[:, 0])):
        if cost_disc[i, 0] < start_year - 0.5:
            cost_disc[i, 1], health_disc[i, 1] = 0, 1
        else:
            cost_disc[i, 1] = (1 - c_disc) ** (cost_disc[i, 0] - (start_year - 0.5))
            health_disc[i, 1] = (1 - h_disc) ** (cost_disc[i, 0] - (start_year - 0.5))

    # Calculate YLL and YPLL (counterfactual method)
    aging = [
        1 / 15,
        1 / 15,
        1 / 20,
        1 / 15,
    ]  # Age bracket length (65+ retire and no longer contribute)
    death_dist = [
        prod_costs["014_prp"].values[0],
        prod_costs["1529_prp"].values[0],
        prod_costs["3049_prp"].values[0],
        prod_costs["5059_prp"].values[0] + 0.5 * prod_costs["6069_prp"].values[0],
        0.5 * prod_costs["6069_prp"].values[0] + prod_costs["70_prp"].values[0],
    ]
    acm = [
        prod_costs["mrate_014"].values[0],
        prod_costs["mrate_1529"].values[0],
        prod_costs["mrate_3049"].values[0],
        prod_costs["mrate_5064"].values[0],
        1 / prod_costs["le_70"].values[0],
    ]

    ghosts = {}
    for scen in scen_name:
        ghosts[scen] = {}
        for i, sim in enumerate(sim_loop):
            ghosts[scen][f"sim_{sim}"] = np.zeros(
                (len(agg_data[scen]["HCV mortality"]), len(acm))
            )

    for scen in scen_name:
        for i, sim in enumerate(sim_loop):
            for t in range(1, len(agg_data[scen]["HCV mortality"])):
                ppl_who_age = ghosts[scen][f"sim_{sim}"][t - 1, 0 : len(aging)] * aging
                ghosts[scen][f"sim_{sim}"][t, 0] = max(
                    0,
                    ghosts[scen][f"sim_{sim}"][t - 1, 0]
                    - ppl_who_age[0]
                    - acm[0] * ghosts[scen][f"sim_{sim}"][t - 1, 0],
                )  # 0-14
                ghosts[scen][f"sim_{sim}"][t, 1] = max(
                    0,
                    ghosts[scen][f"sim_{sim}"][t - 1, 1]
                    - ppl_who_age[1]
                    + ppl_who_age[0]
                    - acm[1] * ghosts[scen][f"sim_{sim}"][t - 1, 1],
                )  # 15-29
                ghosts[scen][f"sim_{sim}"][t, 2] = max(
                    0,
                    ghosts[scen][f"sim_{sim}"][t - 1, 2]
                    - ppl_who_age[2]
                    + ppl_who_age[1]
                    - acm[2] * ghosts[scen][f"sim_{sim}"][t - 1, 2],
                )  # 30-49
                ghosts[scen][f"sim_{sim}"][t, 3] = max(
                    0,
                    ghosts[scen][f"sim_{sim}"][t - 1, 3]
                    - ppl_who_age[3]
                    + ppl_who_age[2]
                    - acm[3] * ghosts[scen][f"sim_{sim}"][t - 1, 3],
                )  # 50-64
                ghosts[scen][f"sim_{sim}"][t, 4] = max(
                    0,
                    ghosts[scen][f"sim_{sim}"][t - 1, 4]
                    + ppl_who_age[3]
                    - acm[4] * ghosts[scen][f"sim_{sim}"][t - 1, 4],
                )  # 65+

                for dage in range(len(death_dist)):
                    ghosts[scen][f"sim_{sim}"][t, dage] = (
                        ghosts[scen][f"sim_{sim}"][t, dage]
                        + agg_data[scen]["HCV mortality"].iloc[t, i + 1]
                        * death_dist[dage]
                    )

    # Collate outcomes from counterfactuals
    lives_lost = {}
    for scen in scen_name:
        lives_lost[scen] = {}

    lives_outs = ["Years Life Lost", "Years Productive Life Lost"]

    for scen in scen_name:
        for out in lives_outs:
            lives_lost[scen][out] = pd.DataFrame(columns=colnames)
            lives_lost[scen][out]["year"] = np.arange(2000.5, 2051.5, 1)

    for scen in scen_name:
        for i, sim in enumerate(sim_loop):
            lives_lost[scen]["Years Life Lost"].iloc[:, i + 1] = np.sum(
                ghosts[scen][f"sim_{sim}"], axis=1
            )
            lives_lost[scen]["Years Productive Life Lost"].iloc[:, i + 1] = np.sum(
                ghosts[scen][f"sim_{sim}"][:, 1:4], axis=1
            )

    # Calculate costs
    econ_ests = {}
    for scen in scen_name:
        econ_ests[scen] = {}

    econ_outs = [
        "Diagnosis Cost",
        "Disease Management Cost",
        "Treatment Cost",
        "Morbidity Productivity",
        "Mortality Productivity",
        "Total Productivity",
        "Total Cost",
        "Vaccine_5USD",
        "Vaccine_10USD",
        "Vaccine_15USD",
        "Vaccine_20USD",
        "YLD",
        "DALYs",
    ]

    for scen in scen_name:
        for out in econ_outs:
            econ_ests[scen][out] = pd.DataFrame(columns=colnames)
            econ_ests[scen][out]["year"] = np.arange(2000.5, 2051.5, 1)

    for scen in scen_name:
        for i, sim in enumerate(sim_loop):
            for t in range(len(econ_ests[scen]["Diagnosis Cost"]["year"])):
                econ_ests[scen]["Diagnosis Cost"].iloc[t, i + 1] = (
                    (
                        util_data[scen]["PWID Ab Tests"].iloc[t, i + 1]
                        + util_data[scen]["Prisoner Ab Tests"].iloc[t, i + 1]
                        + util_data[scen]["Gen Ab Tests"].iloc[t, i + 1]
                    )
                    * (
                        (ab_test[i] * (1 + waste + logistics + prog_costs))
                        + (2 * (hr_salary * (1 + overheads + prog_costs)))
                    )
                    + (
                        util_data[scen]["PWID RNA Tests"].iloc[t, i + 1]
                        + util_data[scen]["Prisoner RNA Tests"].iloc[t, i + 1]
                        + util_data[scen]["Gen RNA Tests"].iloc[t, i + 1]
                    )
                    * (
                        (rna_test[i] * (1 + waste + logistics + prog_costs))
                        + (2 * (hr_salary * (1 + overheads + prog_costs)))
                    )
                ) * cost_disc[t, 1]

                econ_ests[scen]["Disease Management Cost"].iloc[t, i + 1] = (
                    (
                        dis_data[scen]["F0-F2 diagnosed"].iloc[t, i + 1]
                        * f0f2_cost[i]
                        * f0f3_care
                    )
                    + (
                        dis_data[scen]["F3 diagnosed"].iloc[t, i + 1]
                        * f3_cost[i]
                        * f0f3_care
                    )
                    + (
                        dis_data[scen]["F4 diagnosed"].iloc[t, i + 1]
                        * f4_cost[i]
                        * f4dx_care
                    )
                    + (
                        (
                            dis_data[scen]["DC diagnosed"].iloc[t, i + 1]
                            + dis_data[scen]["DC undiagnosed"].iloc[t, i + 1]
                        )
                        * dc_cost[i]
                        * dchcc_care
                    )
                    + (
                        (
                            dis_data[scen]["HCC undiagnosed"].iloc[t, i + 1]
                            + dis_data[scen]["HCC diagnosed"].iloc[t, i + 1]
                        )
                        * hcc_cost[i]
                        * dchcc_care
                    )
                ) * cost_disc[t, 1]

                econ_ests[scen]["Treatment Cost"].iloc[t, i + 1] = (
                    util_data[scen]["Total treatment"].iloc[t, i + 1]
                    * (
                        daa[i] * (1 + waste + logistics + prog_costs)
                        + hr_salary * (1 + overheads + prog_costs)
                    )
                    * cost_disc[t, 1]
                )

                econ_ests[scen]["Morbidity Productivity"].iloc[t, i + 1] = (
                    (gdp / employed)
                    * employed
                    * (
                        dis_data[scen]["F0-F3 infected"].iloc[t, i + 1]
                        * (abs_f + pres_f)
                        + dis_data[scen]["F0-F3 cured"].iloc[t, i + 1]
                        * (abs_fr + pres_fr)
                        + dis_data[scen]["F4+ infected"].iloc[t, i + 1]
                        * (abs_c + pres_c)
                        + dis_data[scen]["F4+ cured"].iloc[t, i + 1]
                        * (abs_cr + pres_cr)
                    )
                    * cost_disc[t, 1]
                )

                econ_ests[scen]["Mortality Productivity"].iloc[t, i + 1] = (
                    (gdp / employed)
                    * employed
                    * lives_lost[scen]["Years Productive Life Lost"].iloc[t, i + 1]
                    * cost_disc[t, 1]
                )

                econ_ests[scen]["Total Productivity"].iloc[t, i + 1] = (
                    econ_ests[scen]["Morbidity Productivity"].iloc[t, i + 1]
                    + econ_ests[scen]["Mortality Productivity"].iloc[t, i + 1]
                )

                econ_ests[scen]["Total Cost"].iloc[t, i + 1] = (
                    econ_ests[scen]["Diagnosis Cost"].iloc[t, i + 1]
                    + econ_ests[scen]["Disease Management Cost"].iloc[t, i + 1]
                    + econ_ests[scen]["Treatment Cost"].iloc[t, i + 1]
                    + econ_ests[scen]["Total Productivity"].iloc[t, i + 1]
                )

                econ_ests[scen]["Vaccine_5USD"].iloc[t, i + 1] = (
                    (
                        util_data[scen]["PWID vaccinations"].iloc[t, i + 1]
                        + util_data[scen]["Prisoner vaccinations"].iloc[t, i + 1]
                        + util_data[scen]["Gen vaccinations"].iloc[t, i + 1]
                    )
                    * (
                        5 * vax_weight * (1 + waste + logistics + prog_costs)
                        + (hr_salary * (1 / 12)) * (1 + overheads + prog_costs)
                    )
                    * cost_disc[t, 1]
                )

                econ_ests[scen]["Vaccine_10USD"].iloc[t, i + 1] = (
                    (
                        util_data[scen]["PWID vaccinations"].iloc[t, i + 1]
                        + util_data[scen]["Prisoner vaccinations"].iloc[t, i + 1]
                        + util_data[scen]["Gen vaccinations"].iloc[t, i + 1]
                    )
                    * (
                        10 * vax_weight * (1 + waste + logistics + prog_costs)
                        + (hr_salary * (1 / 12)) * (1 + overheads + prog_costs)
                    )
                    * cost_disc[t, 1]
                )

                econ_ests[scen]["Vaccine_15USD"].iloc[t, i + 1] = (
                    (
                        util_data[scen]["PWID vaccinations"].iloc[t, i + 1]
                        + util_data[scen]["Prisoner vaccinations"].iloc[t, i + 1]
                        + util_data[scen]["Gen vaccinations"].iloc[t, i + 1]
                    )
                    * (
                        15 * vax_weight * (1 + waste + logistics + prog_costs)
                        + (hr_salary * (1 / 12)) * (1 + overheads + prog_costs)
                    )
                    * cost_disc[t, 1]
                )

                econ_ests[scen]["Vaccine_20USD"].iloc[t, i + 1] = (
                    (
                        util_data[scen]["PWID vaccinations"].iloc[t, i + 1]
                        + util_data[scen]["Prisoner vaccinations"].iloc[t, i + 1]
                        + util_data[scen]["Gen vaccinations"].iloc[t, i + 1]
                    )
                    * (
                        20 * vax_weight * (1 + waste + logistics + prog_costs)
                        + (hr_salary * (1 / 12)) * (1 + overheads + prog_costs)
                    )
                    * cost_disc[t, 1]
                )

                econ_ests[scen]["YLD"].iloc[t, i + 1] = (
                    dis_data[scen]["Total Decompensated"].iloc[t, i + 1] * dc_weight[i]
                    + dis_data[scen]["Total Liver Cancer"].iloc[t, i + 1]
                    * hcc_weight[i]
                    + dis_data[scen]["Total F4"].iloc[t, i + 1] * f4_weight[i]
                )

                econ_ests[scen]["DALYs"].iloc[t, i + 1] = (
                    econ_ests[scen]["YLD"].iloc[t, i + 1]
                    + lives_lost[scen]["Years Life Lost"].iloc[t, i + 1]
                ) * health_disc[t, 1]

    # Aggregate Dictionary (of dictionaries of dictionaries....)
    output_dict = {}
    output_dict["agg_data"] = agg_data
    output_dict["util_data"] = util_data
    output_dict["dis_data"] = dis_data
    output_dict["lives_lost"] = lives_lost
    output_dict["econ_ests"] = econ_ests

    # Save as pkl file
    filename = sc.makefilepath(
        filename=f"{country}_econ eval.pkl",
        folder=results_folder / "scenarios" / "agg_outputs",
        makedirs=True,
    )
    sc.save(filename=filename, obj=output_dict)


def aggregate_ensembles(scens_folder, n_samples):
    """Aggregate epidemiological and economic data from multiple country scenarios.
    
    This function reads scenario data from specified folders, aggregates epidemiological and economic outputs for different regions and scenarios, and saves the aggregated data into pickle files.
    
    Args:
        scens_folder (str or Path): The path to the folder containing scenario outputs.
        n_samples (int): The number of samples to aggregate for each output.
    
    Returns:
        None: The function saves the aggregated data to files and does not return any value.
    
    Raises:
        FileNotFoundError: If the specified scenario folder or data files do not exist.
        ValueError: If the data format is not as expected.
    
    Notes:
        - The function expects specific Excel sheets and structures in the input files.
        - The aggregation is performed for various epidemiological outputs such as incidence and mortality, as well as economic outputs like costs and treatment statistics.
    """

    ref_countries = pd.read_excel(
        str(rootdir) + "/data/flat_datasheet.xlsx",
        sheet_name="Cost - YLL and productivity",
    ).iloc[:, np.r_[1, 4]]
    who_reg = list(pd.unique(ref_countries.WHO_reg)) + ["global", "top10"]
    loop_folder = pathlib.Path(str(scens_folder) + "/agg_outputs/")
    countries = [
        country.stem[:3] for country in loop_folder.iterdir() if country.is_file()
    ]

    scenarios = pd.read_excel(
        str(rootdir) + "/data/progbook_inputs.xlsx", sheet_name="scenarios"
    )
    scen_name = list(pd.unique(scenarios.scenario_name))

    # Produce lists of countries within each region, this can be adjusted pretty easily now
    ref_countries = ref_countries[ref_countries["ISO3"].isin(countries)].reset_index(
        drop=True
    )
    top_10_burden = [
        "PAK",
        "IND",
        "CHN",
        "RUS",
        "USA",
        "IDN",
        "NGA",
        "UKR",
        "UZB",
        "BGD",
    ]

    AFR, AMR, EMR, EUR, SEAR, WPR, TOP = [], [], [], [], [], [], []
    for i in range(len(ref_countries)):
        if ref_countries.iloc[i, 1] == "AFR":
            AFR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 1] == "AMR":
            AMR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 1] == "EMR":
            EMR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 1] == "EUR":
            EUR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 1] == "SEAR":
            SEAR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 1] == "WPR":
            WPR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 0] in top_10_burden:
            TOP.append(ref_countries.iloc[i, 0])

    # Produce empty dataframes to collate results (fillna means can do iterative sums)
    epis_outs = ["HCV incidence", "HCC incidence", "HCV mortality", "DALYs"]

    epi_agg = {}
    for reg in who_reg:
        epi_agg[reg] = {}
        for scen in scen_name:
            epi_agg[reg][scen] = {}
            for out in epis_outs:
                epi_agg[reg][scen][out] = pd.DataFrame(
                    columns=["year", "central"] + [f"run_{i}" for i in range(n_samples)]
                )
                epi_agg[reg][scen][out].year = np.arange(2000.5, 2051.1, 1)
                epi_agg[reg][scen][out] = epi_agg[reg][scen][out].fillna(0)

    # Economic aggregates
    econ_outs = [
        "PWID Ab Tests",
        "PWID RNA Tests",
        "Prisoner Ab Tests",
        "Prisoner RNA Tests",
        "Gen Ab Tests",
        "Gen RNA Tests",
        "Total Ab Tests",
        "Total RNA Tests",
        "PWID Ab Positive",
        "PWID RNA Positive",
        "Prisoner Ab Positive",
        "Prisoner RNA Positive",
        "Gen Ab Positive",
        "Gen RNA Positive",
        "Total Ab Positive",
        "Total RNA Positive",
        "Total treatment",
        "PWID vaccinations",
        "Prisoner vaccinations",
        "Gen vaccinations",
        "Total vaccinations",
        "Diagnosis Cost",
        "Disease Management Cost",
        "Treatment Cost",
        "Morbidity Cost",
        "Mortality Cost",
        "Direct Cost",
        "Productivity Cost",
        "Total Cost",
        "Vaccine 5USD Cost",
        "Vaccine 10USD Cost",
        "Vaccine 15USD Cost",
        "Vaccine 20USD Cost",
    ]

    econ_agg = {}
    for reg in who_reg:
        econ_agg[reg] = {}
        for scen in scen_name:
            econ_agg[reg][scen] = {}
            for out in econ_outs:
                econ_agg[reg][scen][out] = pd.DataFrame(
                    columns=["year", "central"] + [f"run_{i}" for i in range(n_samples)]
                )
                econ_agg[reg][scen][out].year = np.arange(2000.5, 2051.1, 1)
                econ_agg[reg][scen][out] = econ_agg[reg][scen][out].fillna(0)

    for country in countries:
        country_data = sc.load(scens_folder / f"agg_outputs/{country}_econ_eval.pkl")

        if country in AFR:
            adds = ["AFR", "global"]
        if country in AMR:
            adds = ["AMR", "global"]
        if country in EMR:
            adds = ["EMR", "global"]
        if country in EUR:
            adds = ["EUR", "global"]
        if country in SEAR:
            adds = ["SEAR", "global"]
        if country in WPR:
            adds = ["WPR", "global"]
        if country in TOP:
            adds += ["top10"]

        for add in adds:
            for scen in scen_name:
                for i in range(n_samples + 1):
                    # Epidemiological Outs
                    epi_agg[add][scen]["HCV incidence"].iloc[:, i + 1] = (
                        epi_agg[add][scen]["HCV incidence"].iloc[:, i + 1]
                        + country_data["agg_data"][scen]["HCV incidence"].iloc[:, i + 1]
                    )
                    epi_agg[add][scen]["HCC incidence"].iloc[:, i + 1] = (
                        epi_agg[add][scen]["HCC incidence"].iloc[:, i + 1]
                        + country_data["agg_data"][scen]["HCC incidence"].iloc[:, i + 1]
                    )
                    epi_agg[add][scen]["HCV mortality"].iloc[:, i + 1] = (
                        epi_agg[add][scen]["HCV mortality"].iloc[:, i + 1]
                        + country_data["agg_data"][scen]["HCV mortality"].iloc[:, i + 1]
                    )
                    epi_agg[add][scen]["DALYs"].iloc[:, i + 1] = (
                        epi_agg[add][scen]["DALYs"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["DALYs"].iloc[:, i + 1]
                    )

                    # Tests completed
                    econ_agg[add][scen]["PWID Ab Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID Ab Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["PWID Ab Tests"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["PWID RNA Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID RNA Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["PWID RNA Tests"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Prisoner Ab Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Prisoner Ab Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Prisoner Ab Tests"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Prisoner RNA Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Prisoner RNA Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Prisoner RNA Tests"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Gen Ab Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Gen Ab Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Gen Ab Tests"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Gen RNA Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Gen RNA Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Gen RNA Tests"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Total Ab Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID Ab Tests"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Prisoner Ab Tests"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Gen Ab Tests"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Total RNA Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID RNA Tests"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Prisoner RNA Tests"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Gen RNA Tests"].iloc[:, i + 1]
                    )

                    # Positive Tests
                    econ_agg[add][scen]["PWID Ab Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID Ab Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["PWID Ab Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["PWID RNA Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID RNA Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["PWID RNA Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Prisoner Ab Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Prisoner Ab Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Prisoner Ab Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Prisoner RNA Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Prisoner RNA Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Prisoner RNA Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Gen Ab Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Gen Ab Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Gen Ab Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Gen RNA Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Gen RNA Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Gen RNA Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Total Ab Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID Ab Positive"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Prisoner Ab Positive"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Gen Ab Positive"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Total RNA Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID RNA Positive"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Prisoner RNA Positive"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Gen RNA Positive"].iloc[:, i + 1]
                    )

                    # Vaccinations and Treatments Administered
                    econ_agg[add][scen]["Total treatment"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Total treatment"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Total treatment"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["PWID vaccinations"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID vaccinations"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["PWID vaccinations"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Prisoner vaccinations"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Prisoner vaccinations"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Prisoner vaccinations"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Gen vaccinations"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Gen vaccinations"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Gen vaccinations"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Total vaccinations"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID vaccinations"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Prisoner vaccinations"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Gen vaccinations"].iloc[:, i + 1]
                    )

                    # Costs
                    econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["Diagnosis Cost"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen][
                            "Disease Management Cost"
                        ].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["Treatment Cost"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Morbidity Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Morbidity Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen][
                            "Morbidity Productivity"
                        ].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Mortality Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Mortality Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen][
                            "Mortality Productivity"
                        ].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Direct Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Productivity Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Morbidity Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Mortality Cost"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Total Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Productivity Cost"].iloc[:, i + 1]
                    )

                    # Vaccination Cost
                    econ_agg[add][scen]["Vaccine 5USD Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Vaccine 5USD Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["Vaccine_5USD"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Vaccine 10USD Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Vaccine 10USD Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["Vaccine_10USD"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Vaccine 15USD Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Vaccine 15USD Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["Vaccine_15USD"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Vaccine 20USD Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Vaccine 20USD Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["Vaccine_20USD"].iloc[
                            :, i + 1
                        ]
                    )

    filename = sc.makefilepath(
        filename="epi_agg.pkl", folder=scens_folder, makedirs=True
    )
    sc.save(filename=filename, obj=epi_agg)

    filename = sc.makefilepath(
        filename="econ_agg.pkl", folder=scens_folder, makedirs=True
    )
    sc.save(filename=filename, obj=econ_agg)


def econ_analysis(scens_folder, n_samples):
    """Perform economic analysis based on provided scenarios and sample data.
    
    This function loads economic and epidemiological data from specified folders, processes the data to calculate aggregate outcomes, and saves the results into Excel and CSV files. It generates tables for regional aggregate outcomes and comparative outcomes, including cost-effectiveness metrics such as Benefit-Cost Ratios (BCR) and Incremental Cost-Effectiveness Ratios (ICER).
    
    Args:
        scens_folder (pathlib.Path): The folder path containing scenario data files.
        n_samples (int): The number of samples to be used in the analysis.
    
    Returns:
        None: The function saves the results to Excel and CSV files in the specified folder.
    """

    # Global and Regional Summary Table
    reg_data = sc.load(scens_folder / "econ_agg.pkl")
    scenarios = pd.read_excel(
        str(rootdir) + "/data/progbook_inputs.xlsx", sheet_name="scenarios"
    )
    scen_name = list(
        pd.unique(scenarios.scenario_name)
    )  # needed to loop through for data extraction etc

    scen_cfs = scenarios.iloc[:, np.r_[0, 1, 5]].dropna().reset_index(drop=True)
    scens = list(scen_cfs.scenario_name)
    comps = list(scen_cfs.counterfactual)
    scen_label = list(pd.unique(scenarios.description))  # Comparators and Scenarios
    vax_label = scen_label[3:]  # remove counterfactuals
    vax_label_scen = [f"Scenario {i}" for i in range(1, 6)]
    regions = ["AFR", "AMR", "EMR", "EUR", "SEAR", "WPR", "global", "top10"]

    # # Supplemental Figure - Number of tests, treatments, and vaccinations per scenario (global)
    # measures = ["Total Ab Tests", "Total Ab Positive", "Total RNA Tests", "Total RNA Positive","Total treatment", "Total vaccinations"]
    # fig1 = plt.figure(figsize = (20, 20))
    #
    # for idx, out in enumerate(measures):
    #     ax = fig1.add_subplot(2,3, idx+1)
    #     for j, scen in enumerate(scen_name):
    #         year = reg_data["global"][scen][out].iloc[1:,0].astype(float)
    #         lower = np.percentile(reg_data["global"][scen][out].iloc[1:,2:n_samples+1], 2.5, axis=1).astype(float)
    #         upper = np.percentile(reg_data["global"][scen][out].iloc[1:,2:n_samples+1], 97.5, axis=1).astype(float)
    #         ax.plot(reg_data["global"][scen][out]["year"], reg_data["global"][scen][out]["central"], alpha=0.6, label = scen_label[j])
    #         ax.fill_between(year, lower, upper, alpha=0.2)
    #         ax.set_title(out)
    #         ax.set_xlim(2015.5, 2050.5)
    #         ax.legend(loc="best")
    #
    # fig1.tight_layout()

    # Aggregate (sum) outcomes (by region)
    agg_reg_table = {}
    outcomes_epi = ["HCV incidence", "HCC incidence", "HCV mortality", "DALYs"]
    outcomes_econ = [
        "Total Ab Tests",
        "Total RNA Tests",
        "Total treatment",
        "Total vaccinations",
        "Diagnosis Cost",
        "Disease Management Cost",
        "Treatment Cost",
        "Morbidity Cost",
        "Mortality Cost",
        "Direct Cost",
        "Productivity Cost",
        "Total Cost",
        "Vaccine 5USD Cost",
        "Vaccine 10USD Cost",
        "Vaccine 15USD Cost",
        "Vaccine 20USD Cost",
    ]
    outcomes = outcomes_epi + outcomes_econ
    epi_agg = sc.load(str(rootdir) + "results/scenarios/epi_agg.pkl")
    econ_agg = sc.load(str(rootdir) + "results/scenarios/econ_agg.pkl")

    for reg in regions:
        agg_reg_table[reg] = pd.DataFrame(columns=["scenarios"] + outcomes)
        agg_reg_table[reg].scenarios = scen_name
        agg_reg_table[reg] = agg_reg_table[reg].set_index("scenarios", drop=True)

        for scen in scen_name:
            for epi in outcomes_epi:
                # Limit from 2026-2050
                cen_epi = np.round(
                    np.sum(epi_agg[reg][scen][epi].iloc[25:, 1], axis=0), 0
                )
                low_epi = np.round(
                    np.percentile(
                        np.sum(epi_agg[reg][scen][epi].iloc[25:, 2:], axis=0), 2.5
                    ),
                    0,
                )
                upp_epi = np.round(
                    np.percentile(
                        np.sum(epi_agg[reg][scen][epi].iloc[25:, 2:], axis=0), 97.5
                    ),
                    0,
                )
                agg_reg_table[reg].at[
                    scen, epi
                ] = f"{cen_epi:,.0f} \n ({low_epi:,.0f} to {upp_epi:,.0f})"
                # agg_reg_table[reg].at[scen, epi] = f"{cen_epi/1e6:,.2f} \n ({low_epi/1e6:,.2f} to {upp_epi/1e6:,.2f})"

            for econ in outcomes_econ:
                cen_econ = np.sum(econ_agg[reg][scen][econ].iloc[:, 1], axis=0)
                low_econ = np.percentile(
                    np.sum(econ_agg[reg][scen][econ].iloc[25:, 2:], axis=0), 2.5
                )
                upp_econ = np.percentile(
                    np.sum(econ_agg[reg][scen][econ].iloc[25:, 2:], axis=0), 97.5
                )
                # if econ in ["Total vaccinations","Vaccine 5USD Cost", "Vaccine 10USD Cost", "Vaccine 15USD Cost", "Vaccine 20USD Cost"]:
                #     agg_reg_table[reg].at[scen, econ] = f"{cen_econ/1e6:,.2f} \n ({low_econ/1e6:,.2f} to {upp_econ/1e6:,.2f})"
                # else:
                #     agg_reg_table[reg].at[scen, econ] = f"{cen_econ/1e9:,.2f} \n ({low_econ/1e9:,.2f} to {upp_econ/1e9:,.2f})"
                agg_reg_table[reg].at[
                    scen, econ
                ] = f"{cen_econ:,.0f} \n ({low_econ:,.0f} to {upp_econ:,.0f})"

    # Vaccine Attributable Differences (by region)
    vad_reg_table = {}
    ben_cost = ["BCR_5USD", "BCR_10USD", "BCR_15USD", "BCR_20USD"]
    icer_cost = ["ICER_5USD", "ICER_10USD", "ICER_15USD", "ICER_20USD"]
    bcr_denom = [
        "Vaccine 5USD Cost",
        "Vaccine 10USD Cost",
        "Vaccine 15USD Cost",
        "Vaccine 20USD Cost",
    ]
    outcomes_vad = outcomes + ben_cost + icer_cost
    flip_list = [
        "Total vaccinations",
        "Vaccine 5USD Cost",
        "Vaccine 10USD Cost",
        "Vaccine 15USD Cost",
        "Vaccine 20USD Cost",
    ]

    for reg in regions:
        vad_reg_table[reg] = pd.DataFrame(columns=["scenarios"] + outcomes_vad)
        vad_reg_table[reg].scenarios = scens
        vad_reg_table[reg] = vad_reg_table[reg].set_index("scenarios", drop=True)

        for i, scen in enumerate(scens):
            for epi in outcomes_epi:
                cen_epi = np.sum(
                    epi_agg[reg][comps[i]][epi].iloc[25:, 1], axis=0
                ) - np.sum(epi_agg[reg][scen][epi].iloc[25:, 1], axis=0)
                low_epi = np.percentile(
                    np.sum(epi_agg[reg][comps[i]][epi].iloc[25:, 2:], axis=0)
                    - np.sum(epi_agg[reg][scen][epi].iloc[25:, 2:], axis=0),
                    2.5,
                )
                upp_epi = np.percentile(
                    np.sum(epi_agg[reg][comps[i]][epi].iloc[25:, 2:], axis=0)
                    - np.sum(epi_agg[reg][scen][epi].iloc[25:, 2:], axis=0),
                    97.5,
                )
                vad_reg_table[reg].at[
                    scen, epi
                ] = f"{cen_epi:,.0f} \n ({low_epi:,.0f} to {upp_epi:,.0f})"

            for econ in outcomes_econ:
                if econ in flip_list:
                    cen_econ = np.sum(
                        econ_agg[reg][scen][econ].iloc[25:, 1], axis=0
                    ) - np.sum(econ_agg[reg][comps[i]][econ].iloc[25:, 1], axis=0)
                    low_econ = np.percentile(
                        np.sum(econ_agg[reg][scen][econ].iloc[25:, 2:], axis=0)
                        - np.sum(econ_agg[reg][comps[i]][econ].iloc[25:, 2:], axis=0),
                        2.5,
                    )
                    upp_econ = np.percentile(
                        np.sum(econ_agg[reg][scen][econ].iloc[25:, 2:], axis=0)
                        - np.sum(econ_agg[reg][comps[i]][econ].iloc[25:, 2:], axis=0),
                        97.5,
                    )
                    vad_reg_table[reg].at[
                        scen, econ
                    ] = f"{cen_econ:,.0f} \n ({low_econ:,.0f} to {upp_econ:,.0f})"

                else:
                    cen_econ = np.sum(
                        econ_agg[reg][comps[i]][econ].iloc[25:, 1], axis=0
                    ) - np.sum(econ_agg[reg][scen][econ].iloc[25:, 1], axis=0)
                    low_econ = np.percentile(
                        np.sum(econ_agg[reg][comps[i]][econ].iloc[25:, 2:], axis=0)
                        - np.sum(econ_agg[reg][scen][econ].iloc[25:, 2:], axis=0),
                        2.5,
                    )
                    upp_econ = np.percentile(
                        np.sum(econ_agg[reg][comps[i]][econ].iloc[25:, 2:], axis=0)
                        - np.sum(econ_agg[reg][scen][econ].iloc[25:, 2:], axis=0),
                        97.5,
                    )
                    # if econ == "Productivity Cost":
                    #     vad_reg_table[reg].at[scen, econ] = f"{cen_econ/1e9:,.2f} \n ({low_econ/1e9:,.2f} to {upp_econ/1e9:,.2f})"
                    # elif econ == "Direct Cost":
                    #     vad_reg_table[reg].at[scen, econ] = f"{cen_econ/1e6:,.2f} \n ({low_econ/1e6:,.2f} to {upp_econ/1e6:,.2f})"
                    # else:
                    #     vad_reg_table[reg].at[scen, econ] = f"{cen_econ:,.0f} \n ({low_econ:,.0f} to {upp_econ:,.0f})"
                    vad_reg_table[reg].at[
                        scen, econ
                    ] = f"{cen_econ:,.0f} \n ({low_econ:,.0f} to {upp_econ:,.0f})"

            for j, bcr in enumerate(ben_cost):
                bcr_cent = (
                    np.sum(econ_agg[reg][comps[i]]["Total Cost"].iloc[25:, 1], axis=0)
                    - np.sum(econ_agg[reg][scen]["Total Cost"].iloc[25:, 1], axis=0)
                ) / np.sum(econ_agg[reg][scen][bcr_denom[j]].iloc[25:, 1], axis=0)
                bcr_low = np.percentile(
                    (
                        np.sum(
                            econ_agg[reg][comps[i]]["Total Cost"].iloc[25:, 2:], axis=0
                        )
                        - np.sum(
                            econ_agg[reg][scen]["Total Cost"].iloc[25:, 2:], axis=0
                        )
                    )
                    / np.sum(econ_agg[reg][scen][bcr_denom[j]].iloc[25:, 2:], axis=0),
                    2.5,
                )
                bcr_upp = np.percentile(
                    (
                        np.sum(
                            econ_agg[reg][comps[i]]["Total Cost"].iloc[25:, 2:], axis=0
                        )
                        - np.sum(
                            econ_agg[reg][scen]["Total Cost"].iloc[25:, 2:], axis=0
                        )
                    )
                    / np.sum(econ_agg[reg][scen][bcr_denom[j]].iloc[25:, 2:], axis=0),
                    97.5,
                )
                vad_reg_table[reg].at[
                    scen, bcr
                ] = f"{bcr_cent:,.2f} \n ({bcr_low:,.2f} to {bcr_upp:,.2f})"

            for k, icer in enumerate(icer_cost):
                icer_cent = (
                    (
                        np.sum(econ_agg[reg][scen]["Direct Cost"].iloc[25:, 1], axis=0)
                        + np.sum(econ_agg[reg][scen][bcr_denom[k]].iloc[25:, 1], axis=0)
                    )
                    - np.sum(
                        econ_agg[reg][comps[i]]["Direct Cost"].iloc[25:, 1], axis=0
                    )
                ) / (
                    np.sum(epi_agg[reg][comps[i]]["DALYs"].iloc[25:, 1], axis=0)
                    - np.sum(epi_agg[reg][scen]["DALYs"].iloc[25:, 1], axis=0)
                )  # Total Cost (Intervention) - Total Cost (Comparison) / DALYs (Comparsion) - DALYs (Intervention)
                icer_low = np.percentile(
                    (
                        (
                            np.sum(
                                econ_agg[reg][scen]["Direct Cost"].iloc[25:, 2:], axis=0
                            )
                            + np.sum(
                                econ_agg[reg][scen][bcr_denom[k]].iloc[25:, 2:], axis=0
                            )
                        )
                        - np.sum(
                            econ_agg[reg][comps[i]]["Direct Cost"].iloc[25:, 2:], axis=0
                        )
                    )
                    / (
                        np.sum(epi_agg[reg][comps[i]]["DALYs"].iloc[25:, 2:], axis=0)
                        - np.sum(epi_agg[reg][scen]["DALYs"].iloc[25:, 2:], axis=0)
                    ),
                    2.5,
                )
                icer_high = np.percentile(
                    (
                        (
                            np.sum(
                                econ_agg[reg][scen]["Direct Cost"].iloc[25:, 2:], axis=0
                            )
                            + np.sum(
                                econ_agg[reg][scen][bcr_denom[k]].iloc[25:, 2:], axis=0
                            )
                        )
                        - np.sum(
                            econ_agg[reg][comps[i]]["Direct Cost"].iloc[25:, 2:], axis=0
                        )
                    )
                    / (
                        np.sum(epi_agg[reg][comps[i]]["DALYs"].iloc[25:, 2:], axis=0)
                        - np.sum(epi_agg[reg][scen]["DALYs"].iloc[25:, 2:], axis=0)
                    ),
                    97.5,
                )
                vad_reg_table[reg].at[
                    scen, icer
                ] = f"{icer_cent:,.0f} \n ({icer_low:,.0f} to {icer_high:,.0f})"

    # Write to Excel spreadsheet and save
    with pd.ExcelWriter(scens_folder / "regional_aggregate_outcomes.xlsx") as writer:
        for sheet_name, df in agg_reg_table.items():
            df.to_excel(writer, sheet_name=sheet_name, index=True)
        print(
            "Aggregate outcome data table saved: {}".format(
                scens_folder / "regional_aggregate_outcomes.xlsx"
            )
        )

    with pd.ExcelWriter(scens_folder / "regional_comparative_outcomes.xlsx") as writer:
        for sheet_name, df in vad_reg_table.items():
            df.to_excel(writer, sheet_name=sheet_name, index=True)
        print(
            "Aggregate outcome data table saved: {}".format(
                scens_folder / "regional_comparative_outcomes.xlsx"
            )
        )

    # Country Level BCRs for mapping
    loop_folder = pathlib.Path(str(scens_folder) + "/agg_outputs/")
    countries = [
        country.stem[:3] for country in loop_folder.iterdir() if country.is_file()
    ]
    map_bcr = "Vaccine_5USD"

    country_bcr_map = pd.DataFrame(columns=["country"] + vax_label)
    country_bcr_map.country = countries
    country_bcr_map = country_bcr_map.set_index("country", drop=True)

    country_bcr = pd.DataFrame(columns=["country"] + vax_label_scen)
    country_bcr.country = [iso_to_country[country] for country in countries]
    country_bcr = country_bcr.set_index("country", drop=True)

    country_icer = pd.DataFrame(columns=["country"] + vax_label_scen)
    country_icer.country = [iso_to_country[country] for country in countries]
    country_icer = country_icer.set_index("country", drop=True)

    for country in countries:
        data = sc.load(scens_folder / f"agg_outputs/{country}_econ_eval.pkl")

        for i, scen in enumerate(vax_label):
            bcr_cent = (
                np.sum(data["econ_ests"][comps[i]]["Total Cost"].iloc[25:, 1])
                - np.sum(data["econ_ests"][scens[i]]["Total Cost"].iloc[25:, 1])
            ) / np.sum(data["econ_ests"][scens[i]][map_bcr].iloc[25:, 1])
            country_bcr_map.at[country, scen] = round(bcr_cent, 3)
            bcr_low = np.percentile(
                (
                    np.sum(
                        data["econ_ests"][comps[i]]["Total Cost"].iloc[25:, 2:], axis=0
                    )
                    - np.sum(
                        data["econ_ests"][scens[i]]["Total Cost"].iloc[25:, 2:], axis=0
                    )
                )
                / np.sum(data["econ_ests"][scens[i]][map_bcr].iloc[25:, 2:], axis=0),
                2.5,
            )
            bcr_upp = np.percentile(
                (
                    np.sum(
                        data["econ_ests"][comps[i]]["Total Cost"].iloc[25:, 2:], axis=0
                    )
                    - np.sum(
                        data["econ_ests"][scens[i]]["Total Cost"].iloc[25:, 2:], axis=0
                    )
                )
                / np.sum(data["econ_ests"][scens[i]][map_bcr], axis=0),
                97.5,
            )
            country_bcr.at[iso_to_country[country], vax_label_scen[i]] = (
                f"{bcr_cent:,.2f} \n ({bcr_low:,.2f} to {bcr_upp:,.2f})"
            )

            # econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i+1] + econ_agg[add][scen]["Disease Management Cost"].iloc[:, i+1] + econ_agg[add][scen]["Treatment Cost"].iloc[:, i+1]
            # +Vaccine 5USD Cost"
            num_1 = (
                np.sum(data["econ_ests"][scens[i]]["Diagnosis Cost"].iloc[25:, 1])
                + np.sum(
                    data["econ_ests"][scens[i]]["Disease Management Cost"].iloc[25:, 1]
                )
                + np.sum(data["econ_ests"][scens[i]]["Treatment Cost"].iloc[25:, 1])
                + np.sum(data["econ_ests"][scens[i]]["Vaccine_5USD"].iloc[25:, 1])
            )
            num_2 = (
                np.sum(data["econ_ests"][comps[i]]["Diagnosis Cost"].iloc[25:, 1])
                + np.sum(
                    data["econ_ests"][comps[i]]["Disease Management Cost"].iloc[25:, 1]
                )
                + np.sum(data["econ_ests"][comps[i]]["Treatment Cost"].iloc[25:, 1])
            )
            num = num_1 - num_2
            den = np.sum(data["econ_ests"][comps[i]]["DALYs"].iloc[25:, 1]) - np.sum(
                data["econ_ests"][scens[i]]["DALYs"].iloc[25:, 1]
            )
            icer_cent = num / den
            country_icer.at[iso_to_country[country], vax_label_scen[i]] = icer_cent
            # ((np.sum(econ_agg[reg][scen]["Direct Cost"].iloc[25:, 1], axis=0)+np.sum(econ_agg[reg][scen][bcr_denom[k]].iloc[25:, 1], axis=0)) - np.sum(econ_agg[reg][comps[i]]["Direct Cost"].iloc[25:, 1], axis=0)) / (np.sum(epi_agg[reg][comps[i]]["DALYs"].iloc[25:, 1], axis=0) - np.sum(epi_agg[reg][scen]["DALYs"].iloc[25:, 1], axis=0))   # Total Cost (Intervention) - Total Cost (Comparison) / DALYs (Comparsion) - DALYs (Intervention)
            # vad_reg_table[reg].at[scen, 'ICER'] = f"{icer_cent:,.0f}"

    country_bcr_map.to_csv(scens_folder / "bcr_map.csv")
    print("BCR map data saved: {}".format(scens_folder / "bcr_map.csv"))

    country_bcr.to_csv(scens_folder / "bcr_country.csv")
    print("BCR country data saved: {}".format(scens_folder / "bcr_country.csv"))

    country_icer.to_csv(scens_folder / "icer_country.csv")
    print("ICER country data saved: {}".format(scens_folder / "icer_country.csv"))


def write_print_table(scens_folder):
    """Writes and saves outcome data tables to an Excel file.
    
    This function reads scenario and outcome data from specified Excel sheets, processes the data for various regions, and writes the results to a new Excel file. The output includes two sheets: one for all outcomes and another for the benefit-cost ratios (BCRs) for the years 2026 to 2050.
    
    Args:
        scens_folder (Path): The folder path where the regional outcome Excel files are located and where the output Excel file will be saved.
    
    Returns:
        None: The function saves the output directly to an Excel file and does not return any value.
    
    Raises:
        FileNotFoundError: If the specified Excel files cannot be found in the given folder.
        ValueError: If the expected data structure in the Excel sheets does not match the function's requirements.
    """
    regions = ["global", "top10", "AFR", "AMR", "EMR", "EUR", "SEAR", "WPR"]
    regions_spelled = [
        "Global",
        "10 High Burden Countries",
        "African Region",
        "Region of the Americas",
        "Eastern Mediterranean Region",
        "European Region",
        "South East Asia Region",
        "Western Pacific Region",
    ]
    scenarios = pd.read_excel(
        str(rootdir) + "/data/progbook_inputs.xlsx", sheet_name="scenarios"
    )
    scen_name = list(pd.unique(scenarios.scenario_name))[3:]
    scen_desc = list(pd.unique(scenarios.description))[3:]
    # scenarios = {s:d for s,d in zip(scen_name,scen_desc)}
    outcomes = [
        ["Vaccines Administered", "Total vaccinations"],
        ["Vaccine-averted HCV infections", "HCV incidence"],
        ["Vaccine-averted HCC cases", "HCC incidence"],
        ["Vaccine-averted HCV deaths", "HCV mortality"],
        ["Vaccine-averted DALYs", "DALYs"],
        ["Vaccine-averted direct costs", "Direct Cost"],
        ["Vaccine-attributable productivity gains", "Productivity Cost"],
        ["Vaccine Costs", "Vaccine 5USD Cost"],
        ["Cost per DALY averted", "ICER_5USD"],
        ["Benefit-Cost Ratio", "BCR_5USD"],
    ]

    rows = []
    bcrs = []
    for reg, region in zip(regions, regions_spelled):
        comp = pd.read_excel(
            scens_folder / "regional_comparative_outcomes.xlsx", sheet_name=reg
        )
        agg = pd.read_excel(
            scens_folder / "regional_aggregate_outcomes.xlsx", sheet_name=reg
        )
        rows.append([region] + [""] * (len(scen_desc)))
        for out in outcomes:
            out_row = [out[0]]
            for scen in scen_name:
                data = (
                    comp[comp.scenarios == scen]
                    if out[0] != "Vaccine Costs"
                    else agg[agg.scenarios == scen]
                )
                out_row = out_row + [format_structured_string(data[out[1]].values[0])]
            rows.append(out_row)
    df_row = pd.DataFrame(rows, columns=["Outcomes\n2026–2050"] + scen_desc)

    for reg, region in zip(regions, regions_spelled):
        comp = pd.read_excel(
            scens_folder / "regional_comparative_outcomes.xlsx", sheet_name=reg
        )
        out_row = [region]
        for scen in scen_name:
            data = comp[comp.scenarios == scen]
            out_row = out_row + [data["BCR_5USD"].values[0]]
        bcrs.append(out_row)
    df_bcr = pd.DataFrame(bcrs, columns=["Benefit-Cost Ratio\n2026–2050"] + scen_desc)

    with pd.ExcelWriter(scens_folder / "print_tables.xlsx") as writer:
        df_row.to_excel(writer, sheet_name="all outcomes", index=False)
        df_bcr.to_excel(writer, sheet_name="bcrs", index=False)

    print(
        "Print outcome data tables saved: {}".format(scens_folder / "print_tables.xlsx")
    )


# %% Sensitivity analyses
def run_sensitivity_analyses(country, cal_folder, sens_folder, results_folder):
    """Run sensitivity analyses for a specified country and save the results.
    
    This function performs sensitivity analyses based on calibration and program data for a given country. It loads the necessary calibration and program data, defines various population groups, and computes outputs for different scenarios. The results are then saved to specified folders.
    
    Args:
        country (str): The name of the country for which the sensitivity analyses are to be run.
        cal_folder (str or Path): The folder path where calibration files are located.
        sens_folder (str or Path): The folder path where sensitivity analysis results will be saved.
        results_folder (str or Path): The folder path where program results are stored.
    
    Returns:
        None: The function saves the outputs to the specified directory and does not return any value.
    """
    P = project(
        country, load_calibration=True, cal_folder=cal_folder, load_programs=False
    )
    pops = P.data.pops.keys()  # define population

    # Function to define outputs to save
    def mapping_function(x):
        # Outputs to be saved
        outputs_by_pop = [
            "alive",
            "prevalence",
            "total_hcv",
            "inci_all_m",
            "notifications_m",
            "inci_m",
            "ab_tests_m",
            "pcr_tests_m",
        ]
        outputs_tot = [
            "total_hcv",
            "deaths_hcv_total",
            "tx_m",
            "inci_m",
            "notifications_m",
            "ab_tests_m",
            "pcr_tests_m",
            "f4_all",
            "d_cirrhosis_all",
            "inci_all_m",
            "hcc_inci",
            "acute_all",
            "undiag_all",
            "abpos_chronic",
            "pcr",
            "treated",
            "dc_udx",
            "hcc_udx",
            "f0f2_dx",
            "f3_dx",
            "f4_dx",
            "dc_dx",
            "hcc_dx",
            "d_cirrhosis_all",
            "cancer_all",
        ]
        outputs_tot_weighted = [
            "prop_f0f2",
            "prop_f3",
            "prop_f4",
            "prop_dc",
            "prop_hcc",
            "prevalence",
            "inci_100py",
        ]

        outputs_all = at.PlotData(
            x, outputs=outputs_by_pop, pops=pops, t_bins=1
        )  # outputs for each population group
        extra = at.PlotData(
            x, outputs=outputs_tot, pops="total", pop_aggregation="sum", t_bins=1
        )  # total outputs / summed
        outputs_all.series.extend(extra.series)
        extra = at.PlotData(
            x,
            outputs=outputs_tot_weighted,
            pops="total",
            pop_aggregation="weighted",
            t_bins=1,
        )  # total outputs / weighed
        outputs_all.series.extend(extra.series)

        # Aging pars
        transfer_pars = [
            "age_10-17_males_to_18-64_males",
            "age_10-17_females_to_18-64_females",
        ]
        transfer_pops = ["10-17_males", "10-17_females"]
        for par, pop in zip(transfer_pars, transfer_pops):
            extra = at.PlotData(x, outputs=par, pops=pop, t_bins=1)
            outputs_all.series.extend(extra.series)

        # Specific pop aggregations for econ fcts
        pars_tot = [
            "inci_all_m",
            "ab_all_m",
            "pcr_all_m",
            "f0f3_utx",
            "f4_utx",
            "f0f3_c",
            "f4_c",
            "inci_100py",
        ]
        pops_tot = [
            {
                "gen_pop": [
                    "0-9_males",
                    "0-9_females",
                    "10-17_males",
                    "10-17_females",
                    "18-64_males",
                    "18-64_females",
                    "65+_males",
                    "65+_females",
                ]
            },
            {"pwid": ["PWID_males", "PWID_females"]},
            {"prisoners": ["Prisoners_males", "Prisoners_females"]},
            {
                "working_age": [
                    "18-64_males",
                    "18-64_females",
                    "PWID_males",
                    "PWID_females",
                ]
            },
        ]
        for par in pars_tot:
            for pop in pops_tot:
                extra = at.PlotData(
                    x, outputs=par, pops=pop, pop_aggregation="sum", t_bins=1
                )
                outputs_all.series.extend(extra.series)

        pars_tot = ["rna_test_efficiency", "ab_test_efficiency", "prevalence"]
        extra = at.PlotData(
            x, outputs=pars_tot, pops=pops_tot, pop_aggregation="weighted", t_bins=1
        )
        outputs_all.series.extend(extra.series)

        pops_tot = [
            {"children": ["0-9_males", "0-9_females"]},
            {"teens": ["10-17_males", "10-17_females"]},
            {"adults": ["18-64_males", "18-64_females", "65+_males", "65+_females"]},
            {"pwid": ["PWID_males", "PWID_females"]},
            {"prisoners": ["Prisoners_males", "Prisoners_females"]},
        ]
        par = "vaccine_uptake:flow"
        extra = at.PlotData(
            x, outputs=par, pops=pops_tot, pop_aggregation="sum", t_bins=1
        )
        outputs_all.series.extend(extra.series)

        # Weighted spontaneous clearance rate
        tvec = (
            at.PlotData(x, outputs="acute", pops=pops, pop_aggregation="sum", t_bins=1)
            .series[0]
            .tvec
        )
        numerator = np.zeros_like(tvec)
        denominator = np.zeros_like(tvec)
        for pop in pops:
            acute = at.PlotData(x, outputs="acute", pops=pop, t_bins=1).series[0].vals
            acute_vax = (
                at.PlotData(x, outputs="acute_vax_reinfected", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            clear = (
                at.PlotData(x, outputs="prop_clear", pops=pop, t_bins=1).series[0].vals
            )
            clear_vax = (
                at.PlotData(x, outputs="vaccine_prop_clear", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            acute_vax_i = (
                at.PlotData(x, outputs="acute_vax_initial", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            clear_vax_i = (
                at.PlotData(x, outputs="vaccine_prop_clear_initial", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            numerator += (
                acute * clear + acute_vax * clear_vax + acute_vax_i * clear_vax_i
            )
            denominator += acute + acute_vax + acute_vax_i
        with np.errstate(divide="ignore", invalid="ignore"):
            weighted_clearance = np.where(denominator > 0, numerator / denominator, 0)
        extra = Series(
            tvec=tvec,
            vals=weighted_clearance,
            pop="Total",
            output="spontaneous_clearance",
            units="probability",
        )
        outputs_all.series.append(extra)

        return outputs_all

    # Function to save outputs
    def write_outputs(result, i, savedir):
        if isinstance(result, list):
            assert len(result) == 1, "result is a list with more than 1 item"
            result = result[0]
        outputs_all = mapping_function(result)
        name = result.name
        extracted = dict()
        for serie in outputs_all.series:
            tvec = serie.tvec
            pop = serie.pop
            par = serie.output
            vals = serie.vals
            extracted[(par, pop)] = vals
        extracted[("_name")] = name
        extracted[("_tvec")] = tvec

        filename = sc.makefilepath(
            filename=f"{name}_{i}_extracted.pkl", folder=savedir, makedirs=True
        )
        sc.save(filename=filename, obj=extracted)
        print(f"Saved outputs: {filename}")

    # Function to generate progbook
    parset = P.make_parset()
    parset.load_calibration(cal_folder / f"{country}_calibration.xlsx")

    # Run
    pset = P.load_progbook(
        results_folder / "progbooks" / f"progbook_{country}_central.xlsx"
    )
    progset_instructions, scenarios = define_scenarios(P, pset)

    sens_analyses = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="epi"
    )

    # Save central results to separate central folder as well
    for _, row in sens_analyses.iterrows():
        desc = row["description"]
        savedir = sens_folder / desc / country
        savedir.mkdir(parents=True, exist_ok=True)
        par_name = row["par_name"].split(",")
        value = row["value"]
        if isinstance(value, int):
            value = [value] * len(par_name)
        elif isinstance(value, str):
            value = value.split(",")
        for p_i, scen in zip(progset_instructions, scenarios):
            P_scen = P.copy()
            for par, val in zip(par_name, value):
                for pop in pops:
                    P_scen.data.tdve[par].ts[pop].assumption = float(val)
            parset_scen = P_scen.make_parset()
            parset_scen.load_calibration(cal_folder / f"{country}_calibration.xlsx")
            result_central = P_scen.run_sim(
                parset=parset_scen,
                progset=pset,
                progset_instructions=p_i,
                result_name=scen,
            )
            write_outputs(result_central, "central", savedir)

    sens_analyses = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="econ"
    )
    for _, row in sens_analyses.iterrows():
        desc = row["description"]
        savedir = sens_folder / desc / country
        savedir.mkdir(parents=True, exist_ok=True)
        savedir_og = results_folder / "scenarios" / country
        for scen in scenarios:
            filename = f"{scen}_central_extracted.pkl"
            shutil.copy(savedir_og / filename, savedir / filename)


def run_genotype_analyses(country, cal_folder, sens_folder, results_folder):
    """Run genotype analyses for a specified country using calibration and sensitivity data.
    
    Args:
        country (str): The ISO3 code of the country for which the analyses are to be run.
        cal_folder (str): The folder path containing calibration files.
        sens_folder (str): The folder path where sensitivity analysis results will be saved.
        results_folder (str): The folder path containing results files.
    
    Returns:
        None: The function performs analyses and saves the results to the specified folder.
        
    Raises:
        AssertionError: If the result is a list with more than one item during output extraction.
        
    Notes:
        This function loads calibration data, runs simulations based on various scenarios, 
        and writes the outputs to specified directories. It utilizes internal mapping functions 
        to process and aggregate results for different population groups and parameters.
    """
    P = project(
        country, load_calibration=True, cal_folder=cal_folder, load_programs=False
    )
    pops = P.data.pops.keys()  # define population

    # Function to define outputs to save
    def mapping_function(x):
        # Outputs to be saved
        outputs_by_pop = [
            "alive",
            "prevalence",
            "total_hcv",
            "inci_all_m",
            "notifications_m",
            "inci_m",
            "ab_tests_m",
            "pcr_tests_m",
        ]
        outputs_tot = [
            "total_hcv",
            "deaths_hcv_total",
            "tx_m",
            "inci_m",
            "notifications_m",
            "ab_tests_m",
            "pcr_tests_m",
            "f4_all",
            "d_cirrhosis_all",
            "inci_all_m",
            "hcc_inci",
            "acute_all",
            "undiag_all",
            "abpos_chronic",
            "pcr",
            "treated",
            "dc_udx",
            "hcc_udx",
            "f0f2_dx",
            "f3_dx",
            "f4_dx",
            "dc_dx",
            "hcc_dx",
            "d_cirrhosis_all",
            "cancer_all",
        ]
        outputs_tot_weighted = [
            "prop_f0f2",
            "prop_f3",
            "prop_f4",
            "prop_dc",
            "prop_hcc",
            "prevalence",
            "inci_100py",
        ]

        outputs_all = at.PlotData(
            x, outputs=outputs_by_pop, pops=pops, t_bins=1
        )  # outputs for each population group
        extra = at.PlotData(
            x, outputs=outputs_tot, pops="total", pop_aggregation="sum", t_bins=1
        )  # total outputs / summed
        outputs_all.series.extend(extra.series)
        extra = at.PlotData(
            x,
            outputs=outputs_tot_weighted,
            pops="total",
            pop_aggregation="weighted",
            t_bins=1,
        )  # total outputs / weighed
        outputs_all.series.extend(extra.series)

        # Aging pars
        transfer_pars = [
            "age_10-17_males_to_18-64_males",
            "age_10-17_females_to_18-64_females",
        ]
        transfer_pops = ["10-17_males", "10-17_females"]
        for par, pop in zip(transfer_pars, transfer_pops):
            extra = at.PlotData(x, outputs=par, pops=pop, t_bins=1)
            outputs_all.series.extend(extra.series)

        # Specific pop aggregations for econ fcts
        pars_tot = [
            "inci_all_m",
            "ab_all_m",
            "pcr_all_m",
            "f0f3_utx",
            "f4_utx",
            "f0f3_c",
            "f4_c",
            "inci_100py",
        ]
        pops_tot = [
            {
                "gen_pop": [
                    "0-9_males",
                    "0-9_females",
                    "10-17_males",
                    "10-17_females",
                    "18-64_males",
                    "18-64_females",
                    "65+_males",
                    "65+_females",
                ]
            },
            {"pwid": ["PWID_males", "PWID_females"]},
            {"prisoners": ["Prisoners_males", "Prisoners_females"]},
            {
                "working_age": [
                    "18-64_males",
                    "18-64_females",
                    "PWID_males",
                    "PWID_females",
                ]
            },
        ]
        for par in pars_tot:
            for pop in pops_tot:
                extra = at.PlotData(
                    x, outputs=par, pops=pop, pop_aggregation="sum", t_bins=1
                )
                outputs_all.series.extend(extra.series)

        pars_tot = ["rna_test_efficiency", "ab_test_efficiency", "prevalence"]
        extra = at.PlotData(
            x, outputs=pars_tot, pops=pops_tot, pop_aggregation="weighted", t_bins=1
        )
        outputs_all.series.extend(extra.series)

        pops_tot = [
            {"children": ["0-9_males", "0-9_females"]},
            {"teens": ["10-17_males", "10-17_females"]},
            {"adults": ["18-64_males", "18-64_females", "65+_males", "65+_females"]},
            {"pwid": ["PWID_males", "PWID_females"]},
            {"prisoners": ["Prisoners_males", "Prisoners_females"]},
        ]
        par = "vaccine_uptake:flow"
        extra = at.PlotData(
            x, outputs=par, pops=pops_tot, pop_aggregation="sum", t_bins=1
        )
        outputs_all.series.extend(extra.series)

        # Weighted spontaneous clearance rate
        tvec = (
            at.PlotData(x, outputs="acute", pops=pops, pop_aggregation="sum", t_bins=1)
            .series[0]
            .tvec
        )
        numerator = np.zeros_like(tvec)
        denominator = np.zeros_like(tvec)
        for pop in pops:
            acute = at.PlotData(x, outputs="acute", pops=pop, t_bins=1).series[0].vals
            acute_vax = (
                at.PlotData(x, outputs="acute_vax_reinfected", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            clear = (
                at.PlotData(x, outputs="prop_clear", pops=pop, t_bins=1).series[0].vals
            )
            clear_vax = (
                at.PlotData(x, outputs="vaccine_prop_clear", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            acute_vax_i = (
                at.PlotData(x, outputs="acute_vax_initial", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            clear_vax_i = (
                at.PlotData(x, outputs="vaccine_prop_clear_initial", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            numerator += (
                acute * clear + acute_vax * clear_vax + acute_vax_i * clear_vax_i
            )
            denominator += acute + acute_vax + acute_vax_i
        with np.errstate(divide="ignore", invalid="ignore"):
            weighted_clearance = np.where(denominator > 0, numerator / denominator, 0)
        extra = Series(
            tvec=tvec,
            vals=weighted_clearance,
            pop="Total",
            output="spontaneous_clearance",
            units="probability",
        )
        outputs_all.series.append(extra)

        return outputs_all

    # Function to save outputs
    def write_outputs(result, i, savedir):
        if isinstance(result, list):
            assert len(result) == 1, "result is a list with more than 1 item"
            result = result[0]
        outputs_all = mapping_function(result)
        name = result.name
        extracted = dict()
        for serie in outputs_all.series:
            tvec = serie.tvec
            pop = serie.pop
            par = serie.output
            vals = serie.vals
            extracted[(par, pop)] = vals
        extracted[("_name")] = name
        extracted[("_tvec")] = tvec

        filename = sc.makefilepath(
            filename=f"{name}_{i}_extracted.pkl", folder=savedir, makedirs=True
        )
        sc.save(filename=filename, obj=extracted)
        print(f"Saved outputs: {filename}")

    # Function to generate progbook
    parset = P.make_parset()
    parset.load_calibration(cal_folder / f"{country}_calibration.xlsx")

    # Run
    pset = P.load_progbook(
        results_folder / "progbooks" / f"progbook_{country}_central.xlsx"
    )
    progset_instructions, scenarios = define_scenarios(P, pset)

    sens_analyses = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="other"
    )
    sens_analyses = sens_analyses[sens_analyses.code == "genotype"]

    gen_analyses = pd.read_excel(
        str(rootdir) + "/data/flat_datasheet.xlsx",
        sheet_name="Vaccine-Genotype Country",
    )
    gen_analyses = gen_analyses[gen_analyses.ISO3 == country]
    # Assume 100% efficacy in one genotype at a time, 0% in others
    for _, row in sens_analyses.iterrows():
        desc = row["description"]
        savedir = sens_folder / desc / country
        savedir.mkdir(parents=True, exist_ok=True)
        par_name = row["par_name"].split(",")
        if desc == "100% efficacy against G1-G3-G4":
            value = np.sum([gen_analyses[val].values[0] for val in row['value'].split(',')])
        else:
            value = gen_analyses[row['value']].values[0]
        for p_i, scen in zip(progset_instructions, scenarios):
            P_scen = P.copy()
            for pop in pops:
                initial_value = P_scen.data.tdve["prop_clear"].ts[pop].assumption
                if value > initial_value:
                    for par in par_name:
                        P_scen.data.tdve[par].ts[pop].assumption = value
                else:
                    for par in par_name:
                        P_scen.data.tdve[par].ts[pop].assumption = initial_value
            parset_scen = P_scen.make_parset()
            parset_scen.load_calibration(cal_folder / f"{country}_calibration.xlsx")
            result_central = P_scen.run_sim(
                parset=parset_scen,
                progset=pset,
                progset_instructions=p_i,
                result_name=scen,
            )
            write_outputs(result_central, "central", savedir)


def run_coverage_analyses(country, cal_folder, sens_folder):
    """Run coverage analyses for a specified country using calibration and sensitivity folders.
    
    Args:
        country (str): The name of the country for which the analyses are to be run.
        cal_folder (str): The folder path containing calibration data.
        sens_folder (str): The folder path for saving sensitivity analysis results.
    
    Returns:
        None: The function performs analyses and saves the results to specified directories.
    
    Raises:
        AssertionError: If the result is a list with more than one item during output extraction.
        
    Notes:
        This function loads calibration data, runs simulations for various scenarios, and generates output files
        containing the results of the coverage analyses. It also generates program books for the specified scenarios.
    """
    P = project(
        country, load_calibration=True, cal_folder=cal_folder, load_programs=False
    )
    pops = P.data.pops.keys()  # define population

    # Function to define outputs to save
    def mapping_function(x):
        # Outputs to be saved
        outputs_by_pop = [
            "alive",
            "prevalence",
            "total_hcv",
            "inci_all_m",
            "notifications_m",
            "inci_m",
            "ab_tests_m",
            "pcr_tests_m",
        ]
        outputs_tot = [
            "total_hcv",
            "deaths_hcv_total",
            "tx_m",
            "inci_m",
            "notifications_m",
            "ab_tests_m",
            "pcr_tests_m",
            "f4_all",
            "d_cirrhosis_all",
            "inci_all_m",
            "hcc_inci",
            "acute_all",
            "undiag_all",
            "abpos_chronic",
            "pcr",
            "treated",
            "dc_udx",
            "hcc_udx",
            "f0f2_dx",
            "f3_dx",
            "f4_dx",
            "dc_dx",
            "hcc_dx",
            "d_cirrhosis_all",
            "cancer_all",
        ]
        outputs_tot_weighted = [
            "prop_f0f2",
            "prop_f3",
            "prop_f4",
            "prop_dc",
            "prop_hcc",
            "prevalence",
            "inci_100py",
        ]

        outputs_all = at.PlotData(
            x, outputs=outputs_by_pop, pops=pops, t_bins=1
        )  # outputs for each population group
        extra = at.PlotData(
            x, outputs=outputs_tot, pops="total", pop_aggregation="sum", t_bins=1
        )  # total outputs / summed
        outputs_all.series.extend(extra.series)
        extra = at.PlotData(
            x,
            outputs=outputs_tot_weighted,
            pops="total",
            pop_aggregation="weighted",
            t_bins=1,
        )  # total outputs / weighed
        outputs_all.series.extend(extra.series)

        # Aging pars
        transfer_pars = [
            "age_10-17_males_to_18-64_males",
            "age_10-17_females_to_18-64_females",
        ]
        transfer_pops = ["10-17_males", "10-17_females"]
        for par, pop in zip(transfer_pars, transfer_pops):
            extra = at.PlotData(x, outputs=par, pops=pop, t_bins=1)
            outputs_all.series.extend(extra.series)

        # Specific pop aggregations for econ fcts
        pars_tot = [
            "inci_all_m",
            "ab_all_m",
            "pcr_all_m",
            "f0f3_utx",
            "f4_utx",
            "f0f3_c",
            "f4_c",
            "inci_100py",
        ]
        pops_tot = [
            {
                "gen_pop": [
                    "0-9_males",
                    "0-9_females",
                    "10-17_males",
                    "10-17_females",
                    "18-64_males",
                    "18-64_females",
                    "65+_males",
                    "65+_females",
                ]
            },
            {"pwid": ["PWID_males", "PWID_females"]},
            {"prisoners": ["Prisoners_males", "Prisoners_females"]},
            {
                "working_age": [
                    "18-64_males",
                    "18-64_females",
                    "PWID_males",
                    "PWID_females",
                ]
            },
        ]
        for par in pars_tot:
            for pop in pops_tot:
                extra = at.PlotData(
                    x, outputs=par, pops=pop, pop_aggregation="sum", t_bins=1
                )
                outputs_all.series.extend(extra.series)

        pars_tot = ["rna_test_efficiency", "ab_test_efficiency", "prevalence"]
        extra = at.PlotData(
            x, outputs=pars_tot, pops=pops_tot, pop_aggregation="weighted", t_bins=1
        )
        outputs_all.series.extend(extra.series)

        pops_tot = [
            {"children": ["0-9_males", "0-9_females"]},
            {"teens": ["10-17_males", "10-17_females"]},
            {"adults": ["18-64_males", "18-64_females", "65+_males", "65+_females"]},
            {"pwid": ["PWID_males", "PWID_females"]},
            {"prisoners": ["Prisoners_males", "Prisoners_females"]},
        ]
        par = "vaccine_uptake:flow"
        extra = at.PlotData(
            x, outputs=par, pops=pops_tot, pop_aggregation="sum", t_bins=1
        )
        outputs_all.series.extend(extra.series)

        # Weighted spontaneous clearance rate
        tvec = (
            at.PlotData(x, outputs="acute", pops=pops, pop_aggregation="sum", t_bins=1)
            .series[0]
            .tvec
        )
        numerator = np.zeros_like(tvec)
        denominator = np.zeros_like(tvec)
        for pop in pops:
            acute = at.PlotData(x, outputs="acute", pops=pop, t_bins=1).series[0].vals
            acute_vax = (
                at.PlotData(x, outputs="acute_vax_reinfected", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            clear = (
                at.PlotData(x, outputs="prop_clear", pops=pop, t_bins=1).series[0].vals
            )
            clear_vax = (
                at.PlotData(x, outputs="vaccine_prop_clear", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            acute_vax_i = (
                at.PlotData(x, outputs="acute_vax_initial", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            clear_vax_i = (
                at.PlotData(x, outputs="vaccine_prop_clear_initial", pops=pop, t_bins=1)
                .series[0]
                .vals
            )
            numerator += (
                acute * clear + acute_vax * clear_vax + acute_vax_i * clear_vax_i
            )
            denominator += acute + acute_vax + acute_vax_i
        with np.errstate(divide="ignore", invalid="ignore"):
            weighted_clearance = np.where(denominator > 0, numerator / denominator, 0)
        extra = Series(
            tvec=tvec,
            vals=weighted_clearance,
            pop="Total",
            output="spontaneous_clearance",
            units="probability",
        )
        outputs_all.series.append(extra)

        return outputs_all

    # Function to save outputs
    def write_outputs(result, i):
        if isinstance(result, list):
            assert len(result) == 1, "result is a list with more than 1 item"
            result = result[0]
        outputs_all = mapping_function(result)
        name = result.name
        extracted = dict()
        for serie in outputs_all.series:
            tvec = serie.tvec
            pop = serie.pop
            par = serie.output
            vals = serie.vals
            extracted[(par, pop)] = vals
        extracted[("_name")] = name
        extracted[("_tvec")] = tvec

        filename = sc.makefilepath(
            filename=f"{name}_{i}_extracted.pkl", folder=savedir, makedirs=True
        )
        sc.save(filename=filename, obj=extracted)
        print(f"Saved outputs: {filename}")

    # Function to generate progbook
    def gen_pb(result, i, cov):
        from hcv.generate_progbooks import generate_progbook
        if isinstance(result, list):
            assert len(result) == 1, "result is a list with more than 1 item"
            result = result[0]
        generate_progbook(
            country,
            result=result,
            savedir=savedir_pb / f"progbook_{country}_{i}.xlsx",
            cov_scenario=True,
            cov=cov,
        )

    sens_analyses = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="other"
    )
    sens_analyses = sens_analyses[sens_analyses.code == "coverage"]

    for _, row in sens_analyses.iterrows():
        desc = row["description"]
        value = row["value"]
        # save results
        savedir = sens_folder / desc / country
        savedir.mkdir(parents=True, exist_ok=True)
        # save progbook
        savedir_pb = sens_folder / desc / "progbooks"
        savedir_pb.mkdir(parents=True, exist_ok=True)
        parset = P.make_parset()
        parset.load_calibration(cal_folder / f"{country}_calibration.xlsx")
        # Run with no sampling
        result_central = P.run_sim(parset=parset)
        gen_pb(result_central, "central", int(value * 100))
        pset = P.load_progbook(savedir_pb / f"progbook_{country}_central.xlsx")
        progset_instructions, scenarios = define_scenarios(P, pset)
        for p_i, scen in zip(progset_instructions, scenarios):
            result_central = P.run_sim(
                parset=parset, progset=pset, progset_instructions=p_i, result_name=scen
            )
            write_outputs(result_central, "central")


def econ_eval_central(country, sens_folder):
    """Evaluate the economic impact of HCV scenarios for a given country.
    
    This function processes various economic and epidemiological data related to Hepatitis C Virus (HCV) for a specified country and scenario. It aggregates data from multiple sources, performs calculations related to incidence, prevalence, treatment costs, and productivity losses, and outputs the results in a structured format.
    
    Args:
        country (str): The ISO3 code of the country for which the economic evaluation is performed.
        sens_folder (Path): The path to the folder containing sensitivity analysis data.
    
    Returns:
        None: The function saves the aggregated results to a specified file in the provided sensitivity folder.
    
    Raises:
        FileNotFoundError: If the required data files are not found in the specified paths.
        ValueError: If the data format is not as expected or if there are inconsistencies in the data.
    
    Example:
        econ_eval_central("USA", Path("/path/to/sensitivity/folder"))
    """
    # Import scenario names and definitions for loop and mapping
    scenarios = pd.read_excel(
        str(rootdir) + "/data/progbook_inputs.xlsx", sheet_name="scenarios"
    )
    scen_name = list(pd.unique(scenarios.scenario_name))

    sim_loop = ["central"]

    # Aggregate Summary Data (deaths will be used for counterfactual YLL/YPLLs)
    agg_data = {}
    for scen in scen_name:
        agg_data[scen] = {}

    colnames = ["year", "central"]
    agg_outs = [
        "HCV incidence",
        "HCC incidence",
        "HCV incidence per 100py",
        "HCV mortality",
        "HCV prevalence",
        "PWID prevalence",
        "Prisoner prevalence",
    ]
    for scen in scen_name:
        for out in agg_outs:
            agg_data[scen][out] = pd.DataFrame(columns=colnames)
            agg_data[scen][out]["year"] = np.arange(2000.5, 2051.5, 1)

    sens_analyses_sheet_1 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="epi"
    )
    sens_analyses_sheet_2 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="econ"
    )
    sens_analyses_sheet_3 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="other"
    )
    sens_analyses = (
        list(sens_analyses_sheet_1.description.values)
        + list(sens_analyses_sheet_2.description.values)
        + list(sens_analyses_sheet_3.description.values)
    )
    for analysis in sens_analyses:
        spont_clear = {}
        for scen in scen_name:
            spont_clear[scen] = pd.DataFrame(columns=["year", "spontaneous_clearance"])
            spont_clear[scen].year = np.arange(2000.5, 2051.5, 1)
            for i, sim in enumerate(sim_loop):
                data = sc.load(
                    sens_folder / analysis / country / f"{scen}_{sim}_extracted.pkl"
                )
                # Total Population
                agg_data[scen]["HCV incidence"].iloc[:, i + 1] = data[
                    "inci_all_m", "Total"
                ]
                agg_data[scen]["HCC incidence"].iloc[:, i + 1] = data[
                    "hcc_inci", "Total"
                ]
                agg_data[scen]["HCV incidence per 100py"].iloc[:, i + 1] = data[
                    "inci_100py", "Total"
                ]
                agg_data[scen]["HCV mortality"].iloc[:, i + 1] = data[
                    "deaths_hcv_total", "Total"
                ]
                agg_data[scen]["HCV prevalence"].iloc[:, i + 1] = data[
                    "total_hcv", "Total"
                ]
                # Key Populations
                agg_data[scen]["PWID prevalence"].iloc[:, i + 1] = (
                    data["total_hcv", "Prisoners_males"]
                    + data["total_hcv", "Prisoners_females"]
                )
                agg_data[scen]["Prisoner prevalence"].iloc[:, i + 1] = (
                    data["total_hcv", "PWID_males"] + data["total_hcv", "PWID_females"]
                )
                spont_clear[scen]["spontaneous_clearance"] = data[
                    "spontaneous_clearance", "Total"
                ]

        # Diagnostic Testing, Positive Tests, Treatment and Vaccination Coverage (most can be used for costs)
        util_data = {}
        for scen in scen_name:
            util_data[scen] = {}

        util_outs = [
            "PWID Ab Tests",
            "Prisoner Ab Tests",
            "Gen Ab Tests",
            "Total Ab Tests",
            "PWID Ab Positive",
            "Prisoner Ab Positive",
            "Gen Ab Positive",
            "Total Ab Positive",
            "PWID RNA Tests",
            "Prisoner RNA Tests",
            "Gen RNA Tests",
            "Total RNA Tests",
            "PWID RNA Positive",
            "Prisoner RNA Positive",
            "Gen RNA Positive",
            "Total RNA Positive",
            "Total treatment",
            "PWID vaccinations",
            "Prisoner vaccinations",
            "Gen vaccinations",
            "Total vaccinations",
        ]

        for scen in scen_name:
            for out in util_outs:
                util_data[scen][out] = pd.DataFrame(columns=colnames)
                util_data[scen][out]["year"] = np.arange(2000.5, 2051.5, 1)

        for scen in scen_name:
            for i, sim in enumerate(sim_loop):
                data = sc.load(
                    sens_folder / analysis / country / f"{scen}_{sim}_extracted.pkl"
                )

                # Antibody Tests (total)
                util_data[scen]["PWID Ab Tests"].iloc[:, i + 1] = (
                    data["ab_tests_m", "PWID_males"]
                    + data["ab_tests_m", "PWID_females"]
                )
                util_data[scen]["Prisoner Ab Tests"].iloc[:, i + 1] = (
                    data["ab_tests_m", "Prisoners_males"]
                    + data["ab_tests_m", "Prisoners_females"]
                )
                util_data[scen]["Total Ab Tests"].iloc[:, i + 1] = data[
                    "ab_tests_m", "Total"
                ]
                util_data[scen]["Gen Ab Tests"].iloc[:, i + 1] = util_data[scen][
                    "Total Ab Tests"
                ].iloc[:, i + 1] - (
                    util_data[scen]["PWID Ab Tests"].iloc[:, i + 1]
                    + util_data[scen]["Prisoner Ab Tests"].iloc[:, i + 1]
                )  # acts as a double check

                # Antibody Tests (positive)
                util_data[scen]["PWID Ab Positive"].iloc[:, i + 1] = data[
                    "ab_all_m", "pwid"
                ]
                util_data[scen]["Prisoner Ab Positive"].iloc[:, i + 1] = data[
                    "ab_all_m", "prisoners"
                ]
                util_data[scen]["Gen Ab Positive"].iloc[:, i + 1] = data[
                    "ab_all_m", "gen_pop"
                ]
                util_data[scen]["Total Ab Positive"].iloc[:, i + 1] = (
                    util_data[scen]["PWID Ab Positive"].iloc[:, i + 1]
                    + util_data[scen]["Prisoner Ab Positive"].iloc[:, i + 1]
                    + util_data[scen]["Gen Ab Positive"].iloc[:, i + 1]
                )  # acts as a double check

                # RNA tests (total)
                util_data[scen]["PWID RNA Tests"].iloc[:, i + 1] = data[
                    "ab_all_m", "pwid"
                ] * (1 + spont_clear[scen].iloc[:, 1])
                util_data[scen]["Prisoner RNA Tests"].iloc[:, i + 1] = data[
                    "ab_all_m", "prisoners"
                ] * (1 + spont_clear[scen].iloc[:, 1])
                util_data[scen]["Gen RNA Tests"].iloc[:, i + 1] = data[
                    "ab_all_m", "gen_pop"
                ] * (1 + spont_clear[scen].iloc[:, 1])
                util_data[scen]["Total RNA Tests"].iloc[:, i + 1] = (
                    util_data[scen]["PWID RNA Tests"].iloc[:, i + 1]
                    + util_data[scen]["Prisoner RNA Tests"].iloc[:, i + 1]
                    + util_data[scen]["Gen RNA Tests"].iloc[:, i + 1]
                )  # acts as a double check

                # RNA tests (positive)
                util_data[scen]["PWID RNA Positive"].iloc[:, i + 1] = data[
                    "pcr_all_m", "pwid"
                ]
                util_data[scen]["Prisoner RNA Positive"].iloc[:, i + 1] = data[
                    "pcr_all_m", "prisoners"
                ]
                util_data[scen]["Gen RNA Positive"].iloc[:, i + 1] = data[
                    "pcr_all_m", "gen_pop"
                ]
                util_data[scen]["Total RNA Positive"].iloc[:, i + 1] = (
                    util_data[scen]["PWID RNA Positive"].iloc[:, i + 1]
                    + util_data[scen]["Prisoner RNA Positive"].iloc[:, i + 1]
                    + util_data[scen]["Gen RNA Positive"].iloc[:, i + 1]
                )  # acts as a a double check

                # Treatment initiations
                util_data[scen]["Total treatment"].iloc[:, i + 1] = data[
                    "tx_m", "Total"
                ]

                # Vaccinations administered
                util_data[scen]["PWID vaccinations"].iloc[:, i + 1] = data[
                    "vaccine_uptake:flow", "pwid"
                ]
                util_data[scen]["Prisoner vaccinations"].iloc[:, i + 1] = data[
                    "vaccine_uptake:flow", "prisoners"
                ]
                util_data[scen]["Gen vaccinations"].iloc[:, i + 1] = (
                    data["vaccine_uptake:flow", "children"]
                    + data["vaccine_uptake:flow", "teens"] * 1 / 8
                    + data["vaccine_uptake:flow", "adults"]
                )
                util_data[scen]["Total vaccinations"].iloc[:, i + 1] = (
                    util_data[scen]["PWID vaccinations"].iloc[:, i + 1]
                    + util_data[scen]["Prisoner vaccinations"].iloc[:, i + 1]
                    + util_data[scen]["Gen vaccinations"].iloc[:, i + 1]
                )

        # Disease Stages (can be used for productivity, DALYs and management costs)
        dis_data = {}
        for scen in scen_name:
            dis_data[scen] = {}

        dis_outs = [
            "F0-F3 infected",
            "F4+ infected",
            "F0-F3 cured",
            "F4+ cured",
            "F0-F2 diagnosed",
            "F3 diagnosed",
            "F4 diagnosed",
            "DC diagnosed",
            "HCC diagnosed",
            "DC undiagnosed",
            "HCC undiagnosed",
            "Total F4",
            "Total Decompensated",
            "Total Liver Cancer",
        ]

        for scen in scen_name:
            for out in dis_outs:
                dis_data[scen][out] = pd.DataFrame(columns=colnames)
                dis_data[scen][out]["year"] = np.arange(2000.5, 2051.5, 1)

        for scen in scen_name:
            for i, sim in enumerate(sim_loop):
                data = sc.load(
                    sens_folder / analysis / country / f"{scen}_{sim}_extracted.pkl"
                )
                # Productivity
                dis_data[scen]["F0-F3 infected"].iloc[:, i + 1] = data[
                    "f0f3_utx", "working_age"
                ]
                dis_data[scen]["F4+ infected"].iloc[:, i + 1] = data[
                    "f4_utx", "working_age"
                ]
                dis_data[scen]["F0-F3 cured"].iloc[:, i + 1] = data[
                    "f0f3_c", "working_age"
                ]
                dis_data[scen]["F4+ cured"].iloc[:, i + 1] = data["f4_c", "working_age"]

                # Disease Management Costs
                dis_data[scen]["F0-F2 diagnosed"].iloc[:, i + 1] = data[
                    "f0f2_dx", "Total"
                ]
                dis_data[scen]["F3 diagnosed"].iloc[:, i + 1] = data["f3_dx", "Total"]
                dis_data[scen]["F4 diagnosed"].iloc[:, i + 1] = data["f4_dx", "Total"]
                dis_data[scen]["DC diagnosed"].iloc[:, i + 1] = data["dc_dx", "Total"]
                dis_data[scen]["HCC diagnosed"].iloc[:, i + 1] = data["hcc_dx", "Total"]
                dis_data[scen]["DC undiagnosed"].iloc[:, i + 1] = data[
                    "dc_udx", "Total"
                ]
                dis_data[scen]["HCC undiagnosed"].iloc[:, i + 1] = data[
                    "hcc_udx", "Total"
                ]

                # DALYs
                dis_data[scen]["Total F4"].iloc[:, i + 1] = data["f4_all", "Total"]
                dis_data[scen]["Total Decompensated"].iloc[:, i + 1] = data[
                    "d_cirrhosis_all", "Total"
                ]
                dis_data[scen]["Total Liver Cancer"].iloc[:, i + 1] = data[
                    "cancer_all", "Total"
                ]

        # Call data for economic analysis
        global_pars = pd.read_excel(
            str(rootdir) + "/data/flat_datasheet.xlsx",
            sheet_name="Cost - Global Inputs",
        )
        dir_costs = pd.read_excel(
            str(rootdir) + "/data/flat_datasheet.xlsx",
            sheet_name="Cost - Direct Costs",
        )
        dir_costs = dir_costs[dir_costs["ISO3"] == country]
        prod_costs = pd.read_excel(
            str(rootdir) + "/data/flat_datasheet.xlsx",
            sheet_name="Cost - YLL and productivity",
        )
        prod_costs = prod_costs[prod_costs["ISO3"] == country]

        # Not Sampled Economic Pars
        start_year = global_pars.iloc[10, 1]
        c_disc, h_disc = global_pars.iloc[9, 1], global_pars.iloc[8, 1]
        vax_weight = dir_costs.vax_weight.values[0]
        abs_f, abs_fr, pres_f, pres_fr, abs_c, abs_cr, pres_c, pres_cr = global_pars[
            "value"
        ].iloc[0:8]
        waste, logistics, overheads, prog_costs = global_pars["value"].iloc[17:21]
        employed, gdp, hr_salary = (
            prod_costs.ER_total.values[0],
            prod_costs.GDP.values[0],
            prod_costs.hour_wage.values[0],
        )

        # Sampled Economic Pars (central value, n_sims sampled in array)
        f0f2_cost = [dir_costs.f0f2.values[0]]
        f3_cost = [dir_costs.f3.values[0]]
        f4_cost = [dir_costs.f4.values[0]]
        dc_cost = [dir_costs.DC.values[0]]
        hcc_cost = [dir_costs.HCC.values[0]]
        ab_test = [global_pars["value"].iloc[14]]
        rna_test = [global_pars["value"].iloc[15]]
        f4_weight = [global_pars["value"].iloc[21]]
        dc_weight = [global_pars["value"].iloc[11]]
        hcc_weight = [global_pars["value"].iloc[12]]

        # Discount Arrays
        cost_disc, health_disc = np.zeros((len(np.arange(2000, 2051, 1)), 2)), np.zeros(
            (len(np.arange(2000, 2051, 1)), 2)
        )
        cost_disc[:, 0], health_disc[:, 0] = np.arange(2000.5, 2051.5, 1), np.arange(
            2000.5, 2051.5, 1
        )
        for i in range(len(cost_disc[:, 0])):
            if cost_disc[i, 0] < start_year - 0.5:
                cost_disc[i, 1], health_disc[i, 1] = 0, 1
            else:
                cost_disc[i, 1] = (1 - c_disc) ** (cost_disc[i, 0] - (start_year - 0.5))
                health_disc[i, 1] = (1 - h_disc) ** (
                    cost_disc[i, 0] - (start_year - 0.5)
                )

        # Calculate YLL and YPLL (counterfactual method)
        aging = [
            1 / 15,
            1 / 15,
            1 / 20,
            1 / 15,
        ]  # Age bracket length (65+ retire and no longer contribute)
        death_dist = [
            prod_costs["014_prp"].values[0],
            prod_costs["1529_prp"].values[0],
            prod_costs["3049_prp"].values[0],
            prod_costs["5059_prp"].values[0] + 0.5 * prod_costs["6069_prp"].values[0],
            0.5 * prod_costs["6069_prp"].values[0] + prod_costs["70_prp"].values[0],
        ]
        acm = [
            prod_costs["mrate_014"].values[0],
            prod_costs["mrate_1529"].values[0],
            prod_costs["mrate_3049"].values[0],
            prod_costs["mrate_5064"].values[0],
            1 / prod_costs["le_70"].values[0],
        ]

        ghosts = {}
        for scen in scen_name:
            ghosts[scen] = {}
            for i, sim in enumerate(sim_loop):
                ghosts[scen][f"sim_{sim}"] = np.zeros(
                    (len(agg_data[scen]["HCV mortality"]), len(acm))
                )

        for scen in scen_name:
            for i, sim in enumerate(sim_loop):
                for t in range(1, len(agg_data[scen]["HCV mortality"])):
                    ppl_who_age = (
                        ghosts[scen][f"sim_{sim}"][t - 1, 0 : len(aging)] * aging
                    )
                    ghosts[scen][f"sim_{sim}"][t, 0] = max(
                        0,
                        ghosts[scen][f"sim_{sim}"][t - 1, 0]
                        - ppl_who_age[0]
                        - acm[0] * ghosts[scen][f"sim_{sim}"][t - 1, 0],
                    )  # 0-14
                    ghosts[scen][f"sim_{sim}"][t, 1] = max(
                        0,
                        ghosts[scen][f"sim_{sim}"][t - 1, 1]
                        - ppl_who_age[1]
                        + ppl_who_age[0]
                        - acm[1] * ghosts[scen][f"sim_{sim}"][t - 1, 1],
                    )  # 15-29
                    ghosts[scen][f"sim_{sim}"][t, 2] = max(
                        0,
                        ghosts[scen][f"sim_{sim}"][t - 1, 2]
                        - ppl_who_age[2]
                        + ppl_who_age[1]
                        - acm[2] * ghosts[scen][f"sim_{sim}"][t - 1, 2],
                    )  # 30-49
                    ghosts[scen][f"sim_{sim}"][t, 3] = max(
                        0,
                        ghosts[scen][f"sim_{sim}"][t - 1, 3]
                        - ppl_who_age[3]
                        + ppl_who_age[2]
                        - acm[3] * ghosts[scen][f"sim_{sim}"][t - 1, 3],
                    )  # 50-64
                    ghosts[scen][f"sim_{sim}"][t, 4] = max(
                        0,
                        ghosts[scen][f"sim_{sim}"][t - 1, 4]
                        + ppl_who_age[3]
                        - acm[4] * ghosts[scen][f"sim_{sim}"][t - 1, 4],
                    )  # 65+

                    for dage in range(len(death_dist)):
                        ghosts[scen][f"sim_{sim}"][t, dage] = (
                            ghosts[scen][f"sim_{sim}"][t, dage]
                            + agg_data[scen]["HCV mortality"].iloc[t, i + 1]
                            * death_dist[dage]
                        )

        # Collate outcomes from counterfactuals
        lives_lost = {}
        for scen in scen_name:
            lives_lost[scen] = {}

        lives_outs = ["Years Life Lost", "Years Productive Life Lost"]

        for scen in scen_name:
            for out in lives_outs:
                lives_lost[scen][out] = pd.DataFrame(columns=colnames)
                lives_lost[scen][out]["year"] = np.arange(2000.5, 2051.5, 1)

        for scen in scen_name:
            for i, sim in enumerate(sim_loop):
                lives_lost[scen]["Years Life Lost"].iloc[:, i + 1] = np.sum(
                    ghosts[scen][f"sim_{sim}"], axis=1
                )
                lives_lost[scen]["Years Productive Life Lost"].iloc[:, i + 1] = np.sum(
                    ghosts[scen][f"sim_{sim}"][:, 1:4], axis=1
                )

        # Calculate costs
        econ_ests = {}
        for scen in scen_name:
            econ_ests[scen] = {}

        econ_outs = [
            "Diagnosis Cost",
            "Disease Management Cost",
            "Treatment Cost",
            "Morbidity Productivity",
            "Mortality Productivity",
            "Total Productivity",
            "Total Cost",
            "Vaccine Cost",
            "YLD",
            "DALYs",
        ]

        for scen in scen_name:
            for out in econ_outs:
                econ_ests[scen][out] = pd.DataFrame(columns=colnames)
                econ_ests[scen][out]["year"] = np.arange(2000.5, 2051.5, 1)

        vax_cost = 5
        for scen in scen_name:
            if analysis == "DAA Treatment":
                daa = [dir_costs.daa_cost_ub.values[0]]
            else:
                daa = [global_pars["value"].iloc[16]]
            if analysis == "Engagement in care lb":
                f0f3_care, f4dx_care, dchcc_care = (
                    dir_costs.f0f3_dx_care_lb.values[0],
                    dir_costs.f4_dx_care_lb.values[0],
                    dir_costs.dchcc_care_lb.values[0],
                )
            elif analysis == "Engagement in care ub":
                f0f3_care, f4dx_care, dchcc_care = (
                    dir_costs.f0f3_dx_care_ub.values[0],
                    dir_costs.f4_dx_care_ub.values[0],
                    dir_costs.dchcc_care_ub.values[0],
                )
            else:
                f0f3_care, f4dx_care, dchcc_care = (
                    dir_costs.f0f3_dx_care.values[0],
                    dir_costs.f4_dx_care.values[0],
                    dir_costs.dchcc_care.values[0],
                )
            for c in [10, 15, 20, 25, 30]:
                if analysis == f"Vaccine unit cost ${c}":
                    vax_cost = c
                    break
            for i, sim in enumerate(sim_loop):
                for t in range(len(econ_ests[scen]["Diagnosis Cost"]["year"])):
                    econ_ests[scen]["Diagnosis Cost"].iloc[t, i + 1] = (
                        (
                            util_data[scen]["PWID Ab Tests"].iloc[t, i + 1]
                            + util_data[scen]["Prisoner Ab Tests"].iloc[t, i + 1]
                            + util_data[scen]["Gen Ab Tests"].iloc[t, i + 1]
                        )
                        * (
                            (ab_test[i] * (1 + waste + logistics + prog_costs))
                            + (2 * (hr_salary * (1 + overheads + prog_costs)))
                        )
                        + (
                            util_data[scen]["PWID RNA Tests"].iloc[t, i + 1]
                            + util_data[scen]["Prisoner RNA Tests"].iloc[t, i + 1]
                            + util_data[scen]["Gen RNA Tests"].iloc[t, i + 1]
                        )
                        * (
                            (rna_test[i] * (1 + waste + logistics + prog_costs))
                            + (2 * (hr_salary * (1 + overheads + prog_costs)))
                        )
                    ) * cost_disc[t, 1]

                    econ_ests[scen]["Disease Management Cost"].iloc[t, i + 1] = (
                        (
                            dis_data[scen]["F0-F2 diagnosed"].iloc[t, i + 1]
                            * f0f2_cost[i]
                            * f0f3_care
                        )
                        + (
                            dis_data[scen]["F3 diagnosed"].iloc[t, i + 1]
                            * f3_cost[i]
                            * f0f3_care
                        )
                        + (
                            dis_data[scen]["F4 diagnosed"].iloc[t, i + 1]
                            * f4_cost[i]
                            * f4dx_care
                        )
                        + (
                            (
                                dis_data[scen]["DC diagnosed"].iloc[t, i + 1]
                                + dis_data[scen]["DC undiagnosed"].iloc[t, i + 1]
                            )
                            * dc_cost[i]
                            * dchcc_care
                        )
                        + (
                            (
                                dis_data[scen]["HCC undiagnosed"].iloc[t, i + 1]
                                + dis_data[scen]["HCC diagnosed"].iloc[t, i + 1]
                            )
                            * hcc_cost[i]
                            * dchcc_care
                        )
                    ) * cost_disc[t, 1]

                    econ_ests[scen]["Treatment Cost"].iloc[t, i + 1] = (
                        util_data[scen]["Total treatment"].iloc[t, i + 1]
                        * (
                            daa[i] * (1 + waste + logistics + prog_costs)
                            + hr_salary * (1 + overheads + prog_costs)
                        )
                        * cost_disc[t, 1]
                    )

                    econ_ests[scen]["Morbidity Productivity"].iloc[t, i + 1] = (
                        (gdp / employed)
                        * employed
                        * (
                            dis_data[scen]["F0-F3 infected"].iloc[t, i + 1]
                            * (abs_f + pres_f)
                            + dis_data[scen]["F0-F3 cured"].iloc[t, i + 1]
                            * (abs_fr + pres_fr)
                            + dis_data[scen]["F4+ infected"].iloc[t, i + 1]
                            * (abs_c + pres_c)
                            + dis_data[scen]["F4+ cured"].iloc[t, i + 1]
                            * (abs_cr + pres_cr)
                        )
                        * cost_disc[t, 1]
                    )

                    econ_ests[scen]["Mortality Productivity"].iloc[t, i + 1] = (
                        (gdp / employed)
                        * employed
                        * lives_lost[scen]["Years Productive Life Lost"].iloc[t, i + 1]
                        * cost_disc[t, 1]
                    )

                    econ_ests[scen]["Total Productivity"].iloc[t, i + 1] = (
                        econ_ests[scen]["Morbidity Productivity"].iloc[t, i + 1]
                        + econ_ests[scen]["Mortality Productivity"].iloc[t, i + 1]
                    )

                    econ_ests[scen]["Total Cost"].iloc[t, i + 1] = (
                        econ_ests[scen]["Diagnosis Cost"].iloc[t, i + 1]
                        + econ_ests[scen]["Disease Management Cost"].iloc[t, i + 1]
                        + econ_ests[scen]["Treatment Cost"].iloc[t, i + 1]
                        + econ_ests[scen]["Total Productivity"].iloc[t, i + 1]
                    )

                    econ_ests[scen]["Vaccine Cost"].iloc[t, i + 1] = (
                        (
                            util_data[scen]["PWID vaccinations"].iloc[t, i + 1]
                            + util_data[scen]["Prisoner vaccinations"].iloc[t, i + 1]
                            + util_data[scen]["Gen vaccinations"].iloc[t, i + 1]
                        )
                        * (
                            vax_cost * vax_weight * (1 + waste + logistics + prog_costs)
                            + (hr_salary * (1 / 12)) * (1 + overheads + prog_costs)
                        )
                        * cost_disc[t, 1]
                    )

                    econ_ests[scen]["YLD"].iloc[t, i + 1] = (
                        dis_data[scen]["Total Decompensated"].iloc[t, i + 1]
                        * dc_weight[i]
                        + dis_data[scen]["Total Liver Cancer"].iloc[t, i + 1]
                        * hcc_weight[i]
                        + dis_data[scen]["Total F4"].iloc[t, i + 1] * f4_weight[i]
                    )

                    econ_ests[scen]["DALYs"].iloc[t, i + 1] = (
                        econ_ests[scen]["YLD"].iloc[t, i + 1]
                        + lives_lost[scen]["Years Life Lost"].iloc[t, i + 1]
                    ) * health_disc[t, 1]

        # Aggregate Dictionary (of dictionaries of dictionaries....)
        output_dict = {}
        output_dict["agg_data"] = agg_data
        output_dict["util_data"] = util_data
        output_dict["dis_data"] = dis_data
        output_dict["lives_lost"] = lives_lost
        output_dict["econ_ests"] = econ_ests

        # Save as pkl file
        filename = sc.makefilepath(
            filename=f"{country}_econ_eval.pkl",
            folder=sens_folder / "agg_outputs" / analysis,
            makedirs=True,
        )
        sc.save(filename=filename, obj=output_dict)


def aggregate_ensembles_central(sens_folder, results_folder):
    """Aggregate epidemiological and economic data from multiple countries into regional and global summaries (central results only).
    
    This function reads sensitivity analysis data from specified folders, aggregates the results by region and scenario, and saves the aggregated data into specified output files. It processes both epidemiological and economic outputs, including various health metrics and costs associated with disease management.
    
    Args:
        sens_folder (str or Path): The path to the folder containing sensitivity analysis outputs.
        results_folder (str or Path): The path to the folder where aggregated results will be saved.
    
    Returns:
        None: The function saves the aggregated data to files but does not return any values.
    """

    ref_countries = pd.read_excel(
        str(rootdir) + "/data/flat_datasheet.xlsx",
        sheet_name="Cost - YLL and productivity",
    ).iloc[:, np.r_[1, 4]]
    who_reg = list(pd.unique(ref_countries.WHO_reg)) + ["global", "top10"]

    loop_folder = (
        results_folder / "scenarios" / "agg_outputs"
    )  # pathlib.Path(str(sens_folder)+"/agg_outputs/")
    countries = [
        country.stem[:3] for country in loop_folder.iterdir() if country.is_file()
    ]

    scenarios = pd.read_excel(
        str(rootdir) + "/data/progbook_inputs.xlsx", sheet_name="scenarios"
    )
    scen_name = list(pd.unique(scenarios.scenario_name))

    # Produce lists of countries within each region, this can be adjusted pretty easily now
    ref_countries = ref_countries[ref_countries["ISO3"].isin(countries)].reset_index(
        drop=True
    )
    top_10_burden = [
        "PAK",
        "IND",
        "CHN",
        "RUS",
        "USA",
        "IDN",
        "NGA",
        "UKR",
        "UZB",
        "BGD",
    ]

    AFR, AMR, EMR, EUR, SEAR, WPR, TOP = [], [], [], [], [], [], []
    for i in range(len(ref_countries)):
        if ref_countries.iloc[i, 1] == "AFR":
            AFR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 1] == "AMR":
            AMR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 1] == "EMR":
            EMR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 1] == "EUR":
            EUR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 1] == "SEAR":
            SEAR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 1] == "WPR":
            WPR.append(ref_countries.iloc[i, 0])
        if ref_countries.iloc[i, 0] in top_10_burden:
            TOP.append(ref_countries.iloc[i, 0])

    # Produce empty dataframes to collate results (fillna means can do iterative sums)
    epis_outs = ["HCV incidence", "HCC incidence", "HCV mortality", "DALYs"]

    # Economic aggregates
    econ_outs = [
        "PWID Ab Tests",
        "PWID RNA Tests",
        "Prisoner Ab Tests",
        "Prisoner RNA Tests",
        "Gen Ab Tests",
        "Gen RNA Tests",
        "Total Ab Tests",
        "Total RNA Tests",
        "PWID Ab Positive",
        "PWID RNA Positive",
        "Prisoner Ab Positive",
        "Prisoner RNA Positive",
        "Gen Ab Positive",
        "Gen RNA Positive",
        "Total Ab Positive",
        "Total RNA Positive",
        "Total treatment",
        "PWID vaccinations",
        "Prisoner vaccinations",
        "Gen vaccinations",
        "Total vaccinations",
        "Diagnosis Cost",
        "Disease Management Cost",
        "Treatment Cost",
        "Morbidity Cost",
        "Mortality Cost",
        "Direct Cost",
        "Productivity Cost",
        "Total Cost",
        "Vaccine Cost",
    ]

    sens_analyses_sheet_1 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="epi"
    )
    sens_analyses_sheet_2 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="econ"
    )
    sens_analyses_sheet_3 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="other"
    )
    sens_analyses = (
        list(sens_analyses_sheet_1.description.values)
        + list(sens_analyses_sheet_2.description.values)
        + list(sens_analyses_sheet_3.description.values)
    )
    for analysis in sens_analyses:
        epi_agg = {}
        for reg in who_reg:
            epi_agg[reg] = {}
            for scen in scen_name:
                epi_agg[reg][scen] = {}
                for out in epis_outs:
                    epi_agg[reg][scen][out] = pd.DataFrame(columns=["year", "central"])
                    epi_agg[reg][scen][out].year = np.arange(2000.5, 2051.1, 1)
                    epi_agg[reg][scen][out] = epi_agg[reg][scen][out].fillna(0.0)

        econ_agg = {}
        for reg in who_reg:
            econ_agg[reg] = {}
            for scen in scen_name:
                econ_agg[reg][scen] = {}
                for out in econ_outs:
                    econ_agg[reg][scen][out] = pd.DataFrame(
                        {
                            "year": np.arange(2000.5, 2051.1, 1),
                            "central": np.zeros(
                                len(np.arange(2000.5, 2051.1, 1)), dtype=float
                            ),
                        }
                    )
                    # econ_agg[reg][scen][out] = pd.DataFrame(columns=["year", "central"])
                    # econ_agg[reg][scen][out].year = np.arange(2000.5, 2051.1, 1)
                    # econ_agg[reg][scen][out] = econ_agg[reg][scen][out].fillna(0.0)

        for country in countries:
            country_data = sc.load(
                sens_folder / "agg_outputs" / analysis / f"{country}_econ_eval.pkl"
            )

            if country in AFR:
                adds = ["AFR", "global"]
            if country in AMR:
                adds = ["AMR", "global"]
            if country in EMR:
                adds = ["EMR", "global"]
            if country in EUR:
                adds = ["EUR", "global"]
            if country in SEAR:
                adds = ["SEAR", "global"]
            if country in WPR:
                adds = ["WPR", "global"]
            if country in TOP:
                adds += ["top10"]

            for add in adds:
                for scen in scen_name:
                    i = 0
                    # Epidemiological Outs
                    epi_agg[add][scen]["HCV incidence"].iloc[:, i + 1] = (
                        epi_agg[add][scen]["HCV incidence"].iloc[:, i + 1]
                        + country_data["agg_data"][scen]["HCV incidence"].iloc[:, i + 1]
                    )
                    epi_agg[add][scen]["HCC incidence"].iloc[:, i + 1] = (
                        epi_agg[add][scen]["HCC incidence"].iloc[:, i + 1]
                        + country_data["agg_data"][scen]["HCC incidence"].iloc[:, i + 1]
                    )
                    epi_agg[add][scen]["HCV mortality"].iloc[:, i + 1] = (
                        epi_agg[add][scen]["HCV mortality"].iloc[:, i + 1]
                        + country_data["agg_data"][scen]["HCV mortality"].iloc[:, i + 1]
                    )
                    epi_agg[add][scen]["DALYs"].iloc[:, i + 1] = (
                        epi_agg[add][scen]["DALYs"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["DALYs"].iloc[:, i + 1]
                    )

                    # Tests completed
                    econ_agg[add][scen]["PWID Ab Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID Ab Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["PWID Ab Tests"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["PWID RNA Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID RNA Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["PWID RNA Tests"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Prisoner Ab Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Prisoner Ab Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Prisoner Ab Tests"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Prisoner RNA Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Prisoner RNA Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Prisoner RNA Tests"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Gen Ab Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Gen Ab Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Gen Ab Tests"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Gen RNA Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Gen RNA Tests"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Gen RNA Tests"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Total Ab Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID Ab Tests"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Prisoner Ab Tests"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Gen Ab Tests"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Total RNA Tests"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID RNA Tests"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Prisoner RNA Tests"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Gen RNA Tests"].iloc[:, i + 1]
                    )

                    # Positive Tests
                    econ_agg[add][scen]["PWID Ab Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID Ab Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["PWID Ab Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["PWID RNA Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID RNA Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["PWID RNA Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Prisoner Ab Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Prisoner Ab Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Prisoner Ab Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Prisoner RNA Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Prisoner RNA Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Prisoner RNA Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Gen Ab Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Gen Ab Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Gen Ab Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Gen RNA Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Gen RNA Positive"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Gen RNA Positive"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Total Ab Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID Ab Positive"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Prisoner Ab Positive"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Gen Ab Positive"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Total RNA Positive"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID RNA Positive"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Prisoner RNA Positive"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Gen RNA Positive"].iloc[:, i + 1]
                    )

                    # Vaccinations and Treatments Administered
                    econ_agg[add][scen]["Total treatment"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Total treatment"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Total treatment"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["PWID vaccinations"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID vaccinations"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["PWID vaccinations"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Prisoner vaccinations"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Prisoner vaccinations"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Prisoner vaccinations"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Gen vaccinations"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Gen vaccinations"].iloc[:, i + 1]
                        + country_data["util_data"][scen]["Gen vaccinations"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Total vaccinations"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["PWID vaccinations"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Prisoner vaccinations"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Gen vaccinations"].iloc[:, i + 1]
                    )

                    # Costs
                    econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["Diagnosis Cost"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen][
                            "Disease Management Cost"
                        ].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["Treatment Cost"].iloc[
                            :, i + 1
                        ]
                    )
                    econ_agg[add][scen]["Morbidity Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Morbidity Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen][
                            "Morbidity Productivity"
                        ].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Mortality Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Mortality Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen][
                            "Mortality Productivity"
                        ].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Direct Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Productivity Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Morbidity Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Mortality Cost"].iloc[:, i + 1]
                    )
                    econ_agg[add][scen]["Total Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Diagnosis Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Disease Management Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Treatment Cost"].iloc[:, i + 1]
                        + econ_agg[add][scen]["Productivity Cost"].iloc[:, i + 1]
                    )

                    # Vaccination Cost
                    econ_agg[add][scen]["Vaccine Cost"].iloc[:, i + 1] = (
                        econ_agg[add][scen]["Vaccine Cost"].iloc[:, i + 1]
                        + country_data["econ_ests"][scen]["Vaccine Cost"].iloc[:, i + 1]
                    )

        filename = sc.makefilepath(
            filename="epi_agg.pkl",
            folder=sens_folder / "agg_epi_econ" / analysis,
            makedirs=True,
        )
        sc.save(filename=filename, obj=epi_agg)

        filename = sc.makefilepath(
            filename="econ_agg.pkl",
            folder=sens_folder / "agg_epi_econ" / analysis,
            makedirs=True,
        )
        sc.save(filename=filename, obj=econ_agg)


def econ_analysis_central(sens_folder):
    """Perform economic analysis based on sensitivity analyses data (central results only).
    
    This function reads scenario and sensitivity analysis data from Excel files, processes the data to compute aggregate outcomes for various regions, and saves the results in specified output files. It calculates epidemiological and economic outcomes, including cost-benefit ratios (BCR) for different vaccination scenarios.
    
    Args:
        sens_folder (Path): The folder path containing sensitivity analysis data.
    
    Returns:
        None: The function saves the results to Excel and CSV files in the specified output directory.
    """

    # Global and Regional Summary Table
    scenarios = pd.read_excel(
        str(rootdir) + "/data/progbook_inputs.xlsx", sheet_name="scenarios"
    )
    scen_name = list(
        pd.unique(scenarios.scenario_name)
    )  # needed to loop through for data extraction etc

    scen_cfs = scenarios.iloc[:, np.r_[0, 1, 5]].dropna().reset_index(drop=True)
    scens = list(scen_cfs.scenario_name)
    comps = list(scen_cfs.counterfactual)
    scen_label = list(pd.unique(scenarios.description))  # Comparators and Scenarios
    vax_label = scen_label[3:]  # remove counterfactuals
    regions = ["AFR", "AMR", "EMR", "EUR", "SEAR", "WPR", "global", "top10"]

    # Aggregate (sum) outcomes (by region)
    agg_reg_table = {}
    outcomes_epi = ["HCV incidence", "HCC incidence", "HCV mortality", "DALYs"]
    outcomes_econ = [
        "Total Ab Tests",
        "Total RNA Tests",
        "Total treatment",
        "Total vaccinations",
        "Diagnosis Cost",
        "Disease Management Cost",
        "Treatment Cost",
        "Morbidity Cost",
        "Mortality Cost",
        "Direct Cost",
        "Productivity Cost",
        "Total Cost",
        "Vaccine Cost",
    ]
    outcomes = outcomes_epi + outcomes_econ

    # pathlib.Path
    sens_analyses_sheet_1 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="epi"
    )
    sens_analyses_sheet_2 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="econ"
    )
    sens_analyses_sheet_3 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="other"
    )
    sens_analyses = (
        list(sens_analyses_sheet_1.description.values)
        + list(sens_analyses_sheet_2.description.values)
        + list(sens_analyses_sheet_3.description.values)
    )
    for analysis in sens_analyses:
        epi_agg = sc.load(sens_folder / "agg_epi_econ" / analysis / "epi_agg.pkl")
        econ_agg = sc.load(sens_folder / "agg_epi_econ" / analysis / "econ_agg.pkl")

        for reg in regions:
            agg_reg_table[reg] = pd.DataFrame(columns=["scenarios"] + outcomes)
            agg_reg_table[reg].scenarios = scen_name
            agg_reg_table[reg] = agg_reg_table[reg].set_index("scenarios", drop=True)

            for scen in scen_name:
                for epi in outcomes_epi:
                    # Limit from 2026-2050
                    cen_epi = np.round(
                        np.sum(epi_agg[reg][scen][epi].iloc[25:, 1], axis=0), 0
                    )
                    agg_reg_table[reg].at[scen, epi] = f"{cen_epi:,.0f}"
                    # agg_reg_table[reg].at[scen, epi] = f"{cen_epi/1e6:,.2f}"

                for econ in outcomes_econ:
                    cen_econ = np.sum(econ_agg[reg][scen][econ].iloc[:, 1], axis=0)
                    # if econ in ["Total vaccinations", "Vaccine Cost"]:
                    #     agg_reg_table[reg].at[scen, econ] = f"{cen_econ/1e6:,.2f}"
                    # else:
                    #     agg_reg_table[reg].at[scen, econ] = f"{cen_econ/1e9:,.2f}"
                    agg_reg_table[reg].at[scen, econ] = f"{cen_econ:,.0f}"

        # Vaccine Attributable Differences (by region)
        vad_reg_table = {}
        ben_cost = ["BCR"]
        bcr_denom = ["Vaccine Cost"]
        outcomes_vad = outcomes + ben_cost
        flip_list = ["Total vaccinations", "Vaccine Cost"]

        for reg in regions:
            vad_reg_table[reg] = pd.DataFrame(columns=["scenarios"] + outcomes_vad)
            vad_reg_table[reg].scenarios = scens
            vad_reg_table[reg] = vad_reg_table[reg].set_index("scenarios", drop=True)

            for i, scen in enumerate(scens):
                for epi in outcomes_epi:
                    cen_epi = np.sum(
                        epi_agg[reg][comps[i]][epi].iloc[25:, 1], axis=0
                    ) - np.sum(epi_agg[reg][scen][epi].iloc[25:, 1], axis=0)
                    vad_reg_table[reg].at[scen, epi] = f"{cen_epi:,.0f}"

                for econ in outcomes_econ:
                    if econ in flip_list:
                        cen_econ = np.sum(
                            econ_agg[reg][scen][econ].iloc[25:, 1], axis=0
                        ) - np.sum(econ_agg[reg][comps[i]][econ].iloc[25:, 1], axis=0)
                        vad_reg_table[reg].at[scen, econ] = f"{cen_econ:,.0f}"

                    else:
                        cen_econ = np.sum(
                            econ_agg[reg][comps[i]][econ].iloc[25:, 1], axis=0
                        ) - np.sum(econ_agg[reg][scen][econ].iloc[25:, 1], axis=0)
                        # if econ == "Productivity Cost":
                        #     vad_reg_table[reg].at[scen, econ] = f"{cen_econ/1e9:,.2f}"
                        # elif econ == "Direct Cost":
                        #     vad_reg_table[reg].at[scen, econ] = f"{cen_econ/1e6:,.2f}"
                        # else:
                        #     vad_reg_table[reg].at[scen, econ] = f"{cen_econ:,.0f}"
                        vad_reg_table[reg].at[scen, econ] = f"{cen_econ:,.0f}"

                for j, bcr in enumerate(ben_cost):
                    bcr_cent = (
                        np.sum(
                            econ_agg[reg][comps[i]]["Total Cost"].iloc[25:, 1], axis=0
                        )
                        - np.sum(econ_agg[reg][scen]["Total Cost"].iloc[25:, 1], axis=0)
                    ) / np.sum(econ_agg[reg][scen][bcr_denom[j]].iloc[25:, 1], axis=0)
                    vad_reg_table[reg].at[scen, bcr] = f"{bcr_cent:,.2f}"

        # Write to Excel spreadsheet and save
        outcomes_folder = sens_folder / "outcomes" / analysis
        outcomes_folder.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(
            outcomes_folder / "regional_aggregate_outcomes.xlsx"
        ) as writer:
            for sheet_name, df in agg_reg_table.items():
                df.to_excel(writer, sheet_name=sheet_name, index=True)
            print(
                "Aggregate outcome data table saved: {}".format(
                    outcomes_folder / "regional_aggregate_outcomes.xlsx"
                )
            )

        with pd.ExcelWriter(
            outcomes_folder / "regional_comparative_outcomes.xlsx"
        ) as writer:
            for sheet_name, df in vad_reg_table.items():
                df.to_excel(writer, sheet_name=sheet_name, index=True)
            print(
                "Aggregate outcome data table saved: {}".format(
                    outcomes_folder / "regional_comparative_outcomes.xlsx"
                )
            )

        # Country Level BCRs for mapping
        loop_folder = sens_folder / "agg_outputs" / analysis
        countries = [
            country.stem[:3] for country in loop_folder.iterdir() if country.is_file()
        ]
        map_bcr = "Vaccine Cost"

        country_bcr = pd.DataFrame(columns=["country"] + vax_label)
        country_bcr.country = countries
        country_bcr = country_bcr.set_index("country", drop=True)

        for country in countries:
            data = sc.load(
                sens_folder / "agg_outputs" / analysis / f"{country}_econ_eval.pkl"
            )

            for i, scen in enumerate(vax_label):
                country_bcr.at[country, scen] = np.round(
                    (
                        np.sum(data["econ_ests"][comps[i]]["Total Cost"].central)
                        - np.sum(data["econ_ests"][scens[i]]["Total Cost"].central)
                    )
                    / np.sum(data["econ_ests"][scens[i]][map_bcr].central),
                    3,
                )

        country_bcr.to_csv(outcomes_folder / "bcr_map.csv")
        print("BCR map data saved: {}".format(outcomes_folder / "bcr_map.csv"))


def write_sensitivity_table(scens_folder, sens_folder):
    """Writes a sensitivity analysis table to an Excel file.
    
    This function reads scenario and sensitivity analysis data from specified folders, processes the data, and writes the results to an Excel file. The output includes the Benefit-Cost Ratios (BCR) for various regions and scenarios, along with percentage differences from baseline values.
    
    Args:
        scens_folder (str or Path): The path to the folder containing scenario data.
        sens_folder (str or Path): The path to the folder where the sensitivity analysis results will be saved.
    
    Returns:
        None: The function saves the output directly to an Excel file and does not return any value.
    
    Raises:
        FileNotFoundError: If the specified Excel files cannot be found.
        ValueError: If the data in the Excel files is not in the expected format.
    """
    regions = ["global", "top10", "AFR", "AMR", "EMR", "EUR", "SEAR", "WPR"]
    regions_spelled = [
        "Global",
        "10 High Burden Countries",
        "African Region",
        "Region of the Americas",
        "Eastern Mediterranean Region",
        "European Region",
        "South East Asia Region",
        "Western Pacific Region",
    ]
    scenarios = pd.read_excel(
        str(rootdir) + "/data/progbook_inputs.xlsx", sheet_name="scenarios"
    )
    scen_name = list(pd.unique(scenarios.scenario_name))[3:]
    scen_desc = list(pd.unique(scenarios.description))[3:]
    sens_analyses_sheet_1 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="epi"
    )
    sens_analyses_sheet_2 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="econ"
    )
    sens_analyses_sheet_3 = pd.read_excel(
        str(rootdir) + "/data/sensitivity_analyses.xlsx", sheet_name="other"
    )
    sens_analyses = (
        list(sens_analyses_sheet_1.description.values)
        + list(sens_analyses_sheet_2.description.values)
        + list(sens_analyses_sheet_3.description.values)
    )

    with pd.ExcelWriter(sens_folder / "print_sensitivity.xlsx") as writer:
        for reg, region in zip(regions, regions_spelled):
            rows = []

            # Baseline first
            comp = pd.read_excel(
                scens_folder / "regional_comparative_outcomes.xlsx", sheet_name=reg
            )
            out_row = ["Baseline"]
            baseline_bcr = []
            for scen in scen_name:
                data = comp[comp.scenarios == scen]
                bcr = data["BCR_5USD"].values[0].split(" \n")[0]
                bcr = float(bcr)
                out_row = out_row + [f"{bcr:,.2f}"]
                baseline_bcr = baseline_bcr + [bcr]
            rows.append(out_row)

            # Sensitivity analyses
            for analysis in sens_analyses:
                comp = pd.read_excel(
                    sens_folder
                    / "outcomes"
                    / analysis
                    / "regional_comparative_outcomes.xlsx",
                    sheet_name=reg,
                )
                out_row = [analysis]
                for i, scen in enumerate(scen_name):
                    data = comp[comp.scenarios == scen]
                    bcr = data["BCR"].values[0]
                    diff = (bcr - baseline_bcr[i]) / baseline_bcr[i]
                    out_row = (
                        out_row + [f"{bcr:,.2f}\n(+{diff*100:,.1f}%)"]
                        if diff > 0
                        else out_row + [f"{bcr:,.2f}\n({diff*100:,.1f}%)"]
                    )
                rows.append(out_row)
            df_row = pd.DataFrame(
                rows, columns=["Benefit-Cost Ratio\n2026–2050"] + scen_desc
            )
            df_row.to_excel(writer, sheet_name=reg, index=False)

    print(
        "Print sensitivity analyses tables saved: {}".format(
            sens_folder / "print_sensitivity.xlsx"
        )
    )
