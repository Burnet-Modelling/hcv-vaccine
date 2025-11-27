import atomica as at
import numpy as np
import pandas as pd
import hcv.utils as ut

rootdir = ut.get_project_root()

#%% Parameters used throughout
pb_dir = rootdir / 'progbooks'
pb_dir.mkdir(parents=True, exist_ok=True)
pb_temp_dir = rootdir / 'progbook_templates'
pb_temp_dir.mkdir(parents=True, exist_ok=True)

income_cost = {
    'Low income': 'unit_cost_lic',
    'Lower middle income': 'unit_cost_lmc',
    'Upper middle income': 'unit_cost_umc',
    'High income': 'unit_cost_hic',
    }

df_income = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Info-Country Region Allocation')
country_income = dict()
country_region = dict()
for _, row in df_income.iterrows():
        country_income[row['ISO3']] = row['Income Group']
        country_region[row['ISO3']] = row['Region']

#%% Functions
 
# Calculate an Intervention value from relative risk and Baseline Value  
def rr_to_value(baseline, relative_risk):
    return baseline*relative_risk

# Calculate an Intervention value from Odds Ratio and Baseline Value  
def or_to_value(baseline, odds_ratio): 
    return odds_ratio*baseline/((1-baseline) + (odds_ratio*baseline))

# Calculate an Intervention value from probablity and Baseline Value  
def prob_to_value(baseline, probability, min_or_max):
    if min_or_max == 'min':
        val = min(baseline, probability)
    elif min_or_max == 'max':
        val = max(baseline, probability)
    return val

# Convert proportion linked-to-care to proportion lost-to-follow-up using relative risk
def ltc_rr_to_ltfu(baseline, relative_risk=1): 
    return max(0, 1 - (1 - baseline) * relative_risk) 

# Convert proportion linked-to-care to proportion lost-to-follow-up using odds ratio
def ltc_or_to_ltfu(baseline, odds_ratio): 
    return 1 - or_to_value(1 - baseline, odds_ratio) 

# Convert proportion linked-to-care to proportion lost-to-follow-up using probabilities
def ltc_prob_to_ltfu(baseline, probability):
     return min(baseline,1-probability)

#%% Generate progbook function
def generate_progbook(country, cal_folder=None, result=None, savedir=None, cov_scenario=False, cov=None):  
    
    if result is None:
        assert cal_folder is not None, 'cal_folder should not be empty if a result is not passed as an argument'
        P = ut.project(country, load_calibration=True, cal_folder=cal_folder, load_programs=False)
        parset = P.make_parset()
        parset.load_calibration(cal_folder/f'{country}_calibration.xlsx')   
        result = P.run_sim(parset=parset, result_name='Status-quo')
    if cov_scenario:
        assert cov is not None, 'A value of cov should be provided if cov_scenario is True (cov represents the extension for the progbook data sheet)'
        data_progbook = str(rootdir) + f'/data/progbook_inputs_coverage_{cov}.xlsx'
    else:
        data_progbook = str(rootdir) + '/data/progbook_inputs.xlsx'
    region = country_region[country]
    income_group = country_income[country]
    prog_years = [2026]
    sim_dt = result.dt
    stages = ['f{}'.format(s) for s in [0,1,2,3,4]] + ['dc', 'hcc']
    
    comps_ab_negative = ['s0_abneg']
    comps_susceptible = ['s{}'.format(s) for s in [0,1,2,3,4]] + ['sdc','shcc','s0_uw','s0_aw']
    comps_acute = ['{}_acute'.format(s) for s in stages]
    comps_lost_to_care = ['{}_ud_ltc'.format(s) for s in stages]
    comps_undiagnosed = ['{}_ab_ud'.format(s) for s in stages] + ['f0_acute_uw']
    comps_diagnosed_ab = ['{}_ab_d'.format(s) for s in stages] + ['f0_acute_aw']
    comps_diagnosed_rna = ['{}_pcr_d'.format(s) for s in stages]
    comps_on_treatment = ['{}_t1'.format(s) for s in stages] + ['{}_t2'.format(s) for s in stages] + ['{}_ft1'.format(s) for s in stages] + ['{}_ft2'.format(s) for s in stages]
    
    df_interventions = pd.read_excel(pd.ExcelFile(data_progbook), sheet_name='interventions')
    interventions = {}
    for code, label in zip(df_interventions['code'],df_interventions['label']):
        interventions[code] = label

    # save empty progbook
    F, D = ut.return_fw_db(country)
    P = at.Project(framework=F,databook=D, do_run=False)
    P.make_progbook(pb_temp_dir / f'progbook_{country}.xlsx',progs=interventions,data_start=prog_years[0],data_end=prog_years[-1])      
    #%% Parameters
    # load empty progbook
    P = at.ProgramSet.from_spreadsheet(spreadsheet=pb_temp_dir / f'progbook_{country}.xlsx', framework=F, data=D, _allow_missing_data=True)
    
    # targetable parameters
    all_pars = F.pars
    pars_ab_test = ['test_ab_f0f2_und_ic','test_ab_f3f4_und_ic','test_ab_f0f2_und_nic','test_ab_f3f4_und_nic']
    pars_pcr_test = ['test_pcr_f0f2_ab_ic','test_pcr_f3f4_ab_ic']
    pars_poc_test = ['test_pcr_f0f2_und_ic','test_pcr_f0f2_und_nic','test_pcr_f3f4_und_ic','test_pcr_f3f4_und_nic']
    pars_nic_test = ['test_ab_f0f2_und_nic','test_ab_f3f4_und_nic']
    pars_treat = ['tx_{}'.format(s) for s in stages]
    pars_vax_uptake = ['vaccine_uptake']
    pars_vax_infx = ['vaccine_infx_mult_pwid','vaccine_infx_mult_gen']
    
    # baseline year prior to programs start
    baseline_year = prog_years[0]
    # load programs info sheets
    df_target_pops = pd.read_excel(pd.ExcelFile(data_progbook), sheet_name='target_pops', index_col='code').T
    df_target_comps = pd.read_excel(pd.ExcelFile(data_progbook), sheet_name='target_comps', index_col='code').T
    df_spending = pd.read_excel(pd.ExcelFile(data_progbook), sheet_name='spending', index_col='code')
    df_effect = pd.read_excel(pd.ExcelFile(data_progbook), sheet_name='effects', index_col='code')
    
    #%% Loop over interventions
    pars_effect = []
    for intervention in interventions:
        # target populations
        target_pops = df_target_pops[df_target_pops[intervention]=='y'].index.tolist()
        P.programs[intervention].target_pops = target_pops
        
        # target compartments
        target_comps_overall = df_target_comps[df_target_comps[intervention]=='y'].index.tolist()
        target_comps = []
        if 'ab_negative' in target_comps_overall:
            target_comps += comps_ab_negative
        if 'susceptible' in target_comps_overall:
            target_comps += comps_susceptible
        if 'lost_to_care' in target_comps_overall:
            target_comps += comps_lost_to_care
        if 'acute' in target_comps_overall:
            target_comps += comps_acute
        if 'undiagnosed' in target_comps_overall:
            target_comps += comps_undiagnosed
        if 'diagnosed_ab' in target_comps_overall:
            target_comps += comps_diagnosed_ab
        if 'diagnosed_rna' in target_comps_overall:
            target_comps += comps_diagnosed_rna
        if 'on_treatment' in target_comps_overall:
            target_comps += comps_on_treatment
        P.programs[intervention].target_comps = target_comps
        
        # spending data
        units = df_spending.loc[intervention,'units']
        cost_par = income_cost[income_group]
        unit_cost = df_spending.loc[intervention,cost_par]
        P.programs[intervention].unit_cost = at.TimeSeries(assumption=unit_cost, 
                                                           units='$/person/year' if units=='continuous' else '$/person (one-off)')
        P.programs[intervention].spend_data = at.TimeSeries(assumption=0, 
                                                            units='$/year') # make initial spending a small, negligible but non-zero number for optimisation initialisation
        P.programs[intervention].capacity_constraint = at.TimeSeries(assumption=df_spending.loc[intervention,'capacity'] if pd.notnull(df_spending.loc[intervention,'capacity']) else '', 
                                                                     units='people' if units=='continuous' else 'people/year')
        P.programs[intervention].coverage = at.TimeSeries(assumption=df_spending.loc[intervention,'coverage'] if pd.notnull(df_spending.loc[intervention,'coverage']) else '', 
                                                          units='people' if units=='continuous' else 'people/year')
        P.programs[intervention].saturation = at.TimeSeries(assumption=df_spending.loc[intervention,'saturation'] if pd.notnull(df_spending.loc[intervention,'saturation']) else '', 
                                                            units='N.A.')
    
        # program effects
        pars_effect_intervention_all = df_effect.loc[intervention].dropna().index.tolist()
        pars_effect_intervention = [item for item in pars_effect_intervention_all if ('_measure' not in item) and ('interaction' not in item)]
        
        pars_pops=[]
        pars_all=[]
        for out in result.model._exec_order["transition_pars"]:
            pars_all.append((out.name, out.pop.name))
        for par in pars_treat:
            for pop in target_pops:
                pars_pops.append((par, pop))
        indices_match = [pars_all.index(item) for item in pars_pops]
    
        x = 'treat'
        if x in pars_effect_intervention:
            for par in pars_treat:
                for pop in target_pops:
                    # baseline value
                    baseline = result.get_variable(par, pop)[0].vals[list(result.t).index(baseline_year)]
                    if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                        baseline = baseline*sim_dt
                    elif all_pars['format'][par] == 'number':
                        index_parpop = indices_match[pars_pops.index((par, pop))]
                        index_time = list(result.t).index(baseline_year)
                        if result.model._exec_order["transition_pars"][index_parpop].source_popsize(index_time) != 0:
                            baseline = baseline*sim_dt\
                                /result.model._exec_order["transition_pars"][index_parpop].source_popsize(index_time)
                        else:
                            baseline = 0
                        if df_effect.loc[intervention,'{}_measure'.format(x)] == 'PROB':
                            effect_size = df_effect.loc[intervention,x]*sim_dt
                            effect = max(baseline,effect_size)
                    pars_effect.append({'par':par,'baseline':baseline,'imp_interaction':'','pop':pop,intervention:effect})
    
        x = 'diagnosis_ltc'
        if x in pars_effect_intervention:
            if 'poc_testing' in intervention:
                pars_nic_test += pars_poc_test
            for par in pars_nic_test:
                for pop in target_pops:
                    # baseline value
                    baseline = result.get_variable(par, pop)[0].vals[list(result.t).index(baseline_year)]
                    if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                        baseline = baseline*sim_dt
                    elif all_pars['format'][par] == 'number':
                        index_parpop = indices_match[pars_pops.index((par, pop))]
                        index_time = list(result.t).index(baseline_year)
                        baseline = baseline*sim_dt\
                            /result.model._exec_order["transition_pars"][index_parpop].source_popsize(index_time)
                    # effect size
                    if df_effect.loc[intervention,'{}_measure'.format(x)] == 'RR':
                        effect = rr_to_value(baseline, df_effect.loc[intervention,x])
                    elif df_effect.loc[intervention,'{}_measure'.format(x)] == 'OR':
                        effect = or_to_value(baseline, df_effect.loc[intervention,x])
                    elif df_effect.loc[intervention,'{}_measure'.format(x)] == 'PROB':
                        if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                            effect = prob_to_value(baseline, df_effect.loc[intervention,x]*sim_dt,'max')
                        else:
                            effect = prob_to_value(baseline, df_effect.loc[intervention,x],'max')
                    
                    pars_effect.append({'par':par,'baseline':baseline,'imp_interaction':'','pop':pop,intervention:effect})
    
        x = 'diagnosis_undiag'
        if x in pars_effect_intervention:
            if 'poc_testing' in intervention:
                pars_ab_test += pars_poc_test
            for par in pars_ab_test:
                for pop in target_pops:
                    # baseline value
                    baseline = result.get_variable(par, pop)[0].vals[list(result.t).index(baseline_year)]
                    if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                        baseline = baseline*sim_dt
                    elif all_pars['format'][par] == 'number':
                        index_parpop = indices_match[pars_pops.index((par, pop))]
                        index_time = list(result.t).index(baseline_year)
                        baseline = baseline*sim_dt\
                            /result.model._exec_order["transition_pars"][index_parpop].source_popsize(index_time)
                    # effect size
                    if df_effect.loc[intervention,'{}_measure'.format(x)] == 'RR':
                        effect = rr_to_value(baseline, df_effect.loc[intervention,x])
                    elif df_effect.loc[intervention,'{}_measure'.format(x)] == 'OR':
                        effect = or_to_value(baseline, df_effect.loc[intervention,x])
                    elif df_effect.loc[intervention,'{}_measure'.format(x)] == 'PROB':
                        if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                            effect = prob_to_value(baseline, df_effect.loc[intervention,x]*sim_dt,'max')
                        else:
                            effect = prob_to_value(baseline, df_effect.loc[intervention,x],'max')
                    
                    pars_effect.append({'par':par,'baseline':baseline,'imp_interaction':'','pop':pop,intervention:effect})
    
        x = 'diagnosis_diag_ab'
        if x in pars_effect_intervention:
            for par in pars_pcr_test:
                for pop in target_pops:
                    # baseline value
                    baseline = result.get_variable(par, pop)[0].vals[list(result.t).index(baseline_year)]
                    if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                        baseline = baseline*sim_dt
                    elif all_pars['format'][par] == 'number':
                        index_parpop = indices_match[pars_pops.index((par, pop))]
                        index_time = list(result.t).index(baseline_year)
                        baseline = baseline*sim_dt\
                            /result.model._exec_order["transition_pars"][index_parpop].source_popsize(index_time)
                    # effect size
                    if df_effect.loc[intervention,'{}_measure'.format(x)] == 'RR':
                        effect = rr_to_value(baseline, df_effect.loc[intervention,x])
                    elif df_effect.loc[intervention,'{}_measure'.format(x)] == 'OR':
                        effect = or_to_value(baseline, df_effect.loc[intervention,x])
                    elif df_effect.loc[intervention,'{}_measure'.format(x)] == 'PROB':
                        if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                            effect = prob_to_value(baseline, df_effect.loc[intervention,x]*sim_dt,'max')
                        else:
                            effect = prob_to_value(baseline, df_effect.loc[intervention,x],'max')
                    
                    pars_effect.append({'par':par,'baseline':baseline,'imp_interaction':'','pop':pop,intervention:effect})
                                                                    
        
        x = 'vax_uptake'
        if x in pars_effect_intervention:
            for par in pars_vax_uptake:
                for pop in target_pops:
                    # baseline value
                    baseline = result.get_variable(par, pop)[0].vals[list(result.t).index(baseline_year)]
                    if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                        baseline = baseline*sim_dt
                    elif all_pars['format'][par] == 'number':
                        index_parpop = indices_match[pars_pops.index((par, pop))]
                        index_time = list(result.t).index(baseline_year)
                        baseline = baseline*sim_dt\
                            /result.model._exec_order["transition_pars"][index_parpop].source_popsize(index_time)
                    # effect size
                    effect = df_effect.loc[intervention,x]
                    if isinstance(effect, float):
                        effect_size = effect
                    elif isinstance(effect, str):
                        effect_list = effect.split(", ")
                        if len(effect_list) == 1:
                            effect_list = effect_list[0].split(": ")
                            assert(effect_list[0] == 'Parameter'), "Cell not in the correct format: (Parameter: XXX) or (Sheet: XXX, Column: XXX, How: Country / Region / Income)"
                            effect_size = result.get_variable(effect_list[1],pop)[0].vals[list(result.t).index(baseline_year)]
                        else:
                            assert(effect_list[0][:5] == 'Sheet'), "Cell not in the correct format: (Parameter: XXX) or (Sheet: XXX, Column: XXX, How: Country / Region / Income)"
                            sheet = effect_list[0].split("Sheet: ")[-1]
                            col = effect_list[1].split("Column: ")[-1]
                            method = effect_list[2].split("How: ")[-1]
                            df = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name=sheet)
                            if method == 'Country':
                                df = df[df.ISO3 == country]
                            elif method == 'Region':
                                df = df[df.Region == region]
                            elif method == 'Income':
                                df = df[df.Income == income_group]
                            df = df[col].dropna()
                            effect_size = df.values[-1]
                        
                    if df_effect.loc[intervention,'{}_measure'.format(x)] == 'RR':
                        effect = rr_to_value(baseline, effect_size)
                    elif df_effect.loc[intervention,'{}_measure'.format(x)] == 'OR':
                        effect = or_to_value(baseline, effect_size)
                    elif df_effect.loc[intervention,'{}_measure'.format(x)] == 'PROB':
                        if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                            effect = prob_to_value(baseline, effect_size*sim_dt,'max')
                        else:
                            effect = prob_to_value(baseline, effect_size,'max')
                    
                    pars_effect.append({'par':par,'baseline':baseline,'imp_interaction':'','pop':pop,intervention:effect})

        x = 'vax_infx'
        if x in pars_effect_intervention:
            for par in pars_vax_infx:
                for pop in target_pops:
                    # baseline value
                    baseline = result.get_variable(par, pop)[0].vals[list(result.t).index(baseline_year)]
                    if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                        baseline = baseline*sim_dt
                    elif all_pars['format'][par] == 'number':
                        index_parpop = indices_match[pars_pops.index((par, pop))]
                        index_time = list(result.t).index(baseline_year)
                        baseline = baseline*sim_dt\
                            /result.model._exec_order["transition_pars"][index_parpop].source_popsize(index_time)
                    # effect size
                    effect = df_effect.loc[intervention,x]
                    if isinstance(effect, float):
                        effect_size = effect
                    elif isinstance(effect, str):
                        effect_list = effect.split(", ")
                        sheet = effect_list[0].split("Sheet: ")[-1]
                        col = effect_list[1].split("Column: ")[-1]
                        method = effect_list[2].split("How: ")[-1]
                        df = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name=sheet)
                        if method == 'Country':
                            df = df[df.ISO3 == country]
                        elif method == 'Region':
                            df = df[df.Region == region]
                        elif method == 'Income':
                            df = df[df.Income == income_group]
                        df = df[col].dropna()
                        effect_size = df.values[-1]
                        
                    if df_effect.loc[intervention,'{}_measure'.format(x)] == 'RR':
                        effect = rr_to_value(baseline, effect_size)
                    elif df_effect.loc[intervention,'{}_measure'.format(x)] == 'OR':
                        effect = or_to_value(baseline, effect_size)
                    elif df_effect.loc[intervention,'{}_measure'.format(x)] == 'PROB':
                        if all_pars['format'][par] == 'rate' or all_pars['format'][par] == 'probability':
                            effect = prob_to_value(baseline, effect_size*sim_dt,'max')
                        else:
                            effect = prob_to_value(baseline, effect_size,'max')
                    
                    pars_effect.append({'par':par,'baseline':baseline,'imp_interaction':'','pop':pop,intervention:effect})

    # P.programs[intervention].spend_data = at.TimeSeries(assumption=1e-16, units='$/year') # make 1 initial spending a small, negligible but non-zero number for optimisation initialisation
                    
    # write to progbook
    df_effect = pd.DataFrame(pars_effect)
    df = df_effect.groupby(['par','pop'],as_index=False).first()
    excluded_columns = ['par','pop','baseline','imp_interaction']
    prog_columns = [col for col in df.columns if col not in excluded_columns]
    df_prog = df[prog_columns]
    
    for i in np.arange(len(df)):
        row = df.iloc[i]
        par = row.loc['par']
        pop = row.loc['pop']
        baseline = row.loc['baseline']
        imp_interaction = row.loc['imp_interaction']
        progs = {}
        row_prog = df_prog.iloc[i]
        for prog in df_prog.columns:
            if not np.isnan(row_prog[prog]):
                progs[prog] = row_prog[prog]
        P.covouts[(par, pop)] = at.programs.Covout(par=par,pop=pop,cov_interaction='random',imp_interaction=imp_interaction,baseline=baseline,progs=progs)
    
    if savedir is None:
        savedir = pb_dir / f'progbook_{country}.xlsx'
    P.save(savedir)
    
    # return pb_path
