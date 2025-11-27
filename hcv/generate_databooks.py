import numpy as np
import atomica as at
import pandas as pd
import atomica.model
atomica.model.model_settings["tolerance"] = 2e-6
import hcv
import itertools
import hcv.utils as ut

rootdir = ut.get_project_root()

def generate_databook(country, savedir=None):
    print(country)
    
    F = at.ProjectFramework(str(rootdir) + '/framework/hcv_vaccine_framework.xlsx')
    D = at.ProjectData.new(framework = F, tvec = np.arange(2000,2051), pops=hcv.default_pops, transfers=0)
    
    who_regions = pd.read_excel(str(rootdir)+"/data/flat_datasheet.xlsx", sheet_name="Cost - YLL and productivity").iloc[:,np.r_[1,4]]
    region = who_regions[who_regions.ISO3==country].WHO_reg.values[0]
    
    ## Demographic Parameters
    data1 = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name= 'Country-Demographic')
    data1 = data1[data1.ISO3 == country]

    data_dict = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Info-Demographic')
    data_dict = data_dict[data_dict['In Databook'] == True]

    for idx, row in data_dict.iterrows():
        if row['Population'] == 'Total':
            for pop in D.pops:
                D.tdve[row['Code Name']].ts[pop].t = list(data1[~pd.isna(data1[row['Display Name']])].Year)
                D.tdve[row['Code Name']].ts[pop].vals = list(data1[~pd.isna(data1[row['Display Name']])][row['Display Name']])
        else:
            D.tdve[row['Code Name']].ts[row['Population']].t = list(data1[~pd.isna(data1[row['Display Name']])].Year)
            D.tdve[row['Code Name']].ts[row['Population']].vals = list(data1[~pd.isna(data1[row['Display Name']])][row['Display Name']])

    ### Population Data
    for pop in D.pops:
        if 2000 not in D.tdve['alive'].ts[pop].t:
            D.tdve['alive'].ts[pop].t.append(2000)
            if len(D.tdve['alive'].ts[pop].vals) == 0:
                D.tdve['alive'].ts[pop].vals.append(0)
            else:
                D.tdve['alive'].ts[pop].vals.append(D.tdve['alive'].ts[pop].vals[0])
        if pop not in ['0-9_males', '0-9_females']: # Set the population entry rate of the other pops to 0
            D.tdve['birth_nb'].ts[pop].assumption = 0
    for pop in ['PWID_males','PWID_females']:
        if 2010 not in D.tdve['alive'].ts[pop].t:
            D.tdve['alive'].ts[pop].t.append(2010)
            if len(D.tdve['alive'].ts[pop].vals) == 1:
                D.tdve['alive'].ts[pop].vals.append(D.tdve['alive'].ts[pop].vals[0])
            else:
                D.tdve['alive'].ts[pop].vals.append(np.average(D.tdve['alive'].ts[pop].vals))

    
    # Delete not-used interactions
    interactions_idu = list(zip(hcv.default_pops_inter,hcv.default_pops_inter))
    sexes = ['males', 'females']
    for sx1, sx2 in itertools.product(sexes, sexes):
        interactions_idu += [('PWID_' + sx1, 'PWID_' + sx2)]

    for key in [keys for keys in D.interpops[0].ts.keys() if keys not in interactions_idu]: #list of keys excluding idu transfers
        del D.interpops[0].ts[key]
    
    # gen_to = ['18-64_males','18-64_females','65+_males','65+_females']
    # interactions_gen = []
    # for pop_from in D.pops:
    #     for pop_to in gen_to:
    #         interactions_gen += [(pop_from,pop_to)]
            
    # for key in [keys for keys in D.interpops[1].ts.keys() if keys not in interactions_gen]: #list of keys excluding idu transfers
    #     del D.interpops[1].ts[key]
        
    # Generate all transfers
    for idx, transfer_name in enumerate(hcv.default_transfer_codes.keys()): #adds transfer categories based on hcv.default_transfer_codes
        D.add_transfer(transfer_name,hcv.default_transfer_codes[transfer_name],pop_type='human')
        for pair in hcv.default_transfers[transfer_name]: #switches to Y for chosen pops
            D.transfers[idx].ts.append((pair[0],pair[1]), at.TimeSeries())

    ### Default Transfer Values
    # data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name = 'Global-Transfers')

    transfer_order = {'age':0, 'idu_status':1, 'inc':2}
    # for idx, val in data['Value'].items():
    #     D.transfers[transfer_order[data.loc[idx,'Parameter code']]].ts[data.loc[idx,'From'], data.loc[idx,'To']].assumption = data.loc[idx,'Value']
    #     D.transfers[transfer_order[data.loc[idx,'Parameter code']]].ts[data.loc[idx,'From'], data.loc[idx,'To']].units = data.loc[idx,'Units']

    ### IDU transfers    
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name = 'Country-IDU Transfers')
    data = data[data.ISO3 == country]
    D.transfers[1].ts['18-64_males','PWID_males'].units = 'Rate (per year)'
    D.transfers[1].ts['18-64_males','PWID_males'].assumption = np.mean(data['Rate in males'].values)
    D.transfers[1].ts['PWID_males','18-64_males'].units = 'Rate (per year)'
    D.transfers[1].ts['PWID_males','18-64_males'].assumption = np.mean(data['Rate out'].values)
    D.transfers[1].ts['18-64_females','PWID_females'].units = 'Rate (per year)'
    D.transfers[1].ts['18-64_females','PWID_females'].assumption = np.mean(data['Rate in females'].values)
    D.transfers[1].ts['PWID_females','18-64_females'].units = 'Rate (per year)'
    D.transfers[1].ts['PWID_females','18-64_females'].assumption = np.mean(data['Rate out'].values)
    
    ### Prison transfers    
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name = 'Country-Prison Transfers')
    data = data[data.ISO3 == country]
    prison_transfer_in_map = { 
        ('18-64_males', 'Prisoners_males'): 'Male admission rate (general)', 
        ('PWID_males', 'Prisoners_males'): 'Male admission rate (PWID)', 
        ('18-64_females', 'Prisoners_females'): 'Female admission rate (general)', 
        ('PWID_females', 'Prisoners_females'): 'Female admission rate (PWID)'}
    prison_transfer_out_map = { 
        ('Prisoners_males', '18-64_males'): 'Male release rate (general)', 
        ('Prisoners_males', 'PWID_males'): 'Male release rate (PWID)', 
        ('Prisoners_females', '18-64_females'): 'Female release rate (general)', 
        ('Prisoners_females', 'PWID_females'): 'Female release rate (PWID)' }
    
    for (pop_from, pop_to), col_label in prison_transfer_in_map.items():
        D.transfers[2].ts[(pop_from, pop_to)].units = 'Rate (per year)'
        D.transfers[2].ts[(pop_from, pop_to)].t = data.loc[data[col_label].fillna(0)!=0,'Year'].tolist()
        D.transfers[2].ts[(pop_from, pop_to)].vals = data.loc[data[col_label].fillna(0)!=0,col_label].tolist()
        
    # Estimate exit rate based on entry rate and pop size for countries with entry but no exit data (note this will be calibrated any way)
    for (pop_from_out, pop_to_out), (pop_from_in, pop_to_in) in zip(list(prison_transfer_out_map.keys()),list(prison_transfer_in_map.keys())):        
        alive_prison_t = D.tdve['alive'].ts[pop_from_out].t
        alive_prison_vals = D.tdve['alive'].ts[pop_from_out].vals
        alive_gpop_t = D.tdve['alive'].ts[pop_to_out].t
        alive_gpop_vals = D.tdve['alive'].ts[pop_to_out].vals
        if (D.transfers[2].ts[(pop_from_out, pop_to_out)].vals == []) and (D.transfers[2].ts[(pop_from_in, pop_to_in)].vals != []):
            for idx, year in enumerate(D.transfers[2].ts[(pop_from_in, pop_to_in)].t):
                idx_prison = min(range(len(alive_prison_t)), key=lambda i: abs(alive_prison_t[i]-year)) # find index of pop size value closest to "year"
                idx_gpop = min(range(len(alive_gpop_t)), key=lambda i: abs(alive_gpop_t[i]-year)) # find index of pop size value closest to "year"
                D.transfers[2].ts[(pop_from_out, pop_to_out)].t.append(year)
                flow_out = D.transfers[2].ts[(pop_from_in, pop_to_in)].vals[idx]*alive_gpop_vals[idx_gpop]/max(alive_prison_vals[idx_prison],1e-6) # this assumes a stable prison pop size
                D.transfers[2].ts[(pop_from_out, pop_to_out)].vals.append(flow_out)

    for (pop_from, pop_to), col_label in prison_transfer_in_map.items():
        D.transfers[2].ts[(pop_from, pop_to)].units = 'Rate (per year)'
        if not D.transfers[2].ts[(pop_from, pop_to)].vals:
            D.transfers[2].ts[(pop_from, pop_to)].assumption = 0
        else:
            D.transfers[2].ts[(pop_from, pop_to)].assumption = np.mean(D.transfers[2].ts[(pop_from, pop_to)].vals)
        D.transfers[2].ts[(pop_from, pop_to)].t = []
        D.transfers[2].ts[(pop_from, pop_to)].vals = []

    for (pop_from, pop_to), col_label in prison_transfer_out_map.items():
        D.transfers[2].ts[(pop_from, pop_to)].units = 'Rate (per year)'
        if not D.transfers[2].ts[(pop_from, pop_to)].vals:
            D.transfers[2].ts[(pop_from, pop_to)].assumption = 0
        else:
            D.transfers[2].ts[(pop_from, pop_to)].assumption = np.mean(D.transfers[2].ts[(pop_from, pop_to)].vals)
        D.transfers[2].ts[(pop_from, pop_to)].t = []
        D.transfers[2].ts[(pop_from, pop_to)].vals = []
            
    
    ### Migration
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Data-Net Migration')
    data = data[data.ISO3 == country]
    pops_map = {
                '0-9_males': 'Male population aged 0-9 years',
                '0-9_females': 'Female population aged 0-9 years',
                '10-17_males': 'Male population aged 10-17 years',
                '10-17_females': 'Female population aged 10-17 years',
                '18-64_males': 'Male population aged 18-64 years',
                '18-64_females': 'Female population aged 18-64 years',
                '65+_males': 'Male population aged 65+ years',
                '65+_females':'Female population aged 65+ years'
                }
    for pop_code,pop_label in pops_map.items():
        D.tdve['mig_r'].ts[pop_code].t = data['Year'].values
        D.tdve['mig_r'].ts[pop_code].vals = data[pop_label].values
        
    for pop_code in ['PWID_males','PWID_females','Prisoners_males','Prisoners_females']:
        D.tdve['mig_r'].ts[pop_code].assumption = 0
    
    ## Global Parameter Values
    ### Disease Progression
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Global-Disease Progression')

    for idx, row in data.iterrows():
        if row['Parameter Code Name'] == 'prop_clear':
            if row['Population(s)'] == 'female all':
                for pop in hcv.general_populations:
                    D.tdve[row['Parameter Code Name']].ts[pop + '_females'].assumption = row['Default Value']
            elif row['Population(s)'] == 'male all':
                for pop in hcv.general_populations:
                    D.tdve[row['Parameter Code Name']].ts[pop + '_males'].assumption = row['Default Value']
        else:
            for pop in D.pops:
                D.tdve[row['Parameter Code Name']].ts[pop].assumption = row['Default Value']

    ### Other Global Parameters
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Global-Other')
    for idx, row in data.iterrows():
        pop_groups = [populations.strip() for populations in row['Population(s)'].split(',')]
        if pop_groups == ['all']:
            iter_pops = D.pops
        else:
            iter_pops = pop_groups
        if (row['Parameter Code Name'] in ['treat_dur', 'prop_svr']):
            for pop in iter_pops:
                D.tdve[row['Parameter Code Name']].ts[pop].t.append(row['Year'])
                D.tdve[row['Parameter Code Name']].ts[pop].vals.append(row['Default Value'])
        else:
            for pop in iter_pops:
                D.tdve[row['Parameter Code Name']].ts[pop].assumption = row['Default Value']
                
    ### Care Cascade
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Country-PWID Care Cascade')
    data = data[data.ISO3 == country]
    # data_ever_tested = data.dropna(subset=['Ever tested (proportion)'])
    data_12mo_tested = data.dropna(subset=['Tested past 12mo (proportion)'])
    test_year = data_12mo_tested['Year'].values[0]
    test_val = data_12mo_tested['Tested past 12mo (proportion)'].values[0]
    if (data_12mo_tested['Tested 12mo data type'].values[0] == 'Regional estimate') and (region != 'AFR'):
        test_val = test_val/2
    for pop in ['PWID_males', 'PWID_females']:
        D.tdve['test_ab_f0f2_und_ic_1'].ts[pop].t = [2000, test_year]
        D.tdve['test_ab_f0f2_und_ic_1'].ts[pop].vals = [test_val/2, test_val]
        D.tdve['test_ab_f3f4_und_ic_1'].ts[pop].t = [2000, test_year]
        D.tdve['test_ab_f3f4_und_ic_1'].ts[pop].vals = [test_val/2, test_val]
    

    ### Info-Calibrated Parameters
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Info-Calibrated Pars')
    hr_prison_data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Data-Prisons Harm Reduction')
    hr_prison_data = hr_prison_data[hr_prison_data.ISO3 == country]
    infx_gen_data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Data-gPop Transmission')
    infx_gen_data = infx_gen_data[infx_gen_data.ISO3 == country]
    
    for idx, row in data.iterrows():
        par = row['Parameter Code Name']
        value = row['Default Value']
        # defined populations
        pop_groups = [populations.strip() for populations in row['Population(s)'].split(',')]
        if pop_groups == ['all']:
            iter_pops = D.pops
        else:
            iter_pops = pop_groups
        # defined countries
        country_groups = [countries.strip() for countries in row['Countries'].split(',')]
        if (country_groups != ['all']) and (country not in country_groups):
            continue            
        for pop in iter_pops:
            if ('test_' in par) and ('Prisoners_' in pop):
                if (hr_prison_data.empty) or (hr_prison_data['Any HR in prisons?'].values[0] == False):
                    D.tdve[par].ts[pop].assumption = 0
                else:
                    if 'pcr' in par:
                        D.tdve[par].ts[pop].assumption = D.tdve[par].ts['PWID_males'].assumption
                    else:
                        D.tdve[par].ts[pop].t = D.tdve[par].ts['PWID_males'].t
                        D.tdve[par].ts[pop].vals = D.tdve[par].ts['PWID_males'].vals
            elif par == 'infx_gen':
                D.tdve[par].ts[pop].t = list(infx_gen_data.Year)
                D.tdve[par].ts[pop].vals = list(infx_gen_data.infx_gen)
            else:
                D.tdve[par].ts[pop].assumption = value

    for pop in D.pops:
        if not D.tdve['infx_primary'].ts[pop].assumption:
            D.tdve['infx_primary'].ts[pop].assumption = 0 # Set FOI to zero for remaining populations
        if not D.tdve['infx_gen'].ts[pop].assumption and not D.tdve['infx_gen'].ts[pop].vals:
            D.tdve['infx_gen'].ts[pop].assumption = 0 # Set FOI to zero for remaining populations

    ### Disease burden
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Country-Disease Burden')
    data = data[data.ISO3 == country]

    data_dict = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Info-Disease Burden')
    data_dict = data_dict[data_dict['In Databook'] == True]

    for idx, row in data_dict.iterrows():
        if row['Population'] == 'PWID all':
            for pop in ['PWID_males', 'PWID_females']:
                D.tdve[row['Code Name']].ts[pop].t = list(data[~pd.isna(data[row['Display Name']])].Year)
                D.tdve[row['Code Name']].ts[pop].vals = list(data[~pd.isna(data[row['Display Name']])][row['Display Name']])
        elif row['Population'] == 'Prisoners all':
            for pop in ['Prisoners_males', 'Prisoners_females']:
                D.tdve[row['Code Name']].ts[pop].t = list(data[~pd.isna(data[row['Display Name']])].Year)
                D.tdve[row['Code Name']].ts[pop].vals = list(data[~pd.isna(data[row['Display Name']])][row['Display Name']])
        else:
            D.tdve[row['Code Name']].ts[row['Population']].t = list(data[~pd.isna(data[row['Display Name']])].Year)
            D.tdve[row['Code Name']].ts[row['Population']].vals = list(data[~pd.isna(data[row['Display Name']])][row['Display Name']])

    ### Disease Burden Estimation
    data2 = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx',
                         sheet_name='Country-Prison Prevalence')
    data2 = data2[data2.ISO3 == country]

    for pop in ['Prisoners_males', 'Prisoners_females']:  # Burden Estimates (prisoners) (if empty)
        sx = pop[10:-1]
        if len(D.tdve['prevalence'].ts[pop].vals) == 0:
            D.tdve['prevalence'].ts[pop].t = [2000]
            D.tdve['prevalence'].ts[pop].vals = [data2[f'Estimated {sx.title()} Prisoner HCV Prevalence'].values[0]]
            D.tdve['chronic'].ts[pop].t = [2000]
            D.tdve['chronic'].ts[pop].vals = [D.tdve['prevalence'].ts[pop].vals[0] * D.tdve['alive'].ts[pop].vals[0]]

    for pop in ['PWID_males', 'PWID_females']:  # Burden Estimates (prisoners) (if empty)
        if len(D.tdve['prevalence'].ts[pop].vals) == 0:
            D.tdve['prevalence'].ts[pop].t = [2000]
            D.tdve['prevalence'].ts[pop].vals = [data2['HCV RNA prevalence (PWID)'].values[0]]
            D.tdve['chronic'].ts[pop].t = [2000]
            D.tdve['chronic'].ts[pop].vals = [D.tdve['prevalence'].ts[pop].vals[0] * D.tdve['alive'].ts[pop].vals[0]]


    for par in ['chronic', 'prevalence']:
        for pop in D.pops:
            if len(D.tdve[par].ts[pop].vals) == 0:
                D.tdve[par].ts[pop].t = [2000]
                D.tdve[par].ts[pop].vals = [0]
            elif 2000 not in D.tdve[par].ts[pop].t:
                D.tdve[par].ts[pop].t.append(2000)
                D.tdve[par].ts[pop].vals.append(D.tdve[par].ts[pop].vals[0])

    ### Overall Indicators
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Country-Overall Indicators')
    data = data[data.ISO3 == country]

    data_dict = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Info-Overall Indicators')
    data_dict = data_dict[data_dict['In Databook'] == True]

    for idx, row in data_dict.iterrows():
        if row['Code Name'] == 'percent_diagnosed':
            for pop in D.pops:
                D.tdve[row['Code Name']].ts[pop] = at.TimeSeries(t=[2000], vals=[0], units='Fraction')
            D.tdve[row['Code Name']].ts['Total'] = at.TimeSeries(t=[2000], vals=[0], units='Fraction')
        elif row['Code Name'] != 'prop_undiag':
            for pop in D.pops:
                D.tdve[row['Code Name']].ts[pop] = at.TimeSeries(t=[2000], vals=[0], units='Number')
            D.tdve[row['Code Name']].ts['Total'] = at.TimeSeries(t=[2000], vals=[0], units='Number')


        if row['Data Available'] == True:
            if (len(data[~pd.isna(data[row['Display Name']])]) > 0) and (row['Code Name'] != 'prop_undiag'):
                D.tdve[row['Code Name']].ts['Total'].t = list(data[(data.ISO3 == country) & (~pd.isna(data[row['Display Name']]))].Year)
                D.tdve[row['Code Name']].ts['Total'].vals = list(
                    data[(data.ISO3 == country) & (~pd.isna(data[row['Display Name']]))][row['Display Name']])
                if any(isinstance(v,str) for v in D.tdve[row['Code Name']].ts['Total'].vals):
                    D.tdve[row['Code Name']].ts['Total'].t = [2000]
                    D.tdve[row['Code Name']].ts['Total'].vals = [10]
                if 2000 not in D.tdve[row['Code Name']].ts['Total'].t:
                    D.tdve[row['Code Name']].ts['Total'].t.append(2000)
                    if row['Code Name'] == 'percent_diagnosed':
                        D.tdve[row['Code Name']].ts['Total'].vals.append( D.tdve[row['Code Name']].ts['Total'].vals[0]/3)
                    elif row['Code Name'] == 'treat_total' and D.tdve[row['Code Name']].ts['Total'].t[0] > 2014:
                        D.tdve[row['Code Name']].ts['Total'].vals.append( D.tdve[row['Code Name']].ts['Total'].vals[0]/3)
                        D.tdve[row['Code Name']].ts['Total'].t.append(2014)
                        D.tdve[row['Code Name']].ts['Total'].vals.append( D.tdve[row['Code Name']].ts['Total'].vals[0]/3)
                    else:
                        D.tdve[row['Code Name']].ts['Total'].vals.append( D.tdve[row['Code Name']].ts['Total'].vals[0])
            elif row['Code Name'] == 'prop_undiag': # handle par prop_undiag differently as this is just for initialising pop groups
                for pop in D.pops:
                    D.tdve[row['Code Name']].ts[pop].t = [2000]
                    D.tdve[row['Code Name']].ts[pop].vals = [list(data[(data.ISO3 == country) & (~pd.isna(data[row['Display Name']]))][row['Display Name']])[0]]
                
    
    for pop in D.pops:
        D.tdve['treat_overall'].ts[pop].t = D.tdve['treat_total'].ts['Total'].t
        D.tdve['treat_overall'].ts[pop].vals = D.tdve['treat_total'].ts['Total'].vals

    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Country-Age Transfers')
    data = data[data.ISO3 == country]

    data_dict = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Info-Country Transfers')
    data_dict = data_dict[data_dict['In Databook'] == True]

    for idx, row in data_dict.iterrows():
        D.transfers[transfer_order[row['Code Name']]].ts[row['From Population'], row['To Population']].units = 'Number (years)'
        D.transfers[transfer_order[row['Code Name']]].ts[row['From Population'], row['To Population']].t = list(data[~pd.isna(data[row['Display Name']])].Year)
        D.transfers[transfer_order[row['Code Name']]].ts[row['From Population'], row['To Population']].vals = list(
            data[~pd.isna(data[row['Display Name']])][row['Display Name']])
        
    ### Future Projections
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Data-Future Projections')
    data = data[data.ISO3 == country]

    data_dict = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Info-Future Projections')
    data_dict = data_dict[data_dict['In Databook'] == True]

    for idx, row in data_dict.iterrows():
        if row['Code Name'] == 'age':
            D.transfers[transfer_order[row['Code Name']]].ts[row['From Population'], row['To Population']].t += list(data[~pd.isna(data[row['Display Name']])].Year)
            D.transfers[transfer_order[row['Code Name']]].ts[row['From Population'], row['To Population']].vals += list(
            data[pd.notna(data[row['Display Name']])][row['Display Name']])
        else:
            D.tdve[row['Code Name']].ts[row['Population']].t += list(data[~pd.isna(data[row['Display Name']])].Year)
            D.tdve[row['Code Name']].ts[row['Population']].vals += list(data[~pd.isna(data[row['Display Name']])][row['Display Name']])
    
    # Reinitialise total hcv (chronic+acute), noting that prevalence = total_hcv/alive (not chronic/alive), and assuming acute = 0
    
    for pop in D.pops:        
        D.tdve['total_hcv'].ts[pop].t = [2000]
        D.tdve['total_hcv'].ts[pop].vals = [D.tdve['alive'].ts[pop].vals[0] * max(min(D.tdve['prevalence'].ts[pop].vals[0],1),0)]
        D.tdve['chronic'].ts[pop].t = [2000]
        D.tdve['chronic'].ts[pop].vals = [D.tdve['alive'].ts[pop].vals[0] * max(min(D.tdve['prevalence'].ts[pop].vals[0],1),0)]
        D.tdve['acute'].ts[pop].t = [2000]
        D.tdve['acute'].ts[pop].vals = [0]
        D.tdve['prevalence'].ts[pop].t.append(2010)
        D.tdve['prevalence'].ts[pop].vals.append(D.tdve['prevalence'].ts[pop].vals[0])
    
    # Note if sum of plhcv in key pop groups != total hcv, then distribute the rest of plhcv to the remaining 18+ pop groups
    gpop = ['18-64_males','18-64_females','65+_males','65+_females']
    gpop_size = sum([D.tdve['alive'].ts[pop].vals[0] for pop in gpop])
    plhcv_diff = D.tdve['plhcv_total'].ts['Total'].vals[0] - sum([D.tdve['chronic'].ts[pop].vals[0] for pop in D.pops])
    if plhcv_diff > 0:
        for pop in gpop:
            if D.tdve['total_hcv'].ts[pop].vals[0] == 0:
                D.tdve['total_hcv'].ts[pop].vals[0] = plhcv_diff*D.tdve['alive'].ts[pop].vals[0]/gpop_size
                D.tdve['chronic'].ts[pop].vals[0] = plhcv_diff*D.tdve['alive'].ts[pop].vals[0]/gpop_size
                D.tdve['prevalence'].ts[pop].vals[0] = D.tdve['total_hcv'].ts[pop].vals[0]/D.tdve['alive'].ts[pop].vals[0]
                
    for pop in D.pops:
        D.tdve['prev'].ts[pop].t = D.tdve['prevalence'].ts[pop].t # intermediate parameter for calibration
        D.tdve['prev'].ts[pop].vals = D.tdve['prevalence'].ts[pop].vals

    
    # calculate prop treat according to plhcv distribution
    total_hcv = 0
    if (hr_prison_data.empty) or (hr_prison_data['Any HR in prisons?'].values[0] == False):
        pop_groups = [pop for pop in D.pops if 'Prisoners' not in pop]
        D.tdve['prop_treat_pre2015_pwid'].ts['Prisoners_males'].assumption = 0
        D.tdve['prop_treat_pre2015_pwid'].ts['Prisoners_females'].assumption = 0
        D.tdve['prop_treat_post2015_pwid'].ts['Prisoners_males'].assumption = 0
        D.tdve['prop_treat_post2015_pwid'].ts['Prisoners_females'].assumption = 0
    else:
        pop_groups = D.pops
    for pop in pop_groups:
        total_hcv += D.tdve['total_hcv'].ts[pop].vals[0]
    prop_weighted = dict()
    for pop in pop_groups:
        if ('PWID' in pop) or ('Prisoners' in pop):
            treat_weight = 1
        else:
            treat_weight = 1
        prop_weighted[pop] = D.tdve['total_hcv'].ts[pop].vals[0]/total_hcv*treat_weight
    total_weights = np.sum(list(prop_weighted.values()))
    prop = dict()
    for pop in pop_groups:        
        prop[pop] = prop_weighted[pop]/total_weights
        if ('PWID' in pop) or ('Prisoners' in pop):
            par = 'prop_treat_pre2015_pwid'
            D.tdve[par].ts[pop].t = [2000,2014,2015]
            D.tdve[par].ts[pop].vals = [prop[pop],prop[pop],0]
            par = 'prop_treat_post2015_pwid'
            D.tdve[par].ts[pop].t = [2000,2014,2015]
            D.tdve[par].ts[pop].vals = [0,0,prop[pop]]
        else:
            par = 'prop_treat_pre2015'
            D.tdve[par].ts[pop].t = [2000,2014,2015]
            D.tdve[par].ts[pop].vals = [prop[pop],prop[pop],0]
            par = 'prop_treat_post2015'
            D.tdve[par].ts[pop].t = [2000,2014,2015]
            D.tdve[par].ts[pop].vals = [0,0,prop[pop]]
            
    ### Model Initialisation
    stages = ['f0', 'f1', 'f2', 'f3', 'f4', 'dc', 'hcc']
    pops_burden = gpop + ['PWID_males','PWID_females','Prisoners_males','Prisoners_females']
    for pop in D.pops:
        if abs(1 - sum([D.tdve['prop_f0f2'].ts[pop].assumption/3 if j in ['f0', 'f1', 'f2'] else D.tdve['prop_' + j].ts[pop].assumption for j in stages])) > atomica.model.model_settings['tolerance']:
            print('ERROR: Disease stage proportions do not add up to 1')
        
        D.tdve['total_hcv'].ts[pop].t = [2000]
        D.tdve['total_hcv'].ts[pop].vals = [D.tdve['alive'].ts[pop].vals[0] * max(min(D.tdve['prevalence'].ts[pop].vals[0],1),0)]
        
        # initialize susception populations ('sus_all' = never infected + past infections, 'sus' = past infections only)
        if pop not in pops_burden:
            D.tdve['sus_all'].ts[pop].assumption = D.tdve['alive'].ts[pop].vals[0]
            D.tdve['sus'].ts[pop].assumption = 0
        else:
            D.tdve['sus_all'].ts[pop].assumption = D.tdve['alive'].ts[pop].vals[0] - D.tdve['chronic'].ts[pop].vals[0] - D.tdve['acute'].ts[pop].vals[0]
            D.tdve['sus'].ts[pop].assumption = D.tdve['sus_all'].ts[pop].assumption*0.005
        
        # initialize undiagnosed amongst acute infections
        D.tdve['undiag_acute'].ts[pop].assumption = D.tdve['prop_undiag'].ts[pop].vals[0]*D.tdve['acute'].ts[pop].vals[0]
        
        pair_pars = [
            ['undiag_',['prop_undiag']],
            ['undiag_ltc_',['prop_undiag','prop_ud_ltc']],
            ['fail1_',['prop_tx_fail',2/3]],
            ['fail2_',['prop_tx_fail',1/3]],
            ['ltfu_ab_',['prop_ltfu_ab']],
            ['ltfu_pcr_',['prop_ltfu_pcr']],
            ['ltfu_tx_',['prop_ltfu_tx']],
            ]
        
        # initialize undiagnosed, treatment failures, LTFU after antibody, rna tests and treatments amonst the chronic
        for j in stages: # initialize undiagnosed, treatment failures, LTFU after antibody/rna tests and treatments
            if j == 'f0' or j == 'f1':
                coeff = D.tdve['prop_f0f2'].ts[pop].assumption * 4/9
            elif j == 'f2':
                coeff = D.tdve['prop_f0f2'].ts[pop].assumption * 1/9
            else:
                coeff = D.tdve['prop_' + j].ts[pop].assumption
                
            for par in pair_pars:
                multiplier = 1
                for par1 in par[1]:
                    if isinstance(par1, str):
                        if D.tdve[par1].ts[pop].vals:
                            multiplier *= D.tdve[par1].ts[pop].vals[0]
                        else:
                            multiplier *= D.tdve[par1].ts[pop].assumption
                    elif isinstance(par1, float):
                        multiplier *= par1
                D.tdve[par[0] + j].ts[pop].t = [2000]
                D.tdve[par[0] + j].ts[pop].vals = [multiplier*D.tdve['chronic'].ts[pop].vals[0]*coeff]
            D.tdve['tx1_' + j].ts[pop].t = [2000]
            D.tdve['tx1_' + j].ts[pop].vals = [0]
            D.tdve['tx2_' + j].ts[pop].t = [2000]
            D.tdve['tx2_' + j].ts[pop].vals = [0]
            
    for pop_pwid, pop_gp in zip(['PWID_males','PWID_females'],['18-64_males','18-64_females']):
        D.tdve['death_all'].ts[pop_pwid].vals = D.tdve['death_all'].ts[pop_gp].vals
    
    ## Uncertainty
    data = pd.read_excel(str(rootdir) + '/data/flat_datasheet.xlsx', sheet_name='Country-Standard Deviations')
    data = data[data.ISO3 == country]
    for col_name, col_data in data.items():
        if '__' not in col_name:
            continue
        col_name_split = col_name.split('__')
        par = col_name_split[0]
        pop = col_name_split[1]
        if pop == 'PWID_all':
            pops = ['PWID_males','PWID_females']
        elif pop == 'all':
            pops = D.pops
        elif pop == 'females_all':
            pops = [p for p in D.pops if '_females' in p]
        elif pop == 'males_all':
            pops = [p for p in D.pops if '_males' in p]
        else:
            pops = [pop]
        for pop in pops:
            if par == 'treat_overall':
                if pop in prop:
                    D.tdve[par].ts[pop].sigma = col_data.values[0]*prop[pop]*0.8
                else:
                    D.tdve[par].ts[pop].sigma = 0
            else:
                D.tdve[par].ts[pop].sigma = col_data.values[0]
                
    if savedir == None:
        savedir = str(rootdir) + f'/databooks/databook_{country}.xlsx'
    D.save(savedir)