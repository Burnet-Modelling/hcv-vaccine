from hcv import utils as ut

countries = ut.country_scope()
# countries = ['CHN','EGY','ZAF','UKR','BRA','MEX']

with open("job.conf","w") as f:
    for i,c in enumerate(countries):
        f.write(f"{i} python srun_main.py {c}\n")