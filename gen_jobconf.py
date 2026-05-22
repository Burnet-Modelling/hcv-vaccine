from hcv import utils as ut

countries = ut.country_scope()

with open("job.conf","w") as f:
    for i,c in enumerate(countries):
        f.write(f"{i} python srun_main.py {c}\n")