from hcv import utils as ut

countries = ut.country_scope()

with open("job.conf","w") as f:
    for i,c in enumerate(countries):
        f.write(f"{i} python srun_main.py {c}\n")
        
n_samples = [100, 500, 1000]

with open("conv.conf","w") as f:
    for i,c in enumerate(n_samples):
        f.write(f"{i} python run_convergence.py {c}\n")