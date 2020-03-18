Run on cluster with:
```
#with immediate submit
snakemake -j 100 --use-envmodules   --default-resources time=60 mem_mb=4000 --immediate-submit --notemp --cluster="./immediate_submit.py {dependencies}"

#without immediate submit
snakemake  --cluster "sbatch --parsable --account=rrg-akhanf --time={resources.time} --cpus-per-task={threads} --mem={resources.mem_mb}"  -j 4 --use-envmodules   --default-resources time=60 mem_mb=4000 
```
