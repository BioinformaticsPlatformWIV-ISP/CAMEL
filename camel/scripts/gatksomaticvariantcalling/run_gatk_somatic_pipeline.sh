#! /bin/bash
echo 'running GATK somatic variant calling pipeline'
echo $@
# for testing
#. /home/todel/PycharmProjects/camel_3.0/camel_3_env/bin/activate
#export PYTHONPATH=/home/todel/PycharmProjects/camel_3.0/
#python3 /home/todel/PycharmProjects/camel_3.0/camel/scripts/gatksomaticvariantcalling/run_gatk_somatic_pipeline.py "$@"  --wdir $PWD

# for real running with lmod
python3 /usr/local/bin/BIOIT/GATK_somatic_variant_calling/camel_3.0/camel/scripts/gatksomaticvariantcalling/run_gatk_somatic_pipeline.py "$@"
