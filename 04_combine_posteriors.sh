#!/bin/bash

#SBATCH --job-name=combinePosteriors
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --output=error_combinePosteriors.%j

RUN_LOCATION="/scratch/user/sblain/hummer_aims/run_hmm_autosomes"
HMM_INPUT="ancestryHMM_counts_HummerHybrids_f50_autosomes"
RAW_AIMS_MOD="HummerRefPanel.f50.noUnderscores.aims.mod"
run_tag="HummerHybrids"
samplesFile="ancestryHMM_samples_HummerHybrids_f50"

module load GCC/11.2.0  OpenMPI/4.1.1 R/4.2.0


cd $RUN_LOCATION

awk -F"\t" '{ print $1 }' $samplesFile > hmm_inds_"$run_tag"
awk -F"\t" '{print $1,$2}' $HMM_INPUT | sed "s/ /_/" > hmm_aims_"$run_tag".mod_1

#grep lines that match pattern
#w=word matching; f=get pattern from file; F=compare strings (not regular expressions)
grep -wFf hmm_aims_"$run_tag".mod_1 $RAW_AIMS_MOD | sed "s/_/\t/" | awk -F "\t" '{print $1,$2,$3,$4}' | sed "s/ /\t/g" > hmm_aims_"$run_tag"


AIMS_LIST=hmm_aims_"$run_tag"

cp $AIMS_LIST "$run_tag"_posterior_probs_tmp
sed  -i '1i CHROM\tPOS\tREF\tALT' "$run_tag"_posterior_probs_tmp

while read IND

        if [ -z $IND ] #avoid infinite loop if variable is empty
        then
                echo "individual variable empty"
                break
        fi

        do
                awk -F"\t" '{ print $3,$4,$5 }' "$IND".posterior | sed '1d' | sed 's/ /,/g' > "$IND".posterior_only
                sed  -i "1i $IND" "$IND".posterior_only #add individual name as header
        done < hmm_inds_"$run_tag"

paste -d"\t" "$run_tag"_posterior_probs_tmp *.posterior_only > "$run_tag"_posterior_probs

Rscript --vanilla get_posterior_states.R $run_tag


cp $AIMS_LIST "$run_tag"_tmp
sed  -i '1i CHROM\tPOS\tREF\tALT' "$run_tag"_tmp

paste -d"\t" "$run_tag"_tmp *.posterior_state > "$run_tag"_posterior_states
rm *.posterior_state *.posterior_only
rm *tmp