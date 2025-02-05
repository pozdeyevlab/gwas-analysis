# Read in named command line args 
while getopts ":i:w:d:p:o:" opt; do
  case $opt in
    i) input_pgen="$OPTARG"
    ;;
    w) weight_file="$OPTARG"
    ;;
    d) phenotype="$OPTARG"
    ;;
    p) plink="$OPTARG"
    ;;
    o) output_dir="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

# Make output directory if it does not already exist 
$output_dir/$phenotype

#Create temp var list formatted for pgens
cut -f1 $weight_file > $output_dir/$phenotype/$phenotype.varlist
chroms=$(cut -d':' -f1 $output_dir/$phenotype/$phenotype.varlist |grep -v 'X' |sort |uniq)

echo $chroms

sed -i.bak 's/^/chr/g' $output_dir/$phenotype/$phenotype.varlist

# For every chrom 
for CHR in $chroms;
do
  echo ${CHR}
    # Construct the file prefix for the current chromosome
    input_file="${input_pgen}${CHR}"
    echo $input_file 

    # Construct the output file name
    output_file="${output_dir}/$phenotype/${CHR}"

    # Extract variants from files
    ${plink} --pfile ${input_file} \
      --extract $output_dir/$phenotype/$phenotype.varlist \
      --make-pgen \
      --out ${output_file}
    
    # Remove plink log for cleaner directory
    rm ${output_file}.log

    # Check if the command was successful
    if [ $? -ne 0 ]; then
      echo "Error processing chromosome ${CHR}"
      exit 1
    else
      echo "Successfully processed chromosome ${CHR}"
    fi
done

# Merge variant specific pgens
ls $output_dir/$phenotype/*.psam | sed 's/.psam//g' > $output_dir/$phenotype/merge.tsv

${plink} \
  --pmerge-list $output_dir/$phenotype/merge.tsv \
  --export vcf \
  --out ${output_dir}/${phenotype}/${phenotype}

# Make genotype header from psam
echo 'ID' > ${output_dir}/${phenotype}/sample_ids.tsv
tail -n +2 ${output_dir}/${phenotype}/${phenotype}.psam | cut -f1 >> ${output_dir}/${phenotype}/sample_ids.tsv

# Remove temp pgens
rm $output_dir/$phenotype/*.pgen
rm $output_dir/$phenotype/*.pvar
rm $output_dir/$phenotype/*.psam

# Make genotype data
vcf="${output_dir}/${phenotype}/${phenotype}.vcf"

# Get variant id and genotypes
bcftools query -f '%ID[\t%GT]\n' $vcf > ${output_dir}/${phenotype}/gt_data.tsv

count1=`wc -l < ${output_dir}/${phenotype}/gt_data.tsv`

# Use sed to enforce pipe separated genotypes (change \ and / to |)
sed -i.bak 's/\//|/g' ${output_dir}/${phenotype}/gt_data.tsv
sed -i.bak 's/\\/|/g' ${output_dir}/${phenotype}/gt_data.tsv

# Covert GT's with sed -- faster than using logic in python
sed -i.bak 's/1|0/1/g' ${output_dir}/${phenotype}/gt_data.tsv
sed -i.bak 's/0|1/1/g' ${output_dir}/${phenotype}/gt_data.tsv
sed -i.bak 's/1|1/2/g' ${output_dir}/${phenotype}/gt_data.tsv
sed -i.bak 's/0|0/0/g' ${output_dir}/${phenotype}/gt_data.tsv
sed -i.bak 's/.|./NA/g' ${output_dir}/${phenotype}/gt_data.tsv

# Delete intermediate files
rm ${output_dir}/${phenotype}/*.bak
$rm ${vcf}

count1=`wc -l < ${output_dir}/${phenotype}/gt_data.tsv`

if [[ $count1 = 0 ]]; then
    echo "****************************"
    echo "**** ERROR - FILE EMPTY ****"
    echo "****************************"
    exit 0
fi
