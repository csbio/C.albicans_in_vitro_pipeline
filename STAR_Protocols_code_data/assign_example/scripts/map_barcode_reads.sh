# Generate output directory 

ODIR=$(date +'../output/MappingOutput-%F')
if [ ! -d ${ODIR} ]
then
   mkdir -p ${ODIR}
fi

# --

# Annotate reads to barcodes

for sample in $( find ./ -name '*.fastq.gz' ! -name '._*.fastq.gz' -print | xargs -I {} basename {} | sed -E 's/_R1_001.fastq.gz//' | sort | uniq )
do

   oFile=$(echo $sample".map");
   unmapped=$(echo $sample".unmapped");

   echo $sample"->"$oFile;
   fList=`find ./ -name $sample'*.fastq.gz' -print | tr '\n' ' ' | sed -e 's/,$//g'`

   echo "File List:  "$fList >> ${ODIR}/mapping.log 2>&1

   # Run bowtie, using PERL script to find sequence-of-interest in reads
   # zcat is a common command for decompressing .gz files on most Unix and Linux systems
   # macOS users can use either gzcat or zcat if both available
   zcat $fList | bowtie -p 3 -t all_barcodes_v3 -l 20 -m 1 -n 1 --best --strata --un ${ODIR}/$unmapped - ${ODIR}/$oFile >> ${ODIR}/mapping.log 2>&1

   ###
   ## Add in error check for bowtie
   ###

done
