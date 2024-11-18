# Summarize alignments 

# ODIR should be the same output directory as defined in map_barcode_reads.sh
ODIR=$(date +'../output/MappingOutput-%F')

ls -l ${ODIR}/*.map > /dev/null 2>&1    

# Check if .map files exist
if [ "$?" = "0" ]; then
   cd ${ODIR}                          
   # runs a Java program named summarizeAlignments.jar with a maximum heap size of 4096 MB and two arguments (reference txt and .map file)
   java -Xmx4096M -jar ../../Library/summarizeAlignments.jar ../../Library/Candida_barcode_V16_final_annotation_GRACEV1V2.txt `find ./ -name '*.map' -print | tr '\n' ','`         

   if [ "$?" > "0" ]; then
      echo "Program completes, or there could be an error merging files."
      echo "Please check if all alignments are done and the columns & reads in the merged_results.txt file are complete."
      cd -
      exit 1
   fi

   cd -
else
   echo "No files to merge."
fi
