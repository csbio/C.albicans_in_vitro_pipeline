# Summarize alignments 

# ODIR should be the same output directory as defined in map_barcode_reads.sh
ODIR=$(date +'../output/MappingOutput-%F')

ls -l ${ODIR}/*.map > /dev/null 2>&1    
if [ "$?" = "0" ]; then
   cd ${ODIR}                          
   java -Xmx4096M -jar ../../Library/summarizeAlignments.jar ../../Library/Candida_barcode_V16_final_annotation_GRACEV1V2.txt `find ./ -name '*.map' -print | tr '\n' ','`         

   if [ "$?" > "0" ]; then
      echo "There was an error merging files!"
      cd -
      exit 1
   fi

   cd -
else
   echo "No files to merge."
fi
