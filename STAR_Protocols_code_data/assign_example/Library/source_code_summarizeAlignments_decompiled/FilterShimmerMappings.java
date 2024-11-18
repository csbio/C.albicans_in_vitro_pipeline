import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class FilterShimmerMappings {
   private String filterFile = "/Volumes/Bioinf_Data1/SeqData/Project_Ketela_Oct2011/shimmer_hairpin_indexes.txt";
   private String readsFile = "/Volumes/Bioinf_Data1/SeqData/Project_Ketela_Oct2011/Sample_Shimmer_1_ACAA/test.out";
   private String outFile1 = "/Volumes/Bioinf_Data1/SeqData/Project_Ketela_Oct2011/Sample_Shimmer_1_ACAA/inPool.out";
   private String outFile2 = "/Volumes/Bioinf_Data1/SeqData/Project_Ketela_Oct2011/Sample_Shimmer_1_ACAA/outOfPool.out";
   private String outFileQ1 = "/Volumes/Bioinf_Data1/SeqData/Project_Ketela_Oct2011/Sample_Shimmer_1_ACAA/inPool.quals";
   private String outFileQ2 = "/Volumes/Bioinf_Data1/SeqData/Project_Ketela_Oct2011/Sample_Shimmer_1_ACAA/outOfPool.quals";

   public FilterShimmerMappings() {
      Hashtable filter = this.readFilterFile(this.filterFile);
      this.filterReads(this.readsFile, this.outFile1, this.outFile2, this.outFileQ1, this.outFileQ2, filter);
   }

   private void filterReads(String inFile, String outFile1, String outFile2, String outFileQ1, String outFileQ2, Hashtable filter) {
      try {
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));
         String line = null;
         BufferedWriter out1 = new BufferedWriter(new FileWriter(new File(outFile1)));
         BufferedWriter out2 = new BufferedWriter(new FileWriter(new File(outFile2)));
         BufferedWriter outQ1 = new BufferedWriter(new FileWriter(new File(outFileQ1)));
         BufferedWriter outQ2 = new BufferedWriter(new FileWriter(new File(outFileQ2)));

         while((line = in.readLine()) != null) {
            String[] data = line.split("\t");
            if (filter.containsKey(new Integer(data[3]))) {
               out1.write(line);
               out1.newLine();
               outQ1.write(convertQualValues(data[5]));
               outQ1.newLine();
            } else {
               out2.write(line);
               out2.newLine();
               outQ2.write(convertQualValues(data[5]));
               outQ2.newLine();
            }
         }

         in.close();
         out1.flush();
         out1.close();
         outQ1.flush();
         outQ1.close();
         out2.flush();
         out2.close();
         outQ2.flush();
         outQ2.close();
      } catch (FileNotFoundException var14) {
         var14.printStackTrace();
      } catch (IOException var15) {
         var15.printStackTrace();
      }

   }

   private static String convertQualValues(String qualStr) {
      String retStr = "";

      for(int i = 0; i < qualStr.length(); ++i) {
         retStr = retStr + (qualStr.charAt(i) - 33);
         if (i + 1 < qualStr.length()) {
            retStr = retStr + "\t";
         }
      }

      return retStr;
   }

   private Hashtable readFilterFile(String inFile) {
      Hashtable h = new Hashtable();

      try {
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));
         String line = null;

         while((line = in.readLine()) != null) {
            String[] data = line.split("\t");
            h.put(new Integer(data[0]), data);
         }

         in.close();
      } catch (NumberFormatException var6) {
         var6.printStackTrace();
      } catch (FileNotFoundException var7) {
         var7.printStackTrace();
      } catch (IOException var8) {
         var8.printStackTrace();
      }

      return h;
   }

   public static void main(String[] args) {
      new FilterShimmerMappings();
   }
}
