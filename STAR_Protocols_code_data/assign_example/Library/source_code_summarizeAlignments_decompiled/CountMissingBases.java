import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class CountMissingBases {
   private String inFile = "/Users/kbrown/Documents/Work/Scoring of shRNA/Franco's Screen/ABI SOLID Seq - Dec 2010/solid0085_20101022_JM28plex_FC2/bcSample1/results.F1B1/libraries/franco28/primary.20101114062155810/reads/reads.fasta";
   private String outFile = "/Users/kbrown/Documents/Work/Scoring of shRNA/Franco's Screen/ABI SOLID Seq - Dec 2010/solid0085_20101022_JM28plex_FC2/bcSample1/results.F1B1/libraries/franco28/primary.20101114062155810/reads/read_base_freqs.out";
   private int[] counts = new int[51];

   public CountMissingBases() {
      int numSeqs = 0;
      String line = null;

      int i;
      try {
         BufferedReader in;
         for(in = new BufferedReader(new FileReader(new File(this.inFile))); (line = in.readLine()) != null; ++numSeqs) {
            for(i = 0; i < line.length(); ++i) {
               if (line.charAt(i) != '.') {
                  int var10002 = this.counts[i]++;
               }
            }
         }

         in.close();
      } catch (FileNotFoundException var6) {
         var6.printStackTrace();
      } catch (IOException var7) {
         var7.printStackTrace();
      }

      System.out.println("Processed " + numSeqs + " sequences...");

      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(this.outFile)));

         for(i = 0; i < this.counts.length; ++i) {
            out.write(Double.toString((double)this.counts[i] / (double)numSeqs));
            if (i + 1 < this.counts.length) {
               out.write("\t");
            }
         }

         out.newLine();
         out.flush();
         out.close();
      } catch (IOException var5) {
         var5.printStackTrace();
      }

   }

   public static void main(String[] args) {
      new CountMissingBases();
   }
}
