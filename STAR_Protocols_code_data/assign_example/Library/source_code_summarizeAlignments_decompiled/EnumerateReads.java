import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Iterator;

public class EnumerateReads {
   private String inFile = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Run1_June14_2010/blast_results_run1.txt";
   private String outFile = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Run1_June14_2010/blast_results_run1_results.txt";
   private Hashtable results = new Hashtable();

   public EnumerateReads() {
      this.readResults(this.inFile);
      this.writeResults(this.outFile);
   }

   public void writeResults(String outFile) {
      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
         Iterator it = this.results.keySet().iterator();

         while(it.hasNext()) {
            String key = (String)it.next();
            EnumerateReads.Read num = (EnumerateReads.Read)this.results.get(key);
            out.write(key + "\t" + num.getCount());
            out.newLine();
         }

         out.close();
      } catch (IOException var6) {
         var6.printStackTrace();
      }

   }

   public void readResults(String inFile) {
      int numReads = 0;
      new Hashtable();

      try {
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));

         for(String line = null; (line = in.readLine()) != null; ++numReads) {
            String[] data = line.split("\t");
            if (!this.results.containsKey(data[1])) {
               this.results.put(data[1], new EnumerateReads.Read());
            } else {
               ((EnumerateReads.Read)this.results.get(data[1])).incrCount();
            }
         }

         in.close();
      } catch (FileNotFoundException var7) {
         var7.printStackTrace();
      } catch (IOException var8) {
         var8.printStackTrace();
      }

   }

   public static void main(String[] args) {
      new EnumerateReads();
   }

   class Read {
      private int count = 1;

      public Read() {
      }

      public int getCount() {
         return this.count;
      }

      public void incrCount() {
         ++this.count;
      }
   }
}
