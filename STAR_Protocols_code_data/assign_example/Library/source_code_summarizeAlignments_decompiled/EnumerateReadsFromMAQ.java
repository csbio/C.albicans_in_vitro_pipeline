import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Iterator;

public class EnumerateReadsFromMAQ {
   private String inFile = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Run3_June15_2010/s_1_mapped.txt";
   private String outFile = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Run3_June15_2010/s_1_mapped_results.txt";
   private String dir = "/Volumes/Bioinf_Data1/SeqData/Milyavsky080311/L007/L7";
   private Hashtable results = new Hashtable();

   public EnumerateReadsFromMAQ() {
      this.readResults(this.inFile);
      this.writeResults(this.outFile);
   }

   public void writeResults(String outFile) {
      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
         Iterator it = this.results.keySet().iterator();

         while(it.hasNext()) {
            String key = (String)it.next();
            EnumerateReadsFromMAQ.Read num = (EnumerateReadsFromMAQ.Read)this.results.get(key);
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
            if (!this.results.containsKey(data[2])) {
               this.results.put(data[2], new EnumerateReadsFromMAQ.Read());
            } else {
               ((EnumerateReadsFromMAQ.Read)this.results.get(data[2])).incrCount();
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
      new EnumerateReadsFromMAQ();
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
