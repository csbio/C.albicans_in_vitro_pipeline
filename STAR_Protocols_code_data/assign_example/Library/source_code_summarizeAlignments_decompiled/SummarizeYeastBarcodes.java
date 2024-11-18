import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

public class SummarizeYeastBarcodes {
   private Hashtable annotations;
   private Hashtable[] allResults = null;

   public SummarizeYeastBarcodes(String[] args) {
      if (args.length < 2) {
         usage();
      }

      this.readAnnotations(args[0]);
      System.out.println("Read in annotations for " + this.annotations.size() + " hairpins...");
      String[] fileNames = args[1].split(",");
      System.out.println(fileNames.length + " files to parse...");
      this.allResults = new Hashtable[fileNames.length * 2];
      int fNum = 0;

      for(int i = 0; i < fileNames.length; ++i) {
         File f = new File(fileNames[i]);
         String outFile = f.getAbsolutePath().replace(".map", ".counts");
         String outFile2 = f.getAbsolutePath().replace(".map", ".ambig");
         Hashtable results = this.readResults(f.getAbsolutePath(), outFile2);
         Hashtable[] resultsHash = this.writeResults(results, outFile);
         this.allResults[fNum++] = resultsHash[0];
         this.allResults[fNum++] = resultsHash[1];
      }

      this.printMatrix(this.allResults, "merged_results.txt", fileNames);
   }

   public void printMatrix(Hashtable[] results, String outFile, String[] files) {
      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
         out.write("Strain.ID\t" + (String)this.annotations.get("header"));

         for(int i = 0; i < files.length; ++i) {
            File f = new File(files[i]);
            String fName = f.getCanonicalPath();
            fName = fName.substring(fName.lastIndexOf("/") + 1, fName.lastIndexOf("."));
            out.write("\t" + fName + "_UP\t" + fName + "_DOWN");
         }

         out.newLine();
         this.annotations.remove("header");
         Iterator it = this.annotations.keySet().iterator();

         while(it.hasNext()) {
            String hp = (String)it.next();
            out.write(hp + "\t" + (String)this.annotations.get(hp));

            for(int i = 0; i < results.length; ++i) {
               Hashtable result = results[i];
               if (result.containsKey(hp)) {
                  out.write("\t" + (Integer)result.get(hp));
               } else {
                  out.write("\t0");
               }
            }

            out.newLine();
         }

         out.flush();
         out.close();
      } catch (IOException var9) {
         var9.printStackTrace();
      }

   }

   public Hashtable[] writeResults(Hashtable results, String outFile) {
      Hashtable h_up = new Hashtable();
      Hashtable h_down = new Hashtable();
      int failedMatches = false;
      int trashedReads = false;
      System.out.println("Writing to " + outFile);

      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
         Iterator it = results.keySet().iterator();

         while(it.hasNext()) {
            String trcid = (String)it.next();
            SummarizeYeastBarcodes.Read num = (SummarizeYeastBarcodes.Read)results.get(trcid);
            String[] trcids = trcid.split(";");

            for(int i = 0; i < trcids.length; ++i) {
               String clone = trcids[i];
               out.write(clone + "\t" + num.toString());
               out.newLine();
               if (!h_up.containsKey(clone)) {
                  h_up.put(clone, new Integer(num.getUpCount()));
               } else {
                  h_up.put(clone, new Integer((Integer)h_up.get(clone) + num.getUpCount()));
               }

               if (!h_down.containsKey(clone)) {
                  h_down.put(clone, new Integer(num.getDownCount()));
               } else {
                  h_down.put(clone, new Integer((Integer)h_down.get(clone) + num.getDownCount()));
               }
            }
         }

         out.close();
      } catch (IOException var14) {
         var14.printStackTrace();
      }

      Hashtable[] h = new Hashtable[]{h_up, h_down};
      return h;
   }

   public Hashtable readResults(String inFile, String outFile) {
      System.out.println("Reading from " + inFile);
      int numReads = 0;
      new Hashtable();
      Hashtable results = new Hashtable();

      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));

         for(String line = null; (line = in.readLine()) != null; ++numReads) {
            String[] data = line.split("\t");
            String primaryKey = data[2].replaceAll("_(UP|DOWN)", "");
            if (!this.annotations.containsKey(primaryKey)) {
               out.write(line);
               out.newLine();
            }

            if (!results.containsKey(primaryKey)) {
               results.put(primaryKey, new SummarizeYeastBarcodes.Read(data[2]));
            } else {
               ((SummarizeYeastBarcodes.Read)results.get(primaryKey)).incrementCount(data[2]);
            }
         }

         in.close();
         out.close();
      } catch (FileNotFoundException var11) {
         var11.printStackTrace();
      } catch (IOException var12) {
         var12.printStackTrace();
      }

      System.out.println("Read in " + numReads + " representing " + results.size() + " Strain.IDs...");
      return results;
   }

   private void readAnnotations(String inFile) {
      this.annotations = new Hashtable();

      try {
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));
         String line = null;
         int var4 = 0;

         while((line = in.readLine()) != null) {
            String[] data = line.split("\t");
            if (var4++ == 0) {
               this.annotations.put("header", line.substring(line.indexOf("\t", 0) + 1, line.length()));
            } else {
               this.annotations.put(data[0], line.substring(line.indexOf("\t", 0) + 1, line.length()));
            }
         }

         in.close();
      } catch (FileNotFoundException var6) {
         var6.printStackTrace();
      } catch (IOException var7) {
         var7.printStackTrace();
      }

   }

   private static void usage() {
      System.out.println("");
      System.out.println("Usage:  SummarizeYeastBarcodes <annotationsFile> <map1.map,...,mapN.map>");
      System.out.println("        annotationsFile:  TRC_ID assumed to be in first column");
      System.out.println("                          First row is header row");
      System.out.println("");
      System.out.println("        map1.map - Bowtie output files.  ID (col3) assumed to be TRC_ID, column ");
      System.out.println("");
      System.exit(0);
   }

   public static void main(String[] args) {
      new SummarizeYeastBarcodes(args);
   }

   class Read {
      private int upCount = 0;
      private int downCount = 0;

      public Read(String bcID) {
         this.incrementCount(bcID);
      }

      public void incrementCount(String bcID) {
         if (bcID.contains("_UP")) {
            this.incrUpCount();
         } else {
            this.incrDownCount();
         }

      }

      public int getUpCount() {
         return this.upCount;
      }

      public int getDownCount() {
         return this.downCount;
      }

      public void incrUpCount() {
         ++this.upCount;
      }

      public void incrDownCount() {
         ++this.downCount;
      }

      public String toString() {
         return this.upCount + "\t" + this.downCount;
      }
   }

   class Sequence {
      private String seq;
      private Vector clone_id;
      private int subpool;

      public Sequence(String seq, String clone_id) {
         this.seq = seq;
         this.clone_id = new Vector();
         this.addCloneIDs(clone_id);
      }

      private void addCloneIDs(String clone_id) {
         String[] data = clone_id.split(";");

         for(int i = 0; i < data.length; ++i) {
            this.clone_id.add(data[i]);
         }

      }

      public String getSeq() {
         return this.seq;
      }

      public void addClone_ID(String clone_id) {
         if (!this.clone_id.contains(clone_id)) {
            this.clone_id.add(clone_id);
         }

      }

      public String getClone_id() {
         return this.getCloneIDString();
      }

      private String getCloneIDString() {
         String str = "";
         Iterator it = this.clone_id.iterator();

         while(it.hasNext()) {
            str = str.concat((String)it.next());
            if (it.hasNext()) {
               str = str + ";";
            }
         }

         return str;
      }

      public int getSubpool() {
         return this.subpool;
      }

      public boolean equals(Object s2) {
         return this.subpool == ((SummarizeYeastBarcodes.Sequence)s2).getSubpool();
      }
   }
}
