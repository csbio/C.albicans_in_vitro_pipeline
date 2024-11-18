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

public class ExtractBowtieAlignments_V1 {
   private Hashtable idxToTRCID = null;
   private String idxFile = "/Volumes/Bioinf_Data1/SeqData/Ketela_Martin_Sept2011/mouse_genome.idx";
   private String hpAnnotations = "/Users/kbrown/Documents/Work/Scoring of shRNA/Chip Annotation/AllTRC_in_TRPS.txt";
   private String dir = "/Volumes/Bioinf_Data1/SeqData/Ketela_Fine_Sept2011/";
   private Hashtable annotations = null;
   private Hashtable hairpins = new Hashtable();
   private Hashtable[] allResults = null;

   public ExtractBowtieAlignments_V1() {
      this.idxToTRCID = this.readHairpinIndex(this.idxFile);
      this.annotations = this.readAnnotations(this.hpAnnotations);
      System.out.println("Read in " + this.idxToTRCID.size() + " indices...");
      Vector v = this.getFilenames(this.dir);
      System.out.println(v.size() + " files to parse...");
      this.allResults = new Hashtable[v.size()];
      int fNum = 0;

      String outFile;
      Hashtable results;
      for(Iterator it = v.iterator(); it.hasNext(); this.allResults[fNum++] = this.writeResults(results, outFile)) {
         File f = (File)it.next();
         outFile = f.getAbsolutePath().replace(".map", ".counts");
         String outFile2 = f.getAbsolutePath().replace(".map", ".ambig");
         results = this.readResults(f.getAbsolutePath(), outFile2);
      }

      this.printMatrix(this.allResults, this.dir + "/merged_results.txt", v);
   }

   public Hashtable readAnnotations(String inFile) {
      Hashtable h = new Hashtable();

      try {
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));
         String line = null;
         int var5 = 0;

         while((line = in.readLine()) != null) {
            String[] data = line.split("\t");
            if (var5++ == 0) {
               h.put("header", line.substring(line.indexOf("\t", 0) + 1, line.length()));
            } else {
               h.put(data[0], line.substring(line.indexOf("\t", 0) + 1, line.length()));
            }
         }

         in.close();
      } catch (FileNotFoundException var7) {
         var7.printStackTrace();
      } catch (IOException var8) {
         var8.printStackTrace();
      }

      return h;
   }

   public void printMatrix(Hashtable[] results, String outFile, Vector files) {
      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
         out.write("\t" + (String)this.annotations.get("header"));
         Iterator it = files.iterator();

         while(it.hasNext()) {
            File f = (File)it.next();
            String fName = f.getCanonicalPath();
            fName = fName.substring(fName.lastIndexOf("/") + 1, fName.lastIndexOf("."));
            out.write("\t" + fName);
         }

         out.newLine();
         it = this.hairpins.keySet().iterator();

         while(it.hasNext()) {
            String hp = (String)it.next();
            out.write(hp + "\t" + (String)this.annotations.get(hp));

            for(int i = 0; i < files.size(); ++i) {
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

   public Vector getFilenames(String dir) {
      Vector v = new Vector();
      File f = new File(dir);
      File[] files = f.listFiles();

      for(int i = 0; i < files.length; ++i) {
         if (files[i].getName().matches(".+\\.map$")) {
            v.add(files[i].getAbsoluteFile());
            System.out.println(files[i].getAbsolutePath());
         }
      }

      return v;
   }

   public Hashtable writeResults(Hashtable results, String outFile) {
      Hashtable h = new Hashtable();
      int failedMatches = 0;
      int trashedReads = 0;
      System.out.println("Writing to " + outFile);

      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
         Iterator it = results.keySet().iterator();

         while(true) {
            while(it.hasNext()) {
               Integer key = (Integer)it.next();
               ExtractBowtieAlignments_V1.Read num = (ExtractBowtieAlignments_V1.Read)results.get(key);
               Integer keyInt = new Integer(key);
               ExtractBowtieAlignments_V1.Sequence s = (ExtractBowtieAlignments_V1.Sequence)this.idxToTRCID.get(keyInt);
               if (s == null) {
                  if (this.idxToTRCID.containsKey(keyInt - 1)) {
                     s = (ExtractBowtieAlignments_V1.Sequence)this.idxToTRCID.get(keyInt - 1);
                  } else {
                     if (!this.idxToTRCID.containsKey(keyInt + 1)) {
                        System.err.println("Couldn't find idx " + keyInt + " - Reads:  " + num.getCount());
                        ++failedMatches;
                        trashedReads += num.getCount();
                        continue;
                     }

                     s = (ExtractBowtieAlignments_V1.Sequence)this.idxToTRCID.get(keyInt + 1);
                  }
               }

               Iterator it2 = s.clone_id.iterator();

               while(it2.hasNext()) {
                  String clone = (String)it2.next();
                  out.write(clone + "\t" + num.getCount());
                  out.newLine();
                  h.put(clone, new Integer(num.getCount()));
               }
            }

            out.close();
            break;
         }
      } catch (IOException var14) {
         var14.printStackTrace();
      }

      System.err.println(failedMatches + " aligned transcripts failed to match idx.");
      System.err.println("Number of reads trashed due to mismatched alignment idx:  " + trashedReads);
      return h;
   }

   public Hashtable readResults(String inFile, String outFile) {
      System.out.println("Reading from " + inFile);
      int numReads = 0;
      Hashtable reads = new Hashtable();
      Hashtable results = new Hashtable();

      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));

         for(String line = null; (line = in.readLine()) != null; ++numReads) {
            String[] data = line.split("\t");
            if (!this.idxToTRCID.containsKey(new Integer(data[3]))) {
               out.write(line);
               out.newLine();
            }

            if (!results.containsKey(new Integer(data[3]))) {
               results.put(new Integer(data[3]), new ExtractBowtieAlignments_V1.Read());
            } else {
               ((ExtractBowtieAlignments_V1.Read)results.get(new Integer(data[3]))).incrCount();
            }
         }

         in.close();
         out.close();
      } catch (FileNotFoundException var10) {
         var10.printStackTrace();
      } catch (IOException var11) {
         var11.printStackTrace();
      }

      System.out.println("Read in " + numReads + " total alignments from " + reads.size() + " reads...");
      System.out.println("Includes counts on " + results.size() + " hairpins!");
      return results;
   }

   private Hashtable readHairpinIndex(String inFile) {
      Hashtable h = new Hashtable();

      try {
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));
         String line = in.readLine();

         while((line = in.readLine()) != null) {
            String[] data = line.split("\t");
            Integer idx = new Integer(data[0]);
            ExtractBowtieAlignments_V1.Sequence s = new ExtractBowtieAlignments_V1.Sequence(data[1], data[2]);
            h.put(idx, s);
            String[] hps = data[2].split(";");

            for(int i = 0; i < hps.length; ++i) {
               this.hairpins.put(hps[i], new Integer(this.hairpins.size() + 1));
            }
         }

         in.close();
      } catch (IOException var10) {
         var10.printStackTrace();
      }

      System.out.println("Read in " + this.hairpins.size() + " total hairpins...");
      return h;
   }

   public static void main(String[] args) {
      new ExtractBowtieAlignments_V1();
      int t = 2359303;
      System.out.println("Test #:  " + t);
      t = Math.round((float)(t / 10)) * 10;
      System.out.println("Rounded #:  " + t);
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
         return this.subpool == ((ExtractBowtieAlignments_V1.Sequence)s2).getSubpool();
      }
   }
}
