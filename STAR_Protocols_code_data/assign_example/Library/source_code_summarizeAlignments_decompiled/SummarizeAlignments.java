import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Iterator;

public class SummarizeAlignments {
   private Hashtable annotations;
   private Hashtable[] allResults = null;

   public SummarizeAlignments(String[] args) {
      if (args.length < 2) {
         usage();
      }

      this.readAnnotations(args[0]);
      System.out.println("Read in annotations for " + this.annotations.size() + " hairpins...");
      String[] fileNames = args[1].split(",");
      System.out.println(fileNames.length + " files to parse...");
      this.allResults = new Hashtable[fileNames.length];
      int fNum = 0;

      for(int i = 0; i < fileNames.length; ++i) {
         File f = new File(fileNames[i]);
         String outFile = f.getAbsolutePath().replace(".map", ".counts");
         String outFile2 = f.getAbsolutePath().replace(".map", ".ambig");
         Hashtable results = this.readResults(f.getAbsolutePath(), outFile2);
         this.allResults[fNum++] = this.writeResults(results, outFile);
      }

      this.printMatrix(this.allResults, "merged_results.txt", fileNames);
   }

   public void printMatrix(Hashtable[] results, String outFile, String[] files) {
      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
         out.write("TRC.ID\t" + (String)this.annotations.get("header"));

         for(int i = 0; i < files.length; ++i) {
            File f = new File(files[i]);
            String fName = f.getCanonicalPath();
            fName = fName.substring(fName.lastIndexOf("/") + 1, fName.lastIndexOf("."));
            out.write("\t" + fName);
         }

         out.newLine();
         this.annotations.remove("header");
         Iterator it = this.annotations.keySet().iterator();

         while(it.hasNext()) {
            String hp = (String)it.next();
            out.write(hp + "\t" + (String)this.annotations.get(hp));

            for(int i = 0; i < files.length; ++i) {
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

   public Hashtable writeResults(Hashtable results, String outFile) {
      throw new Error("Unresolved compilation problems: \n\tThe method getCount() is undefined for the type SummarizeAlignments.Read\n\tThe method getCount() is undefined for the type SummarizeAlignments.Read\n\tThe method getCount() is undefined for the type SummarizeAlignments.Read\n");
   }

   public Hashtable readResults(String inFile, String outFile) {
      throw new Error("Unresolved compilation problem: \n\tThe method incrCount() is undefined for the type SummarizeAlignments.Read\n");
   }

   private void readAnnotations(String inFile) {
      this.annotations = new Hashtable();

      try {
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));
         String line = null;
         int lNum = 0;

         while((line = in.readLine()) != null) {
            System.out.println("Line " + lNum + ": " + line);
            String[] data = line.split("\t");
            if (lNum++ == 0) {
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
      System.out.println("Usage:  SummarizeAlignments <annotationsFile> <map1.map,...,mapN.map>");
      System.out.println("        annotationsFile:  TRC_ID assumed to be in first column");
      System.out.println("                          First row is header row");
      System.out.println("");
      System.out.println("        map1.map - Bowtie output files.  ID (col3) assumed to be TRC_ID, column ");
      System.out.println("");
      System.exit(0);
   }

   public static void main(String[] args) {
      new SummarizeAlignments(args);
   }
}
