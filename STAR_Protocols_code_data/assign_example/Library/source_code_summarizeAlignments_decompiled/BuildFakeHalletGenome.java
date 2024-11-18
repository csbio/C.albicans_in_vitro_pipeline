import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Random;
import java.util.Vector;

public class BuildFakeHalletGenome {
   private String seqFile = "/Volumes/Bioinf_Data1/SeqData/MikeHallet_July2011/Hallet pool contents.txt";
   private String genomeOutFile = "/Volumes/Bioinf_Data1/SeqData/MikeHallet_July2011/hallet_genome.fasta";
   private String idxOutFile = "/Volumes/Bioinf_Data1/SeqData/MikeHallet_July2011/hallet_genome.idx";
   private String bedOutFile = "/Volumes/Bioinf_Data1/SeqData/MikeHallet_July2011/hallet_genome.bed";
   private Hashtable clone2seq = new Hashtable();

   public BuildFakeHalletGenome() {
      Hashtable seqs = this.readFile(this.seqFile);
      Vector seqV = this.seqHashToVector(seqs);
      this.buildGenome(seqV);
   }

   private void buildFastaFile(Vector sequences) {
      int seqNum = 0;

      try {
         BufferedWriter outGenome = new BufferedWriter(new FileWriter(new File(this.genomeOutFile)));

         for(Iterator it = sequences.iterator(); it.hasNext(); ++seqNum) {
            BuildFakeHalletGenome.Sequence s = (BuildFakeHalletGenome.Sequence)it.next();
            outGenome.write(">trc_hairpin|" + s.getCloneIDString());
            outGenome.newLine();
            outGenome.write(s.getSeq());
            outGenome.newLine();
         }

         outGenome.flush();
         outGenome.close();
      } catch (IOException var6) {
         var6.printStackTrace();
      }

      System.out.println("Sequences written:  " + seqNum);
   }

   private void buildGenome(Vector sequences) {
      try {
         BufferedWriter outGenome = new BufferedWriter(new FileWriter(new File(this.genomeOutFile)));
         BufferedWriter outIdx = new BufferedWriter(new FileWriter(new File(this.idxOutFile)));
         BufferedWriter outBed = new BufferedWriter(new FileWriter(new File(this.bedOutFile)));
         outGenome.write(">gnl|TRC|hairpin_chromosome1");
         outGenome.newLine();

         for(int i = 0; i < 100; ++i) {
            outGenome.write("N");
         }

         outIdx.write("Genome Location\tSequence\tMatching TRC.ID(s)");
         outIdx.newLine();
         outBed.write("track name=hairpins description='Location of Hairpins' itemRgb='On'");
         outBed.newLine();
         new Random();
         String genome = "";
         int length = sequences.size() * 100 + 100;
         int col = true;

         for(int i = 100; i < length; ++i) {
            if (i % 100 != 0) {
               outGenome.write("N");
            } else {
               int idx = (i - 100) / 100;
               BuildFakeHalletGenome.Sequence s = (BuildFakeHalletGenome.Sequence)sequences.get(idx);
               outIdx.write(i + "\t" + s.getSeq() + "\t" + s.getCloneIDString());
               outIdx.newLine();
               outBed.write("gnl|TRC|hairpin_chromosome1\t" + (i + 1) + "\t" + (i + 1 + 20) + "\t" + s.getCloneIDString() + "\t1000\t+\t" + (i + 1) + "\t" + (i + 1 + 20) + "\t255,0,0");
               outBed.newLine();
               outBed.write("gnl|TRC|hairpin_chromosome1\t" + (i + 22) + "\t" + (i + 22 + 21) + "\t" + s.getCloneIDString() + "\t1000\t+\t" + (i + 22) + "\t" + (i + 22 + 21) + "\t0,255,255");
               outBed.newLine();
               outBed.write("gnl|TRC|hairpin_chromosome1\t" + (i + 22 + 21 + 1) + "\t" + (i + 22 + 21 + 27) + "\t" + s.getCloneIDString() + "\t1000\t+\t" + (i + 22 + 21 + 1) + "\t" + (i + 22 + 48) + "\t0,255,255");
               outBed.newLine();

               for(int j = 0; j < s.getSeq().length(); ++j) {
                  outGenome.write(s.getSeq().charAt(j));
               }

               i += s.getSeq().length() - 1;
            }
         }

         outGenome.close();
         outIdx.close();
         outBed.close();
      } catch (IOException var13) {
         var13.printStackTrace();
      }

   }

   private Vector seqHashToVector(Hashtable seq) {
      Vector v = new Vector();
      Iterator it = seq.values().iterator();

      while(it.hasNext()) {
         v.add(it.next());
      }

      return v;
   }

   private Hashtable readFile(String inFile) {
      Hashtable h = new Hashtable();

      try {
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));
         String line = in.readLine();
         in.readLine();

         String[] data;
         for(; (line = in.readLine()) != null; this.clone2seq.put(data[13], data[12])) {
            data = line.split("\t");
            if (!h.containsKey(data[12])) {
               h.put(data[12], new BuildFakeHalletGenome.Sequence(data[12], data[13], "1"));
            } else {
               BuildFakeHalletGenome.Sequence s = (BuildFakeHalletGenome.Sequence)h.get(data[12]);
               s.addClone_ID(data[13]);
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

   public static void main(String[] args) {
      new BuildFakeHalletGenome();
   }

   class Sequence {
      private String seq;
      private Vector clone_id;
      private int subpool;

      public Sequence(String seq, String clone_id, String pool) {
         this.seq = seq + "CTCGAG" + this.reverseComplement(seq.toUpperCase()).substring(0, 5);
         this.subpool = Integer.parseInt(pool);
         this.clone_id = new Vector();
         this.clone_id.add(clone_id);
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
         return this.subpool == ((BuildFakeHalletGenome.Sequence)s2).getSubpool();
      }

      public String reverseComplement(String seq) {
         String retStr = "";

         for(int i = seq.length() - 1; i >= 0; --i) {
            if (seq.charAt(i) == 'A') {
               retStr = retStr + "T";
            } else if (seq.charAt(i) == 'C') {
               retStr = retStr + "G";
            } else if (seq.charAt(i) == 'G') {
               retStr = retStr + "C";
            } else if (seq.charAt(i) == 'T') {
               retStr = retStr + "A";
            }
         }

         return retStr;
      }
   }
}
