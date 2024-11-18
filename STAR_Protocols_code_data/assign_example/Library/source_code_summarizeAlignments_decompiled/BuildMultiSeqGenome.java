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

public class BuildMultiSeqGenome {
   private String seqFile = "/Volumes/Bioinf_Data1/SeqData/Ketela_May2012/Project_Ketela_C/All80k_human.txt";
   private String genomeOutFile = "/Volumes/Bioinf_Data1/SeqData/Ketela_May2012/Project_Ketela_A/Sample_S1R1_High/hairpins_may2012_rev.fasta";
   private String idxOutFile = "/Volumes/Bioinf_Data1/SeqData/Ketela_May2012/Project_Ketela_A/Sample_S1R1_High/hairpins_may2012_rev.idx";
   private String bedOutFile = "/Volumes/Bioinf_Data1/SeqData/Ketela_May2012/Project_Ketela_A/Sample_S1R1_High/hairpins_may2012_rev.bed";
   private Hashtable clone2seq = new Hashtable();

   public BuildMultiSeqGenome() {
      Hashtable seqs = this.readFile(this.seqFile);
      Vector seqV = this.seqHashToVector(seqs);
      this.buildGenome(seqV);
   }

   private void buildFastaFile(Vector sequences) {
      int seqNum = 0;

      try {
         BufferedWriter outGenome = new BufferedWriter(new FileWriter(new File(this.genomeOutFile)));

         for(Iterator it = sequences.iterator(); it.hasNext(); ++seqNum) {
            BuildMultiSeqGenome.Sequence s = (BuildMultiSeqGenome.Sequence)it.next();
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
         new Random();
         String genome = "";
         int length = sequences.size() * 100 + 100;
         int col = true;
         Iterator it = sequences.iterator();

         while(it.hasNext()) {
            BuildMultiSeqGenome.Sequence s = (BuildMultiSeqGenome.Sequence)it.next();
            outGenome.write(">" + s.getClone_id());
            outGenome.newLine();
            outGenome.write(s.getSeq());
            outGenome.newLine();
         }

         outGenome.close();
      } catch (IOException var9) {
         var9.printStackTrace();
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

         while((line = in.readLine()) != null) {
            String[] data = line.split("\t");
            if (data[4].equals("Y")) {
               if (!h.containsKey(data[2])) {
                  h.put(data[2], new BuildMultiSeqGenome.Sequence(data[2], data[0], data[5]));
               } else {
                  BuildMultiSeqGenome.Sequence s = (BuildMultiSeqGenome.Sequence)h.get(data[2]);
                  s.addClone_ID(data[0]);
               }

               this.clone2seq.put(data[0], data[2]);
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
      new BuildMultiSeqGenome();
   }

   class Sequence {
      private String seq;
      private Vector clone_id;
      private int subpool;

      public Sequence(String seq, String clone_id, String pool) {
         this.seq = "G" + seq + "C";
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
         return this.subpool == ((BuildMultiSeqGenome.Sequence)s2).getSubpool();
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
