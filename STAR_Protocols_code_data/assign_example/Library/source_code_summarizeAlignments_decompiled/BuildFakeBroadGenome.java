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

public class BuildFakeBroadGenome {
   private String seqFile = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Broad Sequencing Data/90KRawData/90K_hairpin_list.txt";
   private String genomeOutFile = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Broad Sequencing Data/90KRawData/hairpin_genome.fasta";
   private String idxOutFile = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Broad Sequencing Data/90KRawData/hairpin_genome.idx";
   private Hashtable clone2seq = new Hashtable();

   public BuildFakeBroadGenome() {
      Hashtable seqs = this.readFile(this.seqFile);
      Vector seqV = this.seqHashToVector(seqs);
      this.buildGenome(seqV);
   }

   private void buildGenome(Vector sequences) {
      try {
         BufferedWriter outGenome = new BufferedWriter(new FileWriter(new File(this.genomeOutFile)));
         BufferedWriter outIdx = new BufferedWriter(new FileWriter(new File(this.idxOutFile)));
         outGenome.write(">gnl|TRC|hairpin_chromosome1");
         outGenome.newLine();
         outIdx.write("Genome Location\tSequence\tMatching TRC.ID(s)");
         outIdx.newLine();
         Random rnd = new Random();
         String genome = "";
         int length = sequences.size() * 100 + 100;
         int col = true;

         for(int i = 100; i < length; ++i) {
            int base;
            if (i % 100 != 0) {
               base = rnd.nextInt(4);
               if (base == 0) {
                  outGenome.write("a");
               } else if (base == 1) {
                  outGenome.write("c");
               } else if (base == 2) {
                  outGenome.write("g");
               } else {
                  outGenome.write("t");
               }
            } else {
               base = (i - 100) / 100;
               BuildFakeBroadGenome.Sequence s = (BuildFakeBroadGenome.Sequence)sequences.get(base);
               outIdx.write(i - 99 + "\t" + s.getSeq() + "\t" + s.getCloneIDString());
               outIdx.newLine();

               for(int j = 0; j < s.getSeq().length(); ++j) {
                  outGenome.write(s.getSeq().charAt(j));
               }

               i += s.getSeq().length() - 1;
            }
         }

         outGenome.close();
         outIdx.close();
      } catch (IOException var12) {
         var12.printStackTrace();
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

         String[] data;
         for(String line = in.readLine(); (line = in.readLine()) != null; this.clone2seq.put(data[0], data[1])) {
            data = line.split("\t");
            if (!h.containsKey(data[1])) {
               h.put(data[1], new BuildFakeBroadGenome.Sequence(data[1], data[0]));
            } else {
               BuildFakeBroadGenome.Sequence s = (BuildFakeBroadGenome.Sequence)h.get(data[1]);
               s.addClone_ID(data[0]);
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
      new BuildFakeBroadGenome();
   }

   class Sequence {
      private String seq;
      private Vector clone_id;
      private int subpool;

      public Sequence(String seq, String clone_id) {
         this.seq = seq;
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
         return this.subpool == ((BuildFakeBroadGenome.Sequence)s2).getSubpool();
      }
   }
}
