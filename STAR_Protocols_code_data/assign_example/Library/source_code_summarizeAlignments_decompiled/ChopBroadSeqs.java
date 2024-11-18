import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ChopBroadSeqs {
   private String inFile1 = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Broad Sequencing Data/90KRawData/s_1_sequence.txt";
   private String inFile2 = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Broad Sequencing Data/90KRawData/s_2_sequence.txt";
   private String outFile1 = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Broad Sequencing Data/90KRawData/s_1_sequence_chopped.txt";
   private String outFile2 = "/Users/kbrown/Documents/Work/SequenceOfLibrary/Broad Sequencing Data/90KRawData/s_2_sequence_chopped.txt";

   public ChopBroadSeqs() {
      this.readAndChop(this.inFile1, this.outFile1);
      this.readAndChop(this.inFile2, this.outFile2);
   }

   private void readAndChop(String inFile, String outFile) {
      String line = null;

      try {
         BufferedReader in = new BufferedReader(new FileReader(new File(inFile)));
         BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));

         while((line = in.readLine()) != null) {
            String line2 = in.readLine();
            String line3 = in.readLine();
            String line4 = in.readLine();
            out.write(line);
            out.newLine();
            out.write(line2.substring(12, 33));
            out.newLine();
            out.write(line3);
            out.newLine();
            out.write(line4.substring(12, 33));
            out.newLine();
         }

         out.close();
         in.close();
      } catch (FileNotFoundException var9) {
         var9.printStackTrace();
      } catch (IOException var10) {
         var10.printStackTrace();
      }

   }

   public static void main(String[] args) {
      new ChopBroadSeqs();
   }
}
