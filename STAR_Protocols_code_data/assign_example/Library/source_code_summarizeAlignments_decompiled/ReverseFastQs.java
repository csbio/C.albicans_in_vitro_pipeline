import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class ReverseFastQs {
   public static String reverseComplement(String seq) {
      String retStr = "";

      for(int i = seq.length() - 1; i >= 0; --i) {
         if (seq.charAt(i) == 'A') {
            retStr = retStr.concat("T");
         } else if (seq.charAt(i) == 'C') {
            retStr = retStr.concat("G");
         } else if (seq.charAt(i) == 'G') {
            retStr = retStr.concat("C");
         } else if (seq.charAt(i) == 'T') {
            retStr = retStr.concat("A");
         } else if (seq.charAt(i) == 'N') {
            retStr = retStr.concat("N");
         }
      }

      return retStr;
   }

   public static void main(String[] args) throws IOException {
      int lineno = 0;

      String l;
      for(BufferedReader fp = new BufferedReader(new InputStreamReader(System.in)); (l = fp.readLine()) != null; System.out.println(l)) {
         if (lineno == 1) {
            l = reverseComplement(l);
         }

         if (lineno == 3) {
            l = new String((new StringBuffer(l)).reverse());
            lineno = 0;
         } else {
            ++lineno;
         }
      }

   }
}
