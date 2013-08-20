import java.io.File;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Arrays;
import java.io.InputStreamReader;
import java.io.BufferedReader;

public class RuntimeDemo3 {
	public static void pipeStream(InputStream input, OutputStream output)
	   throws IOException
	{
	   byte buffer[] = new byte[1024];
	   int numRead = 0;

	   do
	   {
		  numRead = input.read(buffer);
		  output.write(buffer, 0, numRead);
	   } while (input.available() > 0);

	   output.flush();
	}

	public static void main(String[] argv)
	{
	   FileInputStream fileIn = null;
	   FileOutputStream fileOut = null;

	   OutputStream procIn = null;
	   InputStream procOut = null;

	   try
	   {
         //fileIn = new FileInputStream("1351_LPAIH7N1_log.txt");
         fileOut = new FileOutputStream("testOut.txt");
         String btcexe = "./btcutils ";
         String bamfile = "-bam SamTestFiles/S1_refHPAI_cons_stampy.bam ";
         String ref = "-ref SamTestFiles/refHPAI_cons.fa ";
         String output = "-stub out";
         // print a message
         String cmd = btcexe + bamfile + ref + output;
         // print a message
         System.out.println("Executing btcutils");
         System.out.println(cmd);
         //System.out.println(Arrays.toString(cmdArray));
         // create a process and execute cmdArray
         Process process = Runtime.getRuntime().exec(cmd);
         String line;
         BufferedReader input = new BufferedReader(new InputStreamReader(process.getInputStream()));
         while ((line = input.readLine()) != null) {
           System.out.println(line);
         }
         input.close();
         // Process process = Runtime.getRuntime().exec ("/Users/josephhughes/Documents/Repo/btctools/btcutils");
         //procIn = process.getOutputStream();
         //procOut = process.getInputStream();
         
         //pipeStream(fileIn, procIn);
         //pipeStream(procOut, fileOut);
	   }
	   catch (IOException ioe)
	   {
		  System.out.println(ioe);
	   }
	}
}