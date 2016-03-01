import Utilities.BitUtilities;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Scanner;

public class FileLoader {

    // Load finger print data into FingerPrint class
    public FingerPrint[] loadFingerprintData(String filename, int size) {
        FingerPrint[] objFP = new FingerPrint[size];
        FileInputStream inputStream = null;
        Scanner sc = null;
        BitUtilities objBU = new BitUtilities();
        String[] splittedLine;
        long start = System.currentTimeMillis();

        int i = 0;

        try {
            inputStream = new FileInputStream(filename);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        sc = new Scanner(inputStream, "UTF-8");
        while(sc.hasNextLine()) {
            ++i;
            FingerPrint obj = new FingerPrint();
            splittedLine = sc.nextLine().split("\\ ");
            obj.idx = i;
            obj.moleculeName = splittedLine[0].trim();
            obj.fp = objBU.getBitSet(splittedLine[1].trim());
            obj.card = (int)objBU.countOnBits(obj.fp);
            obj.bitfreqs = objBU.countBits(obj.fp.getBits());
            obj.hashcode = obj.fp.hashCode();
            obj.bytebitfreqs = objBU.countBitsPerByte(obj.fp.getBits());
            objFP[i-1] = obj;

           // System.err.println(i);
        }
        try {
            inputStream.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


        long end = System.currentTimeMillis();
        System.err.println("Time taken to load : " + (end-start)/1000.0 + " seconds");


        return objFP;
    }


}
