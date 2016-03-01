package Utilities;

import org.apache.lucene.util.OpenBitSet;

import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.List;


public class BitUtilities {


    public int[] countBitsPerByte(long[] arr) {

        int[] ret = new int[256];
        byte[] bytearray = new byte[256];
        byte[] tmp;
        for (int i = 0; i < arr.length; i++) {
            tmp = convertToByteArray(arr[i]);
            int idx = i * 8;
            for (int j = 0; j < tmp.length; j++) {
                ret[j + idx] = Integer.bitCount(tmp[j]);
            }
        }


        return ret;
    }

    public int[] countBits(long[] arr) {
        int[] ret = new int[arr.length];
        for (int i = 0; i < arr.length; i++) {
            ret[i] = Long.bitCount(arr[i]);
        }
        return ret;
    }

    public long countOnBits(String bitstring) throws InterruptedException {


        OpenBitSet bs = new OpenBitSet(2048);

        for (int i = 0; i < 2048; i++) {
            if (bitstring.charAt(i) == '1')
                bs.set(i);
        }
        return countOnBits(bs);
    }

    public long countOnBits(OpenBitSet bs) {
        return bs.cardinality();
    }

    public List<String> setStringIndices(String bitstring) {
        List<String> lst = new ArrayList<String>();
        for (int i = 0; i < 2048; i++) {
            if (bitstring.charAt(i) == '1') {
                lst.add(i + "");
            }
        }
        return lst;
    }

    public OpenBitSet getBitSet(String bitString) {
        OpenBitSet bs = new OpenBitSet(2048);

        for (int i = 0; i < 2048; i++) {
            if (bitString.charAt(i) == '1')
                bs.set(i);
        }

        return bs;

    }

    public static byte[] convertToByteArray(long value) {

        byte[] bytes = new byte[8];
        ByteBuffer buffer = ByteBuffer.allocate(bytes.length);
        buffer.putLong(value);
        return buffer.array();
    }

}
