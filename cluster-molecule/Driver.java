
import org.apache.lucene.util.OpenBitSet;

/**
 * Cluster molecule based on tanimoto similarity. Molecules can be seperated into K parts specified by user, which can then be processed in parallel.
 */
public class Driver {
    public static void main(String[] args) {
        if(args.length != 5) {
            System.out.println("(Usage) java <class_name> fingerprint.fps num_elements threshold");
            System.exit(1);
        }
        String filename = args[0];
        int size = Integer.parseInt(args[1]);
        double threshold = Double.parseDouble(args[2]);
        int start_sample = Integer.parseInt(args[3]);
        int end_sample = Integer.parseInt(args[4]);
        
        //Load finger print data
        FileLoader objFL = new FileLoader();
        FingerPrint[] objData = objFL.loadFingerprintData(filename, size);

        // time the blocking step.
        long start = System.currentTimeMillis();
        // count the number of pairs
        long pairs = 0;


        // perform blocking based on tanimoto similarity
        double intersectCard;
        double t;
        double c = threshold/(1+threshold);
        for(int i = start_sample-1; i < end_sample - 1; i ++) {
            double d1 = objData[i].card;

            for(int j = i + 1; j < size; j ++) {
		        double d2 = objData[j].card;
                t = c*(d1+d2);
                if(objData[i].hashcode == objData[j].hashcode) {
                    pairs ++;
                    continue;
                }
                if(Math.min(d1, d2) >= t) {
                    double s = 0;
                    for(int k = 0; k < 32; k ++) {
                        s = s + (Math.min(objData[i].bitfreqs[k], objData[j].bitfreqs[k]));
                    }
                    if(s >= t) {
		                intersectCard = OpenBitSet.intersectionCount(objData[i].fp, objData[j].fp);
                        if (intersectCard >= t) {
                            System.out.println( i +", " +j);
                            pairs++;
                        }
                    }
                }

            }

        }

        long end = System.currentTimeMillis();
        //System.err.println("Num pairs : " + pairs);
        //System.err.println("Time taken: "+(end-start));





    }
}
