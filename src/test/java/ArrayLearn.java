import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Created on 2015/11/6.
 */
public class ArrayLearn {

    public static void main(String[] args) {

        long start = System.currentTimeMillis();

        Logger logger = LoggerFactory.getLogger(ArrayLearn.class);

        System.out.print(System.currentTimeMillis() - start);

        double[][] arr2d = new double[10][5];

        System.out.print("\nrows " + arr2d.length);

        System.out.print("\ncols " + arr2d[0].length);

    }

}
