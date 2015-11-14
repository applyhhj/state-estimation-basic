package thu.instcloud.app.se.estimator;

import thu.instcloud.app.se.common.EstimationOption;

import static thu.instcloud.app.se.common.Utils.Common.isLinux;

/**
 * Created on 2015/11/7.
 */
public class Executor {

    public static void main(String[] args) throws InterruptedException {

        String fpath;

        if (isLinux()) {

            fpath = "/home/hjh/doc/powersystem/4bus/ieee4cdftest.txt";

        } else {

            fpath = "F:\\projects\\data\\matpower-data-process\\data\\case89pegase.txt";

        }

        EstimationOption option = new EstimationOption();

        option.setVerbose(false);

        option.setDebug(true);

        PowerSystem powerSystem = new PowerSystem(fpath, option);

        int i = 0;

        long start;

        while (i++ <2) {

            start = System.currentTimeMillis();

            powerSystem.run();

            System.out.printf("\nEstimate %d duration: %d ms", i, System.currentTimeMillis() - start);

            powerSystem.printStateInExternalInPolarDegree();

            Thread.sleep(500);

        }

        System.out.print("\nDone!\n");

    }


}
