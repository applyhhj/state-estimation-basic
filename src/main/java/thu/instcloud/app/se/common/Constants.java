package thu.instcloud.app.se.common;

/**
 * Created on 2015/11/6.
 */
public abstract class Constants {

    public static final double E = 1e-5;

    public abstract class CDF {
        public static final String BUS_SECTION = "BUS DATA FOLLOWS";
        public static final String BRANCH_SECTION = "BRANCH DATA FOLLOWS";
        public static final String MEASUREMENT_SECTION = "MEASUREMENT SYSTEM DATA FOLLOWS";
        public static final String BUS_SECTION_END = "-999";
        public static final String BRANCH_SECTION_END = "-999";
        public static final String MEASUREMENT_SECTION_END = "-999";
        public static final int BRANCH_SECTION_NUMBER_OF_COLUMN = 21;
        public static final int MEASUREMENT_SECTION_NUMBER_OF_COLUMN = 6;
    }

    public abstract class MPC {

        public static final String MPC_VERSION = "mpc.version";
        public static final String MPC_VERSION_END = "mpc.version.end";
        public static final String MPC_BASEMVA = "mpc.baseMVA";
        public static final String MPC_BASEMVA_END = "mpc.baseMVA.end";
        public static final String MPC_BUS = "mpc.bus";
        public static final String MPC_BUS_END = "mpc.bus.end";
        public static final String MPC_GEN = "mpc.gen";
        public static final String MPC_GEN_END = "mpc.gen.end";
        public static final String MPC_BRANCH = "mpc.branch";
        public static final String MPC_BRANCH_END = "mpc.branch.end";

        public abstract class BusTypes {

            public static final int PQ = 1;
            public static final int PV = 2;
            public static final int REF = 3;
            public static final int ISO = 4;

        }

    }

    public abstract class ESTIMATOR {

        public static final double ERR_REC = 10;
        public static final double ERR_THETA = 1e-3;
        public static final double ERR_V = 1e-2;

        public static final double TOL = 1e-8;

    }
}
