package thu.instcloud.app.se.common;

/**
 * Created on 2015/11/7.
 */
public class EstimationOption {

    private boolean verbose;

    private boolean debug;

    public EstimationOption() {

        verbose = true;

        debug = false;

    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public boolean isDebug() {
        return debug;
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }
}
