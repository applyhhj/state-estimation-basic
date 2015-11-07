package thu.instcloud.app.se.common;

/**
 * Created on 2015/11/7.
 */
public class EstimationOption {

    private boolean verbose;

    public EstimationOption() {

        verbose = true;

    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

}
