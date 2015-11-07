package thu.instcloud.app.se.mpdata;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

/**
 * Created on 2015/11/6.
 */
public class GeneratorData {

    private static Logger logger = LoggerFactory.getLogger(GeneratorData.class);

    private int[] number;

    private double[] Pg;

    private double[] Qg;

    private double[] Qmax;

    private double[] Qmin;

    private double[] Vg;

    private double[] mBase;

    private int[] status;

    private double[] Pmax;

    private double[] Pmin;

    private double[] Pc1;

    private double[] Pc2;

    private double[] Qc1min;

    private double[] Qc1max;

    private double[] Qc2min;

    private double[] Qc2max;

    private double[] ramp_agc;

    private double[] ramp_10;

    private double[] ramp_30;

    private double[] ramp_q;

    private double[] apf;

    private int n;

    private int paraNum;

    public GeneratorData() {

        paraNum = 21;

        n = 0;

    }

    public boolean loadData(List<String> dataStr) {

        String[] cols;

        int ntmp = dataStr.size();

        int[] numbertmp = new int[ntmp];

        double[] Pgtmp = new double[ntmp];

        double[] Qgtmp = new double[ntmp];

        double[] Qmaxtmp = new double[ntmp];

        double[] Qmintmp = new double[ntmp];

        double[] Vgtmp = new double[ntmp];

        double[] mBasetmp = new double[ntmp];

        int[] statustmp = new int[ntmp];

        double[] Pmaxtmp = new double[ntmp];

        double[] Pmintmp = new double[ntmp];

        double[] Pc1tmp = new double[ntmp];

        double[] Pc2tmp = new double[ntmp];

        double[] Qc1mintmp = new double[ntmp];

        double[] Qc1maxtmp = new double[ntmp];

        double[] Qc2mintmp = new double[ntmp];

        double[] Qc2maxtmp = new double[ntmp];

        double[] ramp_agctmp = new double[ntmp];

        double[] ramp_10tmp = new double[ntmp];

        double[] ramp_30tmp = new double[ntmp];

        double[] ramp_qtmp = new double[ntmp];

        double[] apftmp = new double[ntmp];

        for (int i = 0; i < dataStr.size(); i++) {

            cols = dataStr.get(i).trim().split(" +");

            if (cols.length != paraNum) {

                logger.error("Incorrect data format!");

                return false;

            }

            numbertmp[i] = Integer.parseInt(cols[0]);

            Pgtmp[i] = Double.parseDouble(cols[1]);

            Qgtmp[i] = Double.parseDouble(cols[2]);

            Qmaxtmp[i] = Double.parseDouble(cols[3]);

            Qmintmp[i] = Double.parseDouble(cols[4]);

            Vgtmp[i] = Double.parseDouble(cols[5]);

            mBasetmp[i] = Double.parseDouble(cols[6]);

            statustmp[i] = Integer.parseInt(cols[7]);

            Pmaxtmp[i] = Double.parseDouble(cols[8]);

            Pmintmp[i] = Double.parseDouble(cols[9]);

            Pc1tmp[i] = Double.parseDouble(cols[10]);

            Pc2tmp[i] = Double.parseDouble(cols[11]);

            Qc1mintmp[i] = Double.parseDouble(cols[12]);

            Qc1maxtmp[i] = Double.parseDouble(cols[13]);

            Qc2mintmp[i] = Double.parseDouble(cols[14]);

            Qc2maxtmp[i] = Double.parseDouble(cols[15]);

            ramp_agctmp[i] = Double.parseDouble(cols[16]);

            ramp_10tmp[i] = Double.parseDouble(cols[17]);

            ramp_30tmp[i] = Double.parseDouble(cols[18]);

            ramp_qtmp[i] = Double.parseDouble(cols[19]);

            apftmp[i] = Double.parseDouble(cols[20]);

        }

        setN(ntmp);

        setNumber(numbertmp);

        setPg(Pgtmp);

        setQg(Qgtmp);

        setQmax(Qmaxtmp);

        setQmin(Qmintmp);

        setVg(Vgtmp);

        setmBase(mBasetmp);

        setStatus(statustmp);

        setPmax(Pmaxtmp);

        setPmin(Pmintmp);

        setPc1(Pc1tmp);

        setPc2(Pc2tmp);

        setQc1min(Qc1mintmp);

        setQc1max(Qc1maxtmp);

        setQc2min(Qc2mintmp);

        setQc2max(Qc2maxtmp);

        setRamp_agc(ramp_agctmp);

        setRamp_10(ramp_10tmp);

        setRamp_30(ramp_30tmp);

        setRamp_q(ramp_qtmp);

        setApf(apftmp);

        return true;

    }

    public int[] getNumber() {
        return number;
    }

    public void setNumber(int[] number) {
        this.number = number;
    }

    public double[] getPg() {
        return Pg;
    }

    public void setPg(double[] pg) {
        Pg = pg;
    }

    public double[] getQg() {
        return Qg;
    }

    public void setQg(double[] qg) {
        Qg = qg;
    }

    public double[] getQmax() {
        return Qmax;
    }

    public void setQmax(double[] qmax) {
        Qmax = qmax;
    }

    public double[] getQmin() {
        return Qmin;
    }

    public void setQmin(double[] qmin) {
        Qmin = qmin;
    }

    public double[] getVg() {
        return Vg;
    }

    public void setVg(double[] vg) {
        Vg = vg;
    }

    public double[] getmBase() {
        return mBase;
    }

    public void setmBase(double[] mBase) {
        this.mBase = mBase;
    }

    public int[] getStatus() {
        return status;
    }

    public void setStatus(int[] status) {
        this.status = status;
    }

    public double[] getPmax() {
        return Pmax;
    }

    public void setPmax(double[] pmax) {
        Pmax = pmax;
    }

    public double[] getPmin() {
        return Pmin;
    }

    public void setPmin(double[] pmin) {
        Pmin = pmin;
    }

    public double[] getPc1() {
        return Pc1;
    }

    public void setPc1(double[] pc1) {
        Pc1 = pc1;
    }

    public double[] getPc2() {
        return Pc2;
    }

    public void setPc2(double[] pc2) {
        Pc2 = pc2;
    }

    public double[] getQc1min() {
        return Qc1min;
    }

    public void setQc1min(double[] qc1min) {
        Qc1min = qc1min;
    }

    public double[] getQc1max() {
        return Qc1max;
    }

    public void setQc1max(double[] qc1max) {
        Qc1max = qc1max;
    }

    public double[] getQc2min() {
        return Qc2min;
    }

    public void setQc2min(double[] qc2min) {
        Qc2min = qc2min;
    }

    public double[] getQc2max() {
        return Qc2max;
    }

    public void setQc2max(double[] qc2max) {
        Qc2max = qc2max;
    }

    public double[] getRamp_agc() {
        return ramp_agc;
    }

    public void setRamp_agc(double[] ramp_agc) {
        this.ramp_agc = ramp_agc;
    }

    public double[] getRamp_10() {
        return ramp_10;
    }

    public void setRamp_10(double[] ramp_10) {
        this.ramp_10 = ramp_10;
    }

    public double[] getRamp_30() {
        return ramp_30;
    }

    public void setRamp_30(double[] ramp_30) {
        this.ramp_30 = ramp_30;
    }

    public double[] getRamp_q() {
        return ramp_q;
    }

    public void setRamp_q(double[] ramp_q) {
        this.ramp_q = ramp_q;
    }

    public double[] getApf() {
        return apf;
    }

    public void setApf(double[] apf) {
        this.apf = apf;
    }

    public int getN() {
        return n;
    }

    public void setN(int n) {
        this.n = n;
    }

    public int getParaNum() {
        return paraNum;
    }

    public void setParaNum(int paraNum) {
        this.paraNum = paraNum;
    }
}
