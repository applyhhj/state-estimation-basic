package thu.instcloud.app.se.mpdata;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

public class BranchData {

    private static Logger logger = LoggerFactory.getLogger(BranchData.class);

    private int[] i;

    private int[] j;

    private double[] r;

    private double[] x;

    private double[] b;

    private double[] rateA;

    private double[] rateB;

    private double[] rateC;

    private double[] ratio;

    private double[] angle;

    private int[] status;

    private double[] angmin;

    private double[] angmax;

    private double[] PF;

    private double[] QF;

    private double[] PT;

    private double[] QT;

    private int paraNum;

    private int n;

    public BranchData() {

        paraNum = 17;

        n = 0;

    }

    public boolean loadData(List<String> dataStr) {

        String[] cols;

        int ntmp = dataStr.size();

        int[] itmp = new int[ntmp];

        int[] jtmp = new int[ntmp];

        double[] rtmp = new double[ntmp];

        double[] xtmp = new double[ntmp];

        double[] btmp = new double[ntmp];

        double[] rateAtmp = new double[ntmp];

        double[] rateBtmp = new double[ntmp];

        double[] rateCtmp = new double[ntmp];

        double[] ratiotmp = new double[ntmp];

        double[] angletmp = new double[ntmp];

        int[] statustmp = new int[ntmp];

        double[] angmintmp = new double[ntmp];

        double[] angmaxtmp = new double[ntmp];

        double[] PFtmp = new double[ntmp];

        double[] QFtmp = new double[ntmp];

        double[] PTtmp = new double[ntmp];

        double[] QTtmp = new double[ntmp];


        for (int j = 0; j < dataStr.size(); j++) {

            cols = dataStr.get(j).trim().split(" +");

            if (cols.length != paraNum) {

                logger.error("Incorrect data format!");

                return false;

            }

            itmp[j] = Integer.parseInt(cols[0]);

            jtmp[j] = Integer.parseInt(cols[1]);

            rtmp[j] = Double.parseDouble(cols[2]);

            xtmp[j] = Double.parseDouble(cols[3]);

            btmp[j] = Double.parseDouble(cols[4]);

            rateAtmp[j] = Double.parseDouble(cols[5]);

            rateBtmp[j] = Double.parseDouble(cols[6]);

            rateCtmp[j] = Double.parseDouble(cols[7]);

            ratiotmp[j] = Double.parseDouble(cols[8]);

//            convert to radius
            angletmp[j] = Double.parseDouble(cols[9]) / 180 * Math.PI;

            statustmp[j] = Integer.parseInt(cols[10]);

            angmintmp[j] = Double.parseDouble(cols[11]);

            angmaxtmp[j] = Double.parseDouble(cols[12]);

            PFtmp[j] = Double.parseDouble(cols[13]);

            QFtmp[j] = Double.parseDouble(cols[14]);

            PTtmp[j] = Double.parseDouble(cols[15]);

            QTtmp[j] = Double.parseDouble(cols[16]);

        }

        setN(ntmp);

        setI(itmp);

        setJ(jtmp);

        setR(rtmp);

        setX(xtmp);

        setB(btmp);

        setRateA(rateAtmp);

        setRateB(rateBtmp);

        setRateC(rateCtmp);

        setRatio(ratiotmp);

        setAngle(angletmp);

        setStatus(statustmp);

        setAngmin(angmintmp);

        setAngmax(angmaxtmp);

        setPF(PFtmp);

        setQF(QFtmp);

        setPT(PTtmp);

        setQT(QTtmp);

        return true;

    }

    public int[] getI() {
        return i;
    }

    public void setI(int[] i) {
        this.i = i;
    }

    public int[] getJ() {
        return j;
    }

    public void setJ(int[] j) {
        this.j = j;
    }

    public double[] getR() {
        return r;
    }

    public void setR(double[] r) {
        this.r = r;
    }

    public double[] getX() {
        return x;
    }

    public void setX(double[] x) {
        this.x = x;
    }

    public double[] getB() {
        return b;
    }

    public void setB(double[] b) {
        this.b = b;
    }

    public double[] getRateA() {
        return rateA;
    }

    public void setRateA(double[] rateA) {
        this.rateA = rateA;
    }

    public double[] getRateB() {
        return rateB;
    }

    public void setRateB(double[] rateB) {
        this.rateB = rateB;
    }

    public double[] getRateC() {
        return rateC;
    }

    public void setRateC(double[] rateC) {
        this.rateC = rateC;
    }

    public double[] getRatio() {
        return ratio;
    }

    public void setRatio(double[] ratio) {
        this.ratio = ratio;
    }

    public double[] getAngle() {
        return angle;
    }

    public void setAngle(double[] angle) {
        this.angle = angle;
    }

    public int[] getStatus() {
        return status;
    }

    public void setStatus(int[] status) {
        this.status = status;
    }

    public double[] getAngmin() {
        return angmin;
    }

    public void setAngmin(double[] angmin) {
        this.angmin = angmin;
    }

    public double[] getAngmax() {
        return angmax;
    }

    public void setAngmax(double[] angmax) {
        this.angmax = angmax;
    }

    public double[] getPF() {
        return PF;
    }

    public void setPF(double[] PF) {
        this.PF = PF;
    }

    public double[] getQF() {
        return QF;
    }

    public void setQF(double[] QF) {
        this.QF = QF;
    }

    public double[] getPT() {
        return PT;
    }

    public void setPT(double[] PT) {
        this.PT = PT;
    }

    public double[] getQT() {
        return QT;
    }

    public void setQT(double[] QT) {
        this.QT = QT;
    }

    public int getParaNum() {
        return paraNum;
    }

    public void setParaNum(int paraNum) {
        this.paraNum = paraNum;
    }

    public int getN() {
        return n;
    }

    public void setN(int n) {
        this.n = n;
    }

}
