package thu.instcloud.app.se.mpdata;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import thu.instcloud.app.se.common.Constants;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created on 2015/11/6.
 */
public class BusData {

    private static Logger logger = LoggerFactory.getLogger(BusData.class);

    private Map<Integer, Integer> TIO;

    private Map<Integer, Integer> TOI;

    private Map<Integer, Integer> TOA;

    private int[] number;

    private int[] type;

    private double[] PD;

    private double[] QD;

    private double[] Gs;

    private double[] Bs;

    private int[] area;

    private double[] voltage;

    private double[] angle;

    private double[] VBaseKV;

    private int[] zone;

    private double[] VMax;

    private double[] VMin;

    private int paraNum;

    private int n;

    private int nrefI;

    public BusData() {

        TIO = new HashMap<Integer, Integer>();

        TOI = new HashMap<Integer, Integer>();

        TOA = new HashMap<Integer, Integer>();

        paraNum = 13;

        n = 0;

    }

    public boolean loadData(List<String> dataStr) {

        String[] cols;

        int ntmp = dataStr.size();

        int[] numbertmp = new int[ntmp];

        int[] typetmp = new int[ntmp];

        double[] PDtmp = new double[ntmp];

        double[] QDtmp = new double[ntmp];

        double[] Gstmp = new double[ntmp];

        double[] Bstmp = new double[ntmp];

        int[] areatmp = new int[ntmp];

        double[] voltagetmp = new double[ntmp];

        double[] angletmp = new double[ntmp];

        double[] VBasetmp = new double[ntmp];

        int[] zonetmp = new int[ntmp];

        double[] Vmaxtmp = new double[ntmp];

        double[] Vmintmp = new double[ntmp];

        for (int i = 0; i < dataStr.size(); i++) {

            cols = dataStr.get(i).trim().split(" +");

            if (cols.length != paraNum) {

                logger.error("Incorrect data format!");

                return false;

            }

            numbertmp[i] = Integer.parseInt(cols[0]);

            typetmp[i] = Integer.parseInt(cols[1]);

            PDtmp[i] = Double.parseDouble(cols[2]);

            QDtmp[i] = Double.parseDouble(cols[3]);

            Gstmp[i] = Double.parseDouble(cols[4]);

            Bstmp[i] = Double.parseDouble(cols[5]);

            areatmp[i] = Integer.parseInt(cols[6]);

            voltagetmp[i] = Double.parseDouble(cols[7]);

//            convert to radius
            angletmp[i] = Double.parseDouble(cols[8]) / 180 * Math.PI;

            VBasetmp[i] = Double.parseDouble(cols[9]);

            zonetmp[i] = Integer.parseInt(cols[10]);

            Vmaxtmp[i] = Double.parseDouble(cols[11]);

            Vmintmp[i] = Double.parseDouble(cols[12]);

        }

        setNumber(numbertmp);

        setType(typetmp);

        setPD(PDtmp);

        setQD(QDtmp);

        setGs(Gstmp);

        setBs(Bstmp);

        setArea(areatmp);

        setVoltage(voltagetmp);

        setAngle(angletmp);

        setVBaseKV(VBasetmp);

        setZone(zonetmp);

        setVMax(Vmaxtmp);

        setVMin(Vmintmp);

        setN(ntmp);

        return true;

    }

    //    number starts from 1
    public void reorderBusNumbers(BranchData branchData) {

        int[] i = branchData.getI();

        int[] j = branchData.getJ();

        int[] type = getType();

        int[] buses = getNumber();

        Map<Integer, Integer> busBranchNumberMap = new HashMap<Integer, Integer>();

        int ni, nj, ntmp;

//        compute lines from each bus
        for (int k = 0; k < i.length; k++) {

            ni = i[k];

            nj = j[k];

            if (ni == nj) {

                continue;

            }

            if (!busBranchNumberMap.containsKey(ni)) {

                busBranchNumberMap.put(ni, 1);

            } else {

                busBranchNumberMap.put(ni, busBranchNumberMap.get(ni) + 1);

            }


            if (!busBranchNumberMap.containsKey(nj)) {

                busBranchNumberMap.put(nj, 1);

            } else {

                busBranchNumberMap.put(nj, busBranchNumberMap.get(nj) + 1);

            }

        }

        int idx = busBranchNumberMap.size();

        Map<Integer, Integer> pvbuses = new HashMap<Integer, Integer>();

        for (int k = 0; k < type.length; k++) {

            TOA.put(buses[k], k);

//            swing bus, reference bus
            if (type[k] == Constants.MPC.BusTypes.REF) {

                ntmp = buses[k];

                if (busBranchNumberMap.containsKey(ntmp)) {

                    TIO.put(idx, ntmp);

                    nrefI = idx;

                    busBranchNumberMap.remove(ntmp);

                    idx--;

                } else {

                    logger.error("Reference bus {} does not exist in branch data!", ntmp);

                    return;

                }

            } else if (type[k] == Constants.MPC.BusTypes.PV) {

//                PV bus
                ntmp = buses[k];

                if (busBranchNumberMap.containsKey(ntmp)) {

                    pvbuses.put(ntmp, busBranchNumberMap.get(ntmp));

                    busBranchNumberMap.remove(ntmp);

                }

            }

        }

        idx = addBuses(idx, pvbuses);

        addBuses(idx, busBranchNumberMap);

        for (Map.Entry<Integer, Integer> e : TIO.entrySet()) {

            TOI.put(e.getValue(), e.getKey());

        }

    }

    private int addBuses(int currentIdx, Map<Integer, Integer> busBranchNoMap) {

//        insert from the end
        while (busBranchNoMap.size() > 0) {

            Integer maxKey = null;

            for (Map.Entry<Integer, Integer> e : busBranchNoMap.entrySet()) {

                if (maxKey == null) {

                    maxKey = e.getKey();

                }

                if (busBranchNoMap.get(maxKey) < e.getValue()) {

                    maxKey = e.getKey();

                }

            }

            TIO.put(currentIdx--, maxKey);

            busBranchNoMap.remove(maxKey);

        }

        return currentIdx;

    }


    public Map<Integer, Integer> getTIO() {
        return TIO;
    }

    public void setTIO(Map<Integer, Integer> TIO) {
        this.TIO = TIO;
    }

    public Map<Integer, Integer> getTOI() {
        return TOI;
    }

    public void setTOI(Map<Integer, Integer> TOI) {
        this.TOI = TOI;
    }

    public Map<Integer, Integer> getTOA() {
        return TOA;
    }

    public void setTOA(Map<Integer, Integer> TOA) {
        this.TOA = TOA;
    }

    public int[] getNumber() {
        return number;
    }

    public void setNumber(int[] number) {
        this.number = number;
    }

    public int[] getType() {
        return type;
    }

    public void setType(int[] type) {
        this.type = type;
    }

    public double[] getPD() {
        return PD;
    }

    public void setPD(double[] PD) {
        this.PD = PD;
    }

    public double[] getQD() {
        return QD;
    }

    public void setQD(double[] QD) {
        this.QD = QD;
    }

    public double[] getGs() {
        return Gs;
    }

    public void setGs(double[] gs) {
        Gs = gs;
    }

    public double[] getBs() {
        return Bs;
    }

    public void setBs(double[] bs) {
        Bs = bs;
    }

    public int[] getArea() {
        return area;
    }

    public void setArea(int[] area) {
        this.area = area;
    }

    public double[] getVoltage() {
        return voltage;
    }

    public void setVoltage(double[] voltage) {
        this.voltage = voltage;
    }

    public double[] getAngle() {
        return angle;
    }

    public void setAngle(double[] angle) {
        this.angle = angle;
    }

    public double[] getVBaseKV() {
        return VBaseKV;
    }

    public void setVBaseKV(double[] VBaseKV) {
        this.VBaseKV = VBaseKV;
    }

    public int[] getZone() {
        return zone;
    }

    public void setZone(int[] zone) {
        this.zone = zone;
    }

    public double[] getVMax() {
        return VMax;
    }

    public void setVMax(double[] VMax) {
        this.VMax = VMax;
    }

    public double[] getVMin() {
        return VMin;
    }

    public void setVMin(double[] VMin) {
        this.VMin = VMin;
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

    public int getNrefI() {
        return nrefI;
    }

}
