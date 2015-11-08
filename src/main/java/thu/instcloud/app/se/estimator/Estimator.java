package thu.instcloud.app.se.estimator;

import org.ojalgo.access.Access2D;
import org.ojalgo.function.ComplexFunction;
import org.ojalgo.function.PrimitiveFunction;
import org.ojalgo.matrix.BasicMatrix;
import org.ojalgo.matrix.ComplexMatrix;
import org.ojalgo.scalar.ComplexNumber;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import thu.instcloud.app.se.common.Constants;
import thu.instcloud.app.se.common.OjMatrixManipulator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static thu.instcloud.app.se.common.Utils.MatrixExtension.*;
import static thu.instcloud.app.se.common.Utils.OJ.*;
import static thu.instcloud.app.se.common.Utils.hasDuplicateElement;

/**
 * Created on 2015/11/7.
 */
public class Estimator extends OjMatrixManipulator {

    private static Logger logger = LoggerFactory.getLogger(Estimator.class);

    private BasicMatrix dSbDvmCplx;

    private BasicMatrix dSbDvaCplx;

    private BasicMatrix dSfDvmCplx;

    private BasicMatrix dSfDvaCplx;

    private BasicMatrix dStDvmCplx;

    private BasicMatrix dStDvaCplx;

    private BasicMatrix SfCplx;

    private BasicMatrix StCplx;

    private BasicMatrix HFReal;

    private BasicMatrix WInvReal;

    private BasicMatrix VpfNormCplx;

    private BasicMatrix VpfNormMatrixCplx;

    private BasicMatrix VpfCplx;

    private BasicMatrix VpfMatrixCplx;

    private PowerSystem powerSystem;

    private boolean converged;

    private boolean oneBadAtATime;

    private int ibad;

    private int maxItBadData;

    private int maxIt;

    private double badDataThreshold;

    public Estimator(PowerSystem powerSystem) {

        this.powerSystem = powerSystem;

        computeDSbusDv();

        computeDSbrDv();

        composeFullHMatrix();

        oneBadAtATime = false;

        maxItBadData = 50;

        maxIt = 100;

        badDataThreshold = 6.25;

//        print();

    }

    public void estimate() {

        converged = false;

        WInvReal = powerSystem.getMeasureSystem().getWInvReal();

        ibad = 1;

        boolean hasBadData;

        BasicMatrix zestReal = computeEstimatedMeasurement(powerSystem);

        BasicMatrix deltzReal = powerSystem.getMeasureSystem().getZmReal().subtract(zestReal);

        BasicMatrix normFReal = deltzReal.transpose().multiply(WInvReal).multiply(deltzReal);

        if (powerSystem.getOption().isVerbose()) {

            System.out.printf("\n it     norm( F )       step size");

            System.out.printf("\n----  --------------  --------------");

            System.out.printf("\n%3d    %10.3f      %10.3f", 0, normFReal.get(0, 0), 0.0);

        }

        while (!converged && ibad < maxItBadData) {

            hasBadData = false;

            BasicMatrix HH = getHH(HFReal);

            BasicMatrix WWInv = getWWInv(WInvReal);

            BasicMatrix ddeltz = getDdeltz(deltzReal);

            BasicMatrix VVs = getVVs(powerSystem.getState());

            BasicMatrix VVsa = phaseOfComplexMatrix(VVs);

            BasicMatrix VVsm = absOfComplexMatrix(VVs);

            int i = 0;

            while (!converged && i++ < maxIt) {

                BasicMatrix b = HH.transpose().multiply(WWInv).multiply(ddeltz);

                BasicMatrix A = HH.transpose().multiply(WWInv).multiply(HH);

                BasicMatrix dx = solveLinear(A, b);

                VVsa = VVsa.add(getDdxa(dx));

                VVsm = VVsm.add(getDdxm(dx));

                refreshState(powerSystem, VVsa, VVsm, powerSystem.getMeasureSystem().getVbusExcludeIds());

                zestReal = computeEstimatedMeasurement(powerSystem);

                deltzReal = powerSystem.getMeasureSystem().getZmReal().subtract(zestReal);

                ddeltz = getDdeltz(deltzReal);

                normFReal = ddeltz.transpose().multiply(WWInv).multiply(ddeltz);

                BasicMatrix dx2 = dx.transpose().multiply(dx);

                if (powerSystem.getOption().isVerbose()) {

                    System.out.printf("\n%3d    %10.3f      %10.3e", i, normFReal.get(0, 0), dx2.get(0, 0));

                }

                if (dx2.get(0, 0).doubleValue() < Constants.ESTIMATOR.TOL) {

                    converged = true;

                    if (powerSystem.getOption().isVerbose()) {

                        System.out.printf("\nState estimator converged in %d iterations.\n", i);

                    }

                }

            }

            if (!converged && powerSystem.getOption().isVerbose()) {

                System.out.printf("\nState estimator did not converged in %d iterations.\n", i);

            }

//            checking bad data
            BasicMatrix WW = computeWW(WWInv);

            BasicMatrix HTWHInv = HH.transpose().multiply(WWInv).multiply(HH).invert();

            BasicMatrix WR = WW.subtract(HH.multiply(HTWHInv).multiply(HH.transpose()).multiply(0.95));

            BasicMatrix WRInvDiagVec = getDiagonalInvVector(WR);

            BasicMatrix rN2 = ddeltz.multiplyElements(ddeltz).multiplyElements(WRInvDiagVec);

            double maxBad = maxInMatrix(rN2);

            List<Integer> baddata = badDataRecognition(rN2, oneBadAtATime);

            if (baddata.size() > 0) {

                hasBadData = true;

                converged = false;

                updateZExclude(baddata, powerSystem.getMeasureSystem().getzExcludeIds());

            }

            if (!hasBadData) {

                converged = true;

                if (powerSystem.getOption().isVerbose()) {

                    logger.info("No remaining bad data, after discarding data {} time(s).", ibad - 1);

                }

            }

            ibad++;

        }

    }

    private BasicMatrix computeWW(BasicMatrix WWInv) {

        basicRealMatrixBuilder = basicRealMatrixFactory.getBuilder((int) WWInv.countRows(), (int) WWInv.countColumns());

        if (WWInv.countColumns() != WWInv.countRows()) {

            logger.error("Not a square matrix can not inverse!");

            return null;

        }

        double a;

        for (int i = 0; i < WWInv.countColumns(); i++) {

            a = WWInv.get(i, i).doubleValue();

            if (a == 0) {

                logger.error("Not a full rank matrix, can not inverse!");

                return null;

            }

            basicRealMatrixBuilder.set(i, i, 1 / a);

        }

        return basicRealMatrixBuilder.build();

    }

    //    WARNNING: all exclude should be sorted in ascending order
    private void updateZExclude(List<Integer> baddata, List<Integer> zExclude) {

        List<Integer> badids = new ArrayList<Integer>();

        int originalIdx;

        for (int i = 0; i < baddata.size(); i++) {

            originalIdx = baddata.get(i);

            for (int j = 0; j < zExclude.size(); j++) {

                if (originalIdx >= zExclude.get(j)) {

                    originalIdx++;

                }

            }

            badids.add(originalIdx);

        }

        for (Integer badidx : badids) {

            int i = 0;

            while (i < zExclude.size()) {

                if (badidx == zExclude.get(i)) {

                    break;

                } else {

                    if (badidx < zExclude.get(i)) {

                        zExclude.add(i, badidx);

                        break;

                    }

                }

                i++;

            }

            if (i == zExclude.size()) {

                zExclude.add(badidx);

            }

        }

        hasDuplicateElement(zExclude, "update " + ibad);

    }

    //    just return the index
    private List<Integer> badDataRecognition(BasicMatrix rn2, boolean oneBadAtATime) {

        List<Integer> ret = new ArrayList<Integer>();

        double threshold, tmp;

        if (oneBadAtATime) {

            threshold = maxInMatrix(rn2);

        } else {

            threshold = badDataThreshold;

        }

        for (int i = 0; i < rn2.countRows(); i++) {

            tmp = rn2.get(i, 0).doubleValue();

            if (tmp >= threshold) {

                ret.add(i);

                if (oneBadAtATime) {

                    break;

                }

            }

        }

        return ret;

    }

    private BasicMatrix getDiagonalInvVector(BasicMatrix in) {

        int r = Math.min((int) in.countRows(), (int) in.countColumns());

        basicRealMatrixBuilder = basicRealMatrixFactory.getBuilder(r, 1);

        for (int i = 0; i < r; i++) {

            basicRealMatrixBuilder.set(i, 0, 1 / in.get(i, i).doubleValue());

        }

        return basicRealMatrixBuilder.build();

    }

    private void refreshState(PowerSystem powerSystem, BasicMatrix vvsa, BasicMatrix vvsm, List<Integer> excludeV) {

        BasicMatrix vs = powerSystem.getState();

        Access2D.Builder<ComplexMatrix> newState = basicComplexMatrixFactory.copy(vs).copyToBuilder();

        int vi = 0;

        int vvi = 0;

        int exci = 0;

        double ph, mag;

        while (vi < vs.countRows()) {

            if (vi != excludeV.get(exci)) {

                ph = vvsa.get(vvi, 0).doubleValue();

                mag = vvsm.get(vvi, 0).doubleValue();

                newState.set(vi, 0, new ComplexNumber(mag * Math.cos(ph), mag * Math.sin(ph)));

                vi++;

                vvi++;

            } else {

                vi++;

                exci++;

            }

        }

        powerSystem.setState(newState.build());

    }

    private BasicMatrix getDdxa(BasicMatrix dx) {

        if (dx.countRows() % 2 != 0 || dx.countColumns() != 1) {

            logger.error("Number of state variable should be even and dx should be single column!");

            return null;

        }

        return dx.getRowsRange(0, (int) dx.countRows() / 2);

    }

    private BasicMatrix getDdxm(BasicMatrix dx) {

        if (dx.countRows() % 2 != 0 || dx.countColumns() != 1) {

            logger.error("Number of state variable should be even and dx should be single column!");

            return null;

        }

        return dx.getRowsRange((int) dx.countRows() / 2, (int) dx.countRows());

    }

    private BasicMatrix getHH(BasicMatrix HF) {

        return excludeRowsColumns(
                HF,
                powerSystem.getMeasureSystem().getzExcludeIds(),
                powerSystem.getMeasureSystem().getStateExcludeIds()
        );

    }

    private BasicMatrix getWWInv(BasicMatrix WInv) {

        return excludeRowsColumns(
                WInv,
                powerSystem.getMeasureSystem().getzExcludeIds(),
                powerSystem.getMeasureSystem().getzExcludeIds()
        );

    }

    private BasicMatrix getDdeltz(BasicMatrix deltz) {

        return excludeRowsColumns(
                deltz,
                powerSystem.getMeasureSystem().getzExcludeIds(),
                null
        );

    }

    private BasicMatrix getVVs(BasicMatrix Vs) {

        return excludeRowsColumns(Vs, powerSystem.getMeasureSystem().getVbusExcludeIds(), null);

    }

    public BasicMatrix computeEstimatedMeasurement(PowerSystem powerSystem) {

        BasicMatrix Vsf = getVft(
                powerSystem.getMpData().getBranchData().getI(),
                powerSystem.getState(),
                powerSystem.getMpData().getBusData().getTOI());

        BasicMatrix Vst = getVft(
                powerSystem.getMpData().getBranchData().getJ(),
                powerSystem.getState(),
                powerSystem.getMpData().getBusData().getTOI());

        BasicMatrix Vs = powerSystem.getState();

        BasicMatrix sfe, ste, sbuse;

        BasicMatrix Vsa, Vsm;

        sfe = Vsf.multiplyElements(powerSystem.getyMatrix().getYf().multiply(Vs).modify(ComplexFunction.CONJUGATE));

        ste = Vst.multiplyElements(powerSystem.getyMatrix().getYt().multiply(Vs).modify(ComplexFunction.CONJUGATE));

        sbuse = Vs.multiplyElements(powerSystem.getyMatrix().getYbus().multiply(Vs).modify(ComplexFunction.CONJUGATE));

        Vsa = phaseOfComplexMatrix(Vs);

        Vsm = absOfComplexMatrix(Vs);

        return toMeasurementVector(
                sfe,
                ste,
                sbuse,
                Vsa,
                Vsm);

    }

    public void print() {

        if (!powerSystem.getOption().isVerbose()) {

            return;

        }

        printWithTitle("dSbus_dVm", dSbDvmCplx);

        printWithTitle("dSbus_dVa", dSbDvaCplx);

        printWithTitle("dSf_dVm", dSfDvmCplx);

        printWithTitle("dSf_dVa", dSfDvaCplx);

        printWithTitle("dSt_dVm", dStDvmCplx);

        printWithTitle("dSt_dVa", dStDvaCplx);

        printWithTitle("SfCplx", SfCplx);

        printWithTitle("StCplx", StCplx);

        printWithTitle("Full H matrix", HFReal);

    }

    public void printWithTitle(String title, BasicMatrix matrix) {

        System.out.print("\n" + title + "\n");
        System.out.print(matrix + "\n");

    }

    /*  this is the measurement jacobi matrix
    *          ang1      ang2 ....   angNb   Vm1   Vm2 ....  VmNb
    *  PF1  dPF1/dAng1
    *  PF2
    *  .
    *  .
    *  PFN
    *  PT1
    *  .
    *  .
    *  PTN
    *  PB1
    *  .
    *  .
    *  PBN
    *  ang1
    *  ang2
    *  .
    *  .
    *  angNb
    *  QF1
    *  .
    *  .
    *  QFN
    *  QT1
    *  .
    *  .
    *  QTN
    *  QB1
    *  .
    *  .
    *  QBN
    *  Vm1
    *  .
    *  .
    *  VmNb
    *
    *
    * */
    private void composeFullHMatrix() {

        int nb = powerSystem.getMpData().getnBus();

        BasicMatrix eyenb = basicRealMatrixFactory.makeEye(nb, nb);

        BasicMatrix zeronb = basicRealMatrixFactory.makeZero(nb, nb);

        BasicMatrix HFRealCol1 = cplxMatrixPart(dSfDvaCplx, true)
                .mergeColumns(cplxMatrixPart(dStDvaCplx, true))
                .mergeColumns(cplxMatrixPart(dSbDvaCplx, true))
                .mergeColumns(eyenb)
                .mergeColumns(cplxMatrixPart(dSfDvaCplx, false))
                .mergeColumns(cplxMatrixPart(dStDvaCplx, false))
                .mergeColumns(cplxMatrixPart(dSbDvaCplx, false))
                .mergeColumns(zeronb);

        BasicMatrix HFRealCol2 = cplxMatrixPart(dSfDvmCplx, true)
                .mergeColumns(cplxMatrixPart(dStDvmCplx, true))
                .mergeColumns(cplxMatrixPart(dSbDvmCplx, true))
                .mergeColumns(zeronb)
                .mergeColumns(cplxMatrixPart(dSfDvmCplx, false))
                .mergeColumns(cplxMatrixPart(dStDvmCplx, false))
                .mergeColumns(cplxMatrixPart(dSbDvmCplx, false))
                .mergeColumns(eyenb);

        HFReal = HFRealCol1.mergeRows(HFRealCol2);

    }

    private void computeDSbrDv() {

        BasicMatrix Yf, Yt, If, It, IfMatrix, ItMatrix, Vf, Vt, VfNorm, VtNorm, VfMatrix, VtMatrix;

        Yf = powerSystem.getyMatrix().getYf();

        Yt = powerSystem.getyMatrix().getYt();

        If = Yf.multiply(VpfCplx);

        It = Yt.multiply(VpfCplx);

        IfMatrix = expandVectorToDiagonalMatrix(If, true);

        ItMatrix = expandVectorToDiagonalMatrix(It, true);

        Vf = getVft(powerSystem.getMpData().getBranchData().getI(), VpfCplx, powerSystem.getMpData().getBusData().getTOI());

        Vt = getVft(powerSystem.getMpData().getBranchData().getJ(), VpfCplx, powerSystem.getMpData().getBusData().getTOI());

        VfNorm = getVft(powerSystem.getMpData().getBranchData().getI(), VpfNormCplx, powerSystem.getMpData().getBusData().getTOI());

        VtNorm = getVft(powerSystem.getMpData().getBranchData().getJ(), VpfNormCplx, powerSystem.getMpData().getBusData().getTOI());

        VfMatrix = expandVectorToDiagonalMatrix(Vf, true);

        VtMatrix = expandVectorToDiagonalMatrix(Vt, true);

        int NBranch = powerSystem.getMpData().getnBranch();

        int NBus = powerSystem.getMpData().getnBus();

        int[] idxBranch = new int[NBranch];

        int[] idxBusFInter = new int[NBranch];

        int[] idxBusTInter = new int[NBranch];

        for (int i = 0; i < NBranch; i++) {

            idxBranch[i] = i;

            idxBusFInter[i] = powerSystem.getMpData().getBusData().getTOI().get(
                    powerSystem.getMpData().getBranchData().getI()[i]) - 1;

            idxBusTInter[i] = powerSystem.getMpData().getBusData().getTOI().get(
                    powerSystem.getMpData().getBranchData().getJ()[i]) - 1;

        }

        BasicMatrix VNbrNbF, VNbrNbT, VNbrNbNormF, VNbrNbNormT;

        VNbrNbF = toSpareMatrix(idxBranch, idxBusFInter, Vf, NBranch, NBus);

        VNbrNbT = toSpareMatrix(idxBranch, idxBusTInter, Vt, NBranch, NBus);

        VNbrNbNormF = toSpareMatrix(idxBranch, idxBusFInter, VfNorm, NBranch, NBus);

        VNbrNbNormT = toSpareMatrix(idxBranch, idxBusTInter, VtNorm, NBranch, NBus);

        dSfDvaCplx = IfMatrix.modify(ComplexFunction.CONJUGATE).multiply(VNbrNbF)
                .subtract(VfMatrix.multiply(Yf.multiply(VpfMatrixCplx).modify(ComplexFunction.CONJUGATE)))
                .multiply((Number) new ComplexNumber(0, 1));

        dSfDvmCplx = VfMatrix.multiply(Yf.multiply(VpfNormMatrixCplx).modify(ComplexFunction.CONJUGATE))
                .add(IfMatrix.modify(ComplexFunction.CONJUGATE).multiply(VNbrNbNormF));

        dStDvaCplx = ItMatrix.modify(ComplexFunction.CONJUGATE).multiply(VNbrNbT)
                .subtract(VtMatrix.multiply(Yt.multiply(VpfMatrixCplx).modify(ComplexFunction.CONJUGATE)))
                .multiply((Number) new ComplexNumber(0, 1));

        dStDvmCplx = VtMatrix.multiply(Yt.multiply(VpfNormMatrixCplx).modify(ComplexFunction.CONJUGATE))
                .add(ItMatrix.modify(ComplexFunction.CONJUGATE).multiply(VNbrNbNormT));

        SfCplx = Vf.multiplyElements(If.modify(ComplexFunction.CONJUGATE));

        StCplx = Vt.multiplyElements(It.modify(ComplexFunction.CONJUGATE));

    }

    private void computeDSbusDv() {

        BasicMatrix IbusMatrix;

        VpfCplx = powerSystem.getPowerFlow().getV();

        VpfNormCplx = computeVnorm(VpfCplx);

        VpfMatrixCplx = expandVectorToDiagonalMatrix(VpfCplx, true);

        IbusMatrix = expandVectorToDiagonalMatrix(powerSystem.getyMatrix().getYbus().multiply(VpfCplx), true);

        VpfNormMatrixCplx = expandVectorToDiagonalMatrix(VpfNormCplx, true);

        dSbDvmCplx = VpfMatrixCplx.multiply(powerSystem.getyMatrix().getYbus().multiply(VpfNormMatrixCplx)
                .modify(ComplexFunction.CONJUGATE))
                .add(IbusMatrix.modify(ComplexFunction.CONJUGATE).multiply(VpfNormMatrixCplx));

        dSbDvaCplx = VpfMatrixCplx.multiply((Number) new ComplexNumber(0, 1)).multiply(
                IbusMatrix.subtract(powerSystem.getyMatrix().getYbus().multiply(VpfMatrixCplx))
                        .modify(ComplexFunction.CONJUGATE)
        );

    }

    //    index should start from 0
    private BasicMatrix toSpareMatrix(int[] i, int[] j, BasicMatrix matrix, int rows, int cols) {

        if (i.length != j.length || matrix.countColumns() > 1) {

            logger.error("Invalid input!");

            return null;

        }

        basicComplexMatrixBuilder = basicComplexMatrixFactory.getBuilder(rows, cols);

        for (int k = 0; k < i.length; k++) {

            basicComplexMatrixBuilder.set(i[k], j[k], matrix.get(k, 0));

        }

        return basicComplexMatrixBuilder.build();

    }

    private BasicMatrix expandVectorToDiagonalMatrix(BasicMatrix vec, boolean complex) {

        if (vec.countColumns() != 1) {

            logger.error("Not a vector!");

            return null;

        }

        int rows = (int) vec.countRows();

        if (!complex) {

            basicRealMatrixBuilder = basicRealMatrixFactory.getBuilder(rows, rows);

            for (int i = 0; i < rows; i++) {

                basicRealMatrixBuilder.set(i, i, vec.get(i, 0));

            }

            return basicRealMatrixBuilder.build();

        } else {

            basicComplexMatrixBuilder = basicComplexMatrixFactory.getBuilder(rows, rows);

            for (int i = 0; i < rows; i++) {

                basicComplexMatrixBuilder.set(i, i, vec.get(i, 0));

            }

            return basicComplexMatrixBuilder.build();

        }

    }

    private BasicMatrix computeVnorm(BasicMatrix v) {

        return v.multiplyElements(absOfComplexMatrix(v).modify(PrimitiveFunction.INVERT));

    }

    //    input are external numbers, will convert to internal numbering
    private BasicMatrix getVft(int[] ijExternal, BasicMatrix VvecCplx, Map<Integer, Integer> TOI) {

        if (VvecCplx.countColumns() != 1) {

            logger.error("Not a vector!");

            return null;

        }

        int n = ijExternal.length;

        basicComplexMatrixBuilder = basicComplexMatrixFactory.getBuilder(n, 1);

        int internalIdx;

        for (int i = 0; i < n; i++) {

            internalIdx = TOI.get(ijExternal[i]) - 1;

            basicComplexMatrixBuilder.set(i, 0, VvecCplx.get(internalIdx, 0));

        }

        return basicComplexMatrixBuilder.build();

    }

    public BasicMatrix getdSbDvaCplx() {
        return dSbDvaCplx;
    }

    public BasicMatrix getdSbDvmCplx() {
        return dSbDvmCplx;
    }

    public BasicMatrix getSfCplx() {
        return SfCplx;
    }

    public BasicMatrix getStCplx() {
        return StCplx;
    }

}
