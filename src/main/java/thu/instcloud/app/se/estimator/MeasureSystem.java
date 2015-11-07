package thu.instcloud.app.se.estimator;

import org.ojalgo.access.Access2D;
import org.ojalgo.function.ComplexFunction;
import org.ojalgo.matrix.BasicMatrix;
import org.ojalgo.matrix.PrimitiveMatrix;
import thu.instcloud.app.se.common.OjMatrixManipulator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static thu.instcloud.app.se.common.Utils.Common.comparator;
import static thu.instcloud.app.se.common.Utils.MatrixExtension.toMeasurementVector;
import static thu.instcloud.app.se.common.Utils.OJ.cplxMatrixPart;
import static thu.instcloud.app.se.common.Utils.OJ.newOneRealBasicMatrix;

/**
 * Created on 2015/11/7.
 */
public class MeasureSystem extends OjMatrixManipulator {

    private PowerSystem powerSystem;

    private BasicMatrix sfCplx;

    private BasicMatrix stCplx;

    private BasicMatrix sbusCplx;

    private BasicMatrix VpfmReal;

    private BasicMatrix VpfaReal;

    private BasicMatrix zTrueReal;

    private BasicMatrix zmReal;

    private Access2D.Builder<PrimitiveMatrix> zmRealBuilder;

    private BasicMatrix sigmaReal;

    private BasicMatrix WInvReal;

    private double fullscale;

    private int nz;

    private List<Integer> excludeIdxSf;

    private List<Integer> excludeIdxSt;

    private List<Integer> VbusExcludeIds;

    private List<Integer> zExcludeIds;

    private List<Integer> stateExcludeIds;

    private Random random;

    public MeasureSystem(PowerSystem powerSystem) {

        this.powerSystem = powerSystem;

        fullscale = 30;

        random = new Random();

        excludeIdxSf = new ArrayList<Integer>();

        excludeIdxSt = new ArrayList<Integer>();

        VbusExcludeIds = new ArrayList<Integer>();

        zExcludeIds = new ArrayList<Integer>();

        stateExcludeIds = new ArrayList<Integer>();

        sfCplx = powerSystem.getEstimator().getSfCplx();

        stCplx = powerSystem.getEstimator().getStCplx();

        sbusCplx = powerSystem.getPowerFlow().getSbus();

        VpfmReal = powerSystem.getPowerFlow().getVm();

        VpfaReal = powerSystem.getPowerFlow().getVa();

        nz = 4 * powerSystem.getMpData().getnBranch() + 4 * powerSystem.getMpData().getnBus();

        zmRealBuilder = basicRealMatrixFactory.getBuilder(nz, 1);

        importTrueMeasurement();

        generateSigma();

        computeWInv();

        computeExcludeIndices();

//        print();

    }

    public void print() {

        if (powerSystem.getOption().isVerbose()) {

            System.out.print("\nReal measurement:\n" + zTrueReal.toString());

        }

    }

    private void importTrueMeasurement() {

        zTrueReal = toMeasurementVector(sfCplx, stCplx, sbusCplx, VpfaReal, VpfmReal);

    }

    private void generateSigma() {

        sigmaReal = cplxMatrixPart(sfCplx.modify(ComplexFunction.ABS), true).multiply(0.02).add(0.0052 * fullscale)
                .mergeColumns(cplxMatrixPart(stCplx.modify(ComplexFunction.ABS), true).multiply(0.02).add(0.0052 * fullscale))
                .mergeColumns(cplxMatrixPart(sbusCplx.modify(ComplexFunction.ABS), true).multiply(0.02).add(0.0052 * fullscale))
                .mergeColumns(newOneRealBasicMatrix(powerSystem.getMpData().getnBus(), 1).multiply(0.2 * Math.PI / 180 * 3))
                .mergeColumns(cplxMatrixPart(sfCplx.modify(ComplexFunction.ABS), true).multiply(0.02).add(0.0052 * fullscale))
                .mergeColumns(cplxMatrixPart(stCplx.modify(ComplexFunction.ABS), true).multiply(0.02).add(0.0052 * fullscale))
                .mergeColumns(cplxMatrixPart(sbusCplx.modify(ComplexFunction.ABS), true).multiply(0.02).add(0.0052 * fullscale))
                .mergeColumns(VpfmReal.multiply(0.02).add(0.0052 * 1.1))
                .divide(3);

    }

    public void measure() {

        for (int i = 0; i < nz; i++) {

            zmRealBuilder.set(i, 0, getMeasureI(i));

        }

        zmReal = zmRealBuilder.build();

    }

    private double getMeasureI(int i) {

        return random.nextGaussian() * sigmaReal.get(i, 0).doubleValue() + zTrueReal.get(i, 0).doubleValue();

    }

    private void computeWInv() {

        basicRealMatrixBuilder = basicRealMatrixFactory.getBuilder(nz, nz);

        double sig;

        for (int i = 0; i < nz; i++) {

            sig = sigmaReal.get(i, 0).doubleValue();

            basicRealMatrixBuilder.set(i, i, 1 / sig / sig);

        }

        WInvReal = basicRealMatrixBuilder.build();

    }

    private void computeExcludeIndices() {

        int nbr = powerSystem.getMpData().getnBranch();

        int nb = powerSystem.getMpData().getnBus();

        int refNumI = powerSystem.getMpData().getBusData().getNrefI();

        int[] I = powerSystem.getMpData().getBranchData().getI();

        int[] J = powerSystem.getMpData().getBranchData().getJ();

        Map<Integer, Integer> TOI = powerSystem.getMpData().getBusData().getTOI();

        for (int i = 0; i < nbr; i++) {

            if (TOI.get(I[i]) == refNumI) {

                excludeIdxSf.add(i);

            }

            if (TOI.get(J[i]) == refNumI) {

                excludeIdxSt.add(i);

            }

        }

        VbusExcludeIds.add(refNumI - 1);

        zExcludeIds.clear();

        stateExcludeIds.clear();

        int exIdx;

        for (int i = 0; i < excludeIdxSf.size(); i++) {

            exIdx = excludeIdxSf.get(i);

//            Pf
            zExcludeIds.add(exIdx);

//            Qf
            zExcludeIds.add(exIdx + 2 * (nbr + nb));

        }

        for (int i = 0; i < excludeIdxSt.size(); i++) {

            exIdx = excludeIdxSt.get(i);

//            Pt
            zExcludeIds.add(exIdx + nbr);

//            Qt
            zExcludeIds.add(exIdx + 3 * nbr + 2 * nb);

        }

        for (int i = 0; i < VbusExcludeIds.size(); i++) {

            exIdx = VbusExcludeIds.get(i);

//            Pbus
            zExcludeIds.add(exIdx + 2 * nbr);

//            Qbus
            zExcludeIds.add(exIdx + 4 * nbr + 2 * nb);

//            Va
            zExcludeIds.add(exIdx + 2 * nbr + nb);

//            Vm
            zExcludeIds.add(exIdx + 4 * nbr + 3 * nb);

            stateExcludeIds.add(exIdx);

            stateExcludeIds.add(nb + exIdx);

        }

        zExcludeIds.sort(comparator);

        stateExcludeIds.sort(comparator);

        VbusExcludeIds.sort(comparator);

    }

    public int getNz() {
        return nz;
    }

    public List<Integer> getzExcludeIds() {
        return zExcludeIds;
    }

    public List<Integer> getStateExcludeIds() {
        return stateExcludeIds;
    }

    public BasicMatrix getWInvReal() {
        return WInvReal;
    }

    public BasicMatrix getZmReal() {
        return zmReal;
    }

    public BasicMatrix getzTrueReal() {
        return zTrueReal;
    }

    public List<Integer> getVbusExcludeIds() {
        return VbusExcludeIds;
    }

}
