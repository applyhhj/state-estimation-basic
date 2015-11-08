package thu.instcloud.app.se.mpdata;

import org.ojalgo.access.Access2D;
import org.ojalgo.function.ComplexFunction;
import org.ojalgo.matrix.BasicMatrix;
import org.ojalgo.matrix.ComplexMatrix;
import org.ojalgo.matrix.PrimitiveMatrix;
import org.ojalgo.scalar.ComplexNumber;
import thu.instcloud.app.se.common.OjMatrixManipulator;

/**
 * Created on 2015/11/6.
 */
public class PowerFlow extends OjMatrixManipulator {

    private MPData mpData;

    private YMatrix yMatrix;
//
//    private BasicMatrix Vr;
//
//    private BasicMatrix Vi;

    private BasicMatrix PF;

    private BasicMatrix QF;

    private BasicMatrix PT;

    private BasicMatrix QT;

//    private BasicMatrix SbusP;
//
//    private BasicMatrix SbusQ;

    private BasicMatrix V;

    private BasicMatrix Vm;

    private BasicMatrix Va;

    private BasicMatrix Sbus;

    public PowerFlow(MPData mpData, YMatrix yMatrix) {

        this.mpData = mpData;

        this.yMatrix = yMatrix;

        importV();

        importPQ();

        computeSbus();

//        print();

    }

    //    internal bus numbering
    private void computeSbus() {

        Sbus = V.multiplyElements(yMatrix.getYbus().multiply(V).modify(ComplexFunction.CONJUGATE));

    }

    //    internal bus numbering
    private void importV() {

        Access2D.Builder<PrimitiveMatrix> VmBuilder = basicRealMatrixFactory.getBuilder(mpData.getnBus(), 1);

        Access2D.Builder<PrimitiveMatrix> VaBuilder = basicRealMatrixFactory.getBuilder(mpData.getnBus(), 1);

        Access2D.Builder<ComplexMatrix> VBuilder = basicComplexMatrixFactory.getBuilder(mpData.getnBus(), 1);

        int idx;

        double vm, va;

        for (int i = 0; i < mpData.getBusData().getN(); i++) {

            idx = mpData.getBusData().getTOA().get(mpData.getBusData().getTIO().get(i + 1));

            vm = mpData.getBusData().getVoltage()[idx];

            va = mpData.getBusData().getAngle()[idx];

//            convert to internal bus number
            VmBuilder.set(i, 0, vm);

            VaBuilder.set(i, 0, va);

            VBuilder.set(i, 0, new ComplexNumber(vm * Math.cos(va), vm * Math.sin(va)));

        }

        Va = VaBuilder.build();

        Vm = VmBuilder.build();

        V = VBuilder.build();

    }

    //    for comparison
    private void importPQ() {

        for (int j = 0; j < 4; j++) {

            basicRealMatrixBuilder = basicRealMatrixFactory.getBuilder(mpData.getnBranch(), 1);

            for (int i = 0; i < mpData.getnBranch(); i++) {

                switch (j) {

                    case 0:
                        basicRealMatrixBuilder.set(i, 0, mpData.getBranchData().getPF()[i]);
                        break;

                    case 1:
                        basicRealMatrixBuilder.set(i, 0, mpData.getBranchData().getQF()[i]);
                        break;

                    case 2:
                        basicRealMatrixBuilder.set(i, 0, mpData.getBranchData().getPT()[i]);
                        break;

                    case 3:
                        basicRealMatrixBuilder.set(i, 0, mpData.getBranchData().getQT()[i]);
                        break;

                }

            }

            switch (j) {

                case 0:
                    PF = basicRealMatrixBuilder.build();
                    break;

                case 1:
                    QF = basicRealMatrixBuilder.build();
                    break;

                case 2:
                    PT = basicRealMatrixBuilder.build();
                    break;

                case 3:
                    QT = basicRealMatrixBuilder.build();
                    break;

            }

        }


    }

    private void print() {

        System.out.print("V\n" + V.toString() + "\n");

//        System.out.print("Vi\n" + Vi.toString() + "\n");

        System.out.print("Sbus\n" + Sbus.toString() + "\n");

//        System.out.print("SbusQ\n" + SbusQ.toString() + "\n");

    }

    public BasicMatrix getV() {
        return V;
    }

    public BasicMatrix getSbus() {
        return Sbus;
    }

    public BasicMatrix getVm() {
        return Vm;
    }

    public BasicMatrix getVa() {
        return Va;
    }
}
