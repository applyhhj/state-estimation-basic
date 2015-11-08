package thu.instcloud.app.se.mpdata;

import org.ojalgo.matrix.BasicMatrix;
import org.ojalgo.scalar.ComplexNumber;
import thu.instcloud.app.se.common.OjMatrixManipulator;

import java.util.ArrayList;
import java.util.List;

/**
 * Created on 2015/11/6.
 */
public class YMatrix extends OjMatrixManipulator {

    private List<ComplexNumber> Ytt;

    private List<ComplexNumber> Yff;

    private List<ComplexNumber> Ytf;

    private List<ComplexNumber> Yft;

    private BasicMatrix cf;

    private BasicMatrix ct;

    private BasicMatrix YSh;

//    private BasicMatrix YShB;

    private BasicMatrix Yf;

//    private BasicMatrix YfB;

    private BasicMatrix Yt;

//    private BasicMatrix YtB;

    private BasicMatrix Ybus;

//    private Matrix YB;

//    private ComplexMatrix Ybus;

    private MPData mpData;

    public YMatrix(MPData mpData) {

        this.mpData = mpData;

        Ytt = new ArrayList<ComplexNumber>();

        Yff = new ArrayList<ComplexNumber>();

        Ytf = new ArrayList<ComplexNumber>();

        Yft = new ArrayList<ComplexNumber>();

        computeYMatrix();

//        print();

    }

    private void computeYMatrix() {

        double Gtt, Btt, Gff, Bff, Gft, Bft, Gtf, Btf, Gs, Bs, r, x, zm2, Bc, t, tsh;

        for (int i = 0; i < mpData.getBranchData().getN(); i++) {

            r = mpData.getBranchData().getR()[i];

            x = mpData.getBranchData().getX()[i];

            Bc = mpData.getBranchData().getB()[i];

            t = mpData.getBranchData().getRatio()[i];

            tsh = mpData.getBranchData().getAngle()[i];

            if (t <= 0) {

                t = 1;

            }

            zm2 = r * r + x * x;

            Gs = r / zm2;

            Bs = -x / zm2;

            Gtt = Gs;

            Btt = Bc / 2 + Bs;

            Gff = Gtt / t / t;

            Bff = Btt / t / t;

            if (tsh == 0) {

                Gtf = Gft = -Gs / t;

                Btf = Bft = -Bs / t;

            } else {

                Gtf = -(Gs * Math.cos(tsh) + Bs * Math.sin(tsh)) / t;

                Btf = (Gs * Math.sin(tsh) - Bs * Math.cos(tsh)) / t;

                Gft = (Bs * Math.sin(tsh) - Gs * Math.cos(tsh)) / t;

                Bft = -(Gs * Math.sin(tsh) + Bs * Math.cos(tsh)) / t;

            }

            Ytt.add(new ComplexNumber(Gtt, Btt));

            Yff.add(new ComplexNumber(Gff, Bff));

            Yft.add(new ComplexNumber(Gft, Bft));

            Ytf.add(new ComplexNumber(Gtf, Btf));

        }

        getConnectionMatrix();

        getYfYt();

        getYSparseSh();

        Ybus = cf.transpose().multiply(Yf).add(ct.transpose().multiply(Yt)).add(YSh);

    }

    private void print() {

        System.out.print("cf\n" + cf.toString() + "\n");

        System.out.print("ct\n" + ct.toString() + "\n");

        System.out.print("YSh\n" + YSh.toString() + "\n");

//        System.out.print("YShB\n" + YShB.toString() + "\n");

        System.out.print("Yf\n" + Yf.toString() + "\n");

//        System.out.print("YfB\n" + YfB.toString() + "\n");

        System.out.print("Yt\n" + Yt.toString() + "\n");

//        System.out.print("YtB\n" + YtB.toString() + "\n");

        System.out.print("Ybus\n" + Ybus.toString() + "\n");

//        System.out.print("YB\n" + YB.toString() + "\n");

    }

    private void getYSparseSh() {

        int nbu = mpData.getBusData().getN();

        int idx;

        basicComplexMatrixBuilder = basicComplexMatrixFactory.getBuilder(nbu, nbu);

        for (int i = 0; i < nbu; i++) {

//            this is the index, however i is index should convert to internal bus number
            idx = mpData.getBusData().getTOA().get(mpData.getBusData().getTIO().get(i + 1));

            ComplexNumber ysh = new ComplexNumber(
                    mpData.getBusData().getGs()[idx] / mpData.getSbase(),
                    mpData.getBusData().getBs()[idx] / mpData.getSbase()
            );

            basicComplexMatrixBuilder.set(i, i, ysh);

        }

        YSh = basicComplexMatrixBuilder.build();

    }

    private void getYfYt() {

        int nbr = mpData.getBranchData().getN();

        int nbu = mpData.getBusData().getN();

        int idxi, idxj;

        basicComplexMatrixBuilder = basicComplexMatrixFactory.getBuilder(nbr, nbu);

        for (int i = 0; i < nbr; i++) {

//            convert to internal number
            idxi = mpData.getBusData().getTOI().get(mpData.getBranchData().getI()[i]) - 1;

            basicComplexMatrixBuilder.set(i, idxi, Yff.get(i));

            idxj = mpData.getBusData().getTOI().get(mpData.getBranchData().getJ()[i]) - 1;

            basicComplexMatrixBuilder.set(i, idxj, Yft.get(i));

        }

        Yf = basicComplexMatrixBuilder.build();

        basicComplexMatrixBuilder = basicComplexMatrixFactory.getBuilder(nbr, nbu);

        for (int i = 0; i < nbr; i++) {

//            convert to internal number
            idxi = mpData.getBusData().getTOI().get(mpData.getBranchData().getI()[i]) - 1;

            basicComplexMatrixBuilder.set(i, idxi, Ytf.get(i));

            idxj = mpData.getBusData().getTOI().get(mpData.getBranchData().getJ()[i]) - 1;

            basicComplexMatrixBuilder.set(i, idxj, Ytt.get(i));

        }

        Yt = basicComplexMatrixBuilder.build();

    }

    private void getConnectionMatrix() {

        int nbr = mpData.getBranchData().getN();

        int nbu = mpData.getBusData().getN();

        basicComplexMatrixBuilder = basicComplexMatrixFactory.getBuilder(nbr, nbu);

        for (int i = 0; i < nbr; i++) {

//            convert external bus number to internal bus number, convert to index
            basicComplexMatrixBuilder.set(i, mpData.getBusData().getTOI().get(mpData.getBranchData().getI()[i]) - 1,
                    new ComplexNumber(1, 0));

        }

        cf = basicComplexMatrixBuilder.build();

        basicComplexMatrixBuilder = basicComplexMatrixFactory.getBuilder(nbr, nbu);

        for (int i = 0; i < nbr; i++) {

//            convert external bus number to internal bus number, convert to index
            basicComplexMatrixBuilder.set(i, mpData.getBusData().getTOI().get(mpData.getBranchData().getJ()[i]) - 1,
                    new ComplexNumber(1, 0));

        }

        ct = basicComplexMatrixBuilder.build();

    }

    public BasicMatrix getYbus() {
        return Ybus;
    }

    public BasicMatrix getYt() {
        return Yt;
    }

    public BasicMatrix getYf() {
        return Yf;
    }
}
