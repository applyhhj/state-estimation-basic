package thu.instcloud.app.se.common;

import org.ojalgo.access.Access2D;
import org.ojalgo.function.ComplexFunction;
import org.ojalgo.matrix.BasicMatrix;
import org.ojalgo.matrix.PrimitiveMatrix;
import org.ojalgo.matrix.decomposition.DecompositionStore;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.matrix.store.PhysicalStore;
import org.ojalgo.matrix.task.SolverTask;
import org.ojalgo.matrix.task.TaskException;
import org.ojalgo.scalar.ComplexNumber;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import static thu.instcloud.app.se.common.Utils.OJ.cplxMatrixPart;

/**
 * Created on 2015/11/6.
 */
public class Utils {

    public static class Common {

        public static IntComparator comparator = new IntComparator();

        public static boolean isLinux() {

            return System.getProperty("os.name").toLowerCase().indexOf("linux") >= 0;

        }

        public static List<String> readStringFromFile(String FILE_IN) {

            List<String> ret = new ArrayList<String>();

            File file = new File(FILE_IN);

            try {

                FileInputStream is = new FileInputStream(file);

                InputStreamReader isr = new InputStreamReader(is);

                BufferedReader in = new BufferedReader(isr);

                String line = null;

                while ((line = in.readLine()) != null) {
                    ret.add(line.toString());

                }

                in.close();

                is.close();

            } catch (Exception e) {

                // TODO Auto-generated catch block

                e.printStackTrace();

            }

            return ret;

        }

        public static class IntComparator implements Comparator<Integer> {

            public int compare(Integer o1, Integer o2) {
                if (o1 == o2) {

                    return 0;

                } else if (o1 > o2) {

                    return 1;

                } else {

                    return -1;

                }
            }

        }

    }

    public static class MatrixExtension {

        public static double maxInMatrix(BasicMatrix matrix) {

            double ele, ret = Double.MIN_VALUE;

            Iterator<Number> iterator = matrix.iterator();

            while (iterator.hasNext()) {

                ele = iterator.next().doubleValue();

                if (ret < ele) {

                    ret = ele;

                }

            }

            return ret;

        }

        public static BasicMatrix solveLinear(BasicMatrix matA, BasicMatrix matB) {

            MatrixStore<Double> matAStore = matA.toPrimitiveStore();

            MatrixStore<Double> matBStore = matB.toPrimitiveStore();

            MatrixStore<Double> result = null;

            SolverTask<Double> tmpSolver = SolverTask.PRIMITIVE.make(matAStore, matBStore, false);

            final DecompositionStore<Double> tmpAlloc = tmpSolver.preallocate(matAStore, matBStore);

            try {

                result = tmpSolver.solve(matAStore, matBStore, tmpAlloc);

            } catch (TaskException ex) {

                ex.printStackTrace();

            }

            BasicMatrix.Factory<PrimitiveMatrix> basicRealMatrix2dFactory = PrimitiveMatrix.FACTORY;

            return result == null ? null : basicRealMatrix2dFactory.copy(result);

        }

        public static BasicMatrix excludeRowsColumns(BasicMatrix thisMatrix, List<Integer> excRows, List<Integer> excCols) {

            BasicMatrix ret = null;

            int j;

            if (excRows != null && excRows.size() > 0) {

                j = 0;

                int rn = (int) thisMatrix.countRows();

                int[] selrs = new int[rn - excRows.size()];

                for (int i = 0; i < rn; i++) {

                    if (excRows.contains(i)) {

                        continue;

                    } else {

                        selrs[j++] = i;

                    }

                }

                ret = thisMatrix.selectRows(selrs);

            }

            if (excCols != null && excCols.size() > 0) {

                int cn = (int) thisMatrix.countColumns();

                int[] selcs = new int[cn - excCols.size()];

                j = 0;

                for (int i = 0; i < cn; i++) {

                    if (excCols.contains(i)) {

                        continue;

                    } else {

                        selcs[j++] = i;

                    }

                }

                if (ret != null) {

                    ret = ret.selectColumns(selcs);

                } else {

                    ret = thisMatrix.selectColumns(selcs);

                }

            }

            if (ret == null) {

                ret = thisMatrix.copyToBuilder().build();

            }

            return ret;

        }

        public static BasicMatrix toMeasurementVector(BasicMatrix sf_, BasicMatrix st_, BasicMatrix sbus_,
                                                      BasicMatrix Va_, BasicMatrix Vm_) {

            return cplxMatrixPart(sf_, true)
                    .mergeColumns(cplxMatrixPart(st_, true))
                    .mergeColumns(cplxMatrixPart(sbus_, true))
                    .mergeColumns(Va_)
                    .mergeColumns(cplxMatrixPart(sf_, false))
                    .mergeColumns(cplxMatrixPart(st_, false))
                    .mergeColumns(cplxMatrixPart(sbus_, false))
                    .mergeColumns(Vm_);

        }

    }

    public static class OJ {

        private static final BasicMatrix.Factory<PrimitiveMatrix> basicRealMatrix2dFactory = PrimitiveMatrix.FACTORY;

        private static Access2D.Builder<PrimitiveMatrix> basicRealMatrix2dBuilder;

        public static void printOjMatrix(BasicMatrix matrix) {

            System.out.print("\n");

            for (int i = 0; i < matrix.countRows(); i++) {

                for (int j = 0; j < matrix.countColumns(); j++) {

                    System.out.print(matrix.get(i, j) + "    ");

                }

                System.out.print("\n");

            }

        }

        public static BasicMatrix newRealBasicMatrix(double[][] array2d) {

            int rows = array2d.length;

            int cols = array2d[0].length;

            basicRealMatrix2dBuilder = basicRealMatrix2dFactory.getBuilder(rows, cols);

            for (int i = 0; i < rows; i++) {

                for (int j = 0; j < cols; j++) {

                    basicRealMatrix2dBuilder.set(i, j, array2d[i][j]);

                }

            }

            return basicRealMatrix2dBuilder.build();

        }

        public static BasicMatrix newZeroRealBasicMatrix(int rows, int cols) {

            return basicRealMatrix2dFactory.makeZero(rows, cols);

        }

        public static BasicMatrix newOneRealBasicMatrix(int rows, int cols) {

            basicRealMatrix2dBuilder = basicRealMatrix2dFactory.getBuilder(rows, cols);

            basicRealMatrix2dBuilder.fillAll(1);

            return basicRealMatrix2dBuilder.build();

        }

        public static BasicMatrix absOfComplexMatrix(final BasicMatrix matrix) {

            return cplxMatrixPart(matrix.modify(ComplexFunction.ABS), true);

        }

        public static BasicMatrix phaseOfComplexMatrix(final BasicMatrix matrix) {

            PhysicalStore<ComplexNumber> matrixCplx = matrix.toComplexStore();

            int rows = (int) matrix.countRows();

            int cols = (int) matrix.countColumns();

            basicRealMatrix2dBuilder = basicRealMatrix2dFactory.getBuilder(rows, cols);

            for (int i = 0; i < rows; i++) {

                for (int j = 0; j < cols; j++) {

                    basicRealMatrix2dBuilder.set(i, j, matrixCplx.get(i, j).phase());

                }

            }

            return basicRealMatrix2dBuilder.build();

        }

        public static BasicMatrix cplxMatrixPart(final BasicMatrix matrix, boolean real) {

            final Access2D<ComplexNumber> arg = matrix.toComplexStore();

            int rows = (int) arg.countRows();

            int cols = (int) arg.countColumns();

            basicRealMatrix2dBuilder = basicRealMatrix2dFactory.getBuilder(rows, cols);

            for (int i = 0; i < rows; i++) {

                for (int j = 0; j < cols; j++) {

                    if (real) {

                        basicRealMatrix2dBuilder.set(i, j, arg.get(i, j).getReal());

                    } else {

                        basicRealMatrix2dBuilder.set(i, j, arg.get(i, j).getImaginary());

                    }

                }

            }

            return basicRealMatrix2dBuilder.build();

        }

    }

}
