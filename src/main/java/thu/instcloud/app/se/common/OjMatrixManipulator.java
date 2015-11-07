package thu.instcloud.app.se.common;

import org.ojalgo.access.Access2D;
import org.ojalgo.matrix.BasicMatrix;
import org.ojalgo.matrix.ComplexMatrix;
import org.ojalgo.matrix.PrimitiveMatrix;

/**
 * Created on 2015/11/6.
 */
public class OjMatrixManipulator {

    protected final BasicMatrix.Factory<PrimitiveMatrix> basicRealMatrixFactory = PrimitiveMatrix.FACTORY;

    protected final BasicMatrix.Factory<ComplexMatrix> basicComplexMatrixFactory = ComplexMatrix.FACTORY;

    protected Access2D.Builder<ComplexMatrix> basicComplexMatrixBuilder;

    protected Access2D.Builder<PrimitiveMatrix> basicRealMatrixBuilder;

}
