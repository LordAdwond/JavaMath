package JavaMath;

import javax.management.InvalidAttributeValueException;

public class linAlg
{
    public static double[] linSpace(double a, double b, int count) throws InvalidAttributeValueException
    {
        if(a>b)
        {
            throw new InvalidAttributeValueException("Parameter 'a' must be equal or greater than parameter 'b'");
        }

        double[] res = new double[count];
        double h = (b-a)/count;

        for(int i=0; i<count; i++)
        {
            res[i] = a+i*h;
        }

        return res;
    }

    public static Matrix solve(Matrix A, Vector b) throws InvalidAttributeValueException
    {
        return Matrix.matMul(Matrix.pinv(A), b.toMatrix(VectorAsMatrixType.COLUMN));
    }

    public static double[][] solve(double[][] A, double[] b) throws InvalidAttributeValueException
    {
        Matrix newA = new Matrix(A);
        Vector newB = new Vector(b);
        return Matrix.matMul(Matrix.pinv(newA), newB.toMatrix(VectorAsMatrixType.COLUMN)).getArray();
    }
}
