package JavaMath;

import javax.management.InvalidAttributeValueException;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.Function;
import java.util.stream.IntStream;

public class Matrix implements Sequence
{
    double[][] array;
    int i, j;

    public Matrix(int rows, int columns) throws InvalidAttributeValueException {
        if(rows<=0 || columns<=0)
        {
            throw new InvalidAttributeValueException("Rows number and columns number of vector must be natural");
        }

        this.array = new double[rows][columns];
    }

    public Matrix(double[][] arr) throws InvalidAttributeValueException
    {
        if(arr.length<1)
        {
            throw new InvalidAttributeValueException("This array can't be interpreted as matrix");
        }

        if(Matrix.isArrayCorrect(arr))
        {
            // Використовуйте глибоке копіювання, щоб уникнути зовнішніх змін масиву
            this.array = new double[arr.length][];
            for(int i = 0; i < arr.length; i++) {
                this.array[i] = arr[i].clone();
            }
        }
        else
        {
            throw new InvalidAttributeValueException("Rows in this array have different length");
        }
    }

    public void apply(Function<Double, Double> f)
    {
        Arrays.stream(this.array).parallel().forEach(x -> {
            for(i=0; i<x.length; i++)
            {
                f.apply(x[i]);
            }
        });
    }

    public double[][] getArray()
    {
        return this.array;
    }

    public void setArray(double[][] array) throws InvalidAttributeValueException
    {
        if(Matrix.isArrayCorrect(array))
        {
            this.array = array;
        }
        else
        {
            throw new InvalidAttributeValueException("Rows in this array have different length");
        }
    }

    public void add(Matrix other) throws InvalidAttributeValueException
    {
        if(this.array.length!=other.array.length || this.array[0].length!=other.array[0].length)
        {
            throw new InvalidAttributeValueException("Matrices must have same number of row and number of columns");
        }

        int i, j;
        for(i=0; i<this.array.length; i++)
        {
            for(j=0; j<this.array[0].length; j++)
            {
                this.array[i][j] += other.array[i][j];
            }
        }
    }

    public void subtract(Matrix other) throws InvalidAttributeValueException
    {
        if(this.array.length!=other.array.length || this.array[0].length!=other.array[0].length)
        {
            throw new InvalidAttributeValueException("Matrices must have same number of row and number of columns");
        }

        int i, j;
        for(i=0; i<this.array.length; i++)
        {
            for(j=0; j<this.array[0].length; j++)
            {
                this.array[i][j] -= other.array[i][j];
            }
        }
    }

    public void multiply(double a)
    {
        int i, j;
        for(i=0; i<this.array.length; i++)
        {
            for(j=0; j<this.array[0].length; j++)
            {
                this.array[i][j] *= a;
            }
        }
    }

    public void divide(double a) throws InvalidAttributeValueException
    {
        if(a==0)
        {
            throw new InvalidAttributeValueException("Division by zero");
        }

        int i, j;
        for(i=0; i<this.array.length; i++)
        {
            for(j=0; j<this.array[0].length; j++)
            {
                this.array[i][j] /= a;
            }
        }
    }

    @Override
    public double norm() {
        double L = 0;

        for(i=0; i<this.array.length; i++)
        {
            for(j=0; j<this.array[0].length; j++)
            {
                L += Math.pow(this.array[i][j], 2);
            }
        }

        return Math.sqrt(L);
    }

    @Override
    public double minElement() {
        double res = this.array[0][0];

        for(i=0; i<this.array.length; i++)
        {
            for(j=0; j<this.array[i].length; j++)
            {
                res = Math.min(res, this.array[i][j]);
            }
        }

        return res;
    }

    @Override
    public double maxElement() {
        double res = this.array[0][0];

        for(i=0; i<this.array.length; i++)
        {
            for(j=0; j<this.array[i].length; j++)
            {
                res = Math.max(res, this.array[i][j]);
            }
        }

        return res;
    }

    public static Matrix matMul(Matrix A, Matrix B) throws InvalidAttributeValueException
    {
        if(A.array[0].length!= B.array.length)
        {
            throw new InvalidAttributeValueException("Matrices can't be multiplied");
        }

        double[][] C = new double[A.array.length][B.array[0].length];
        IntStream.range(0, C.length).parallel().forEach(i -> {
            int j, k;
            for(j=0; j<B.array[0].length; j++)
            {
                C[i][j] = 0;
            }
            for(j=0; j<B.array[0].length; j++)
            {
                for(k=0; k<A.array[0].length; k++)
                {
                    C[i][j] += A.array[i][k]*B.array[k][j];
                }
            }
        });

        return new Matrix(C);
    }

    public static Matrix matMul(double alpha, Matrix A) throws InvalidAttributeValueException
    {
        double[][] newArray = new double[A.array.length][A.array[0].length];
        IntStream.range(0, A.array.length).parallel().forEach(i -> {
            int j;
            for(j=0; j<A.array[0].length; j++)
            {
                newArray[i][j] = alpha*A.array[i][j];
            }
        });

        return new Matrix(newArray);
    }

    public static Matrix add(Matrix A, Matrix B) throws InvalidAttributeValueException
    {
        if(A.array.length!=B.array.length || A.array[0].length!=B.array[0].length)
        {
            throw new InvalidAttributeValueException("Matrices must have same number of row and number of columns");
        }

        double[][] newArray = new double[A.array.length][A.array[0].length];
        int i, j;
        for(i=0; i<A.array.length; i++)
        {
            for(j=0; j<A.array[0].length; j++)
            {
                newArray[i][j] = A.array[i][j] + B.array[i][j];
            }
        }

        return new Matrix(newArray);
    }

    public static Matrix subtract(Matrix A, Matrix B) throws InvalidAttributeValueException
    {
        if(A.array.length!=B.array.length || A.array[0].length!=B.array[0].length)
        {
            throw new InvalidAttributeValueException("Matrices must have same number of row and number of columns");
        }

        double[][] newArray = new double[A.array.length][A.array[0].length];
        int i, j;
        for(i=0; i<A.array.length; i++)
        {
            for(j=0; j<A.array[0].length; j++)
            {
                newArray[i][j] = A.array[i][j] - B.array[i][j];
            }
        }

        return new Matrix(newArray);
    }

    public static Matrix divide(Matrix A, double alpha) throws InvalidAttributeValueException
    {
        if(alpha==0)
        {
            throw new InvalidAttributeValueException("Division by zero is not allowed");
        }

        double[][] newArray = new double[A.array.length][A.array[0].length];
        IntStream.range(0, A.array.length).parallel().forEach(i -> {
            int j;
            for(j=0; j<A.array[0].length; j++)
            {
                newArray[i][j] = A.array[i][j]/alpha;
            }
        });

        return new Matrix(newArray);
    }

    public static double det(Matrix A) throws InvalidAttributeValueException
    {
        if(!Matrix.isMatrixSquare(A))
        {
            throw new InvalidAttributeValueException("Matrix isn't square");
        }
        int N = A.array.length;
        if(N==1)
        {
            return A.array[0][0];
        }

        double determinant = 0;

        // Розкладання по першому рядку
        for(int k=0; k<N; k++)
        {
            determinant += A.array[0][k] * cofactor(A, 0, k);
        }

        return determinant;
    }

    // checks
    private static boolean isArrayCorrect(double[][] arr)
    {
        int firstRowLen = arr[0].length;
        AtomicBoolean areLensOfAllRowsSame = new AtomicBoolean(true);
        IntStream.range(0, arr.length).forEach(i -> areLensOfAllRowsSame.set(areLensOfAllRowsSame.get() && arr[i].length==firstRowLen));

        return areLensOfAllRowsSame.get() && arr[0].length>0;
    }

    public static boolean isMatrixSquare(Matrix A)
    {
        return A.array.length == A.array[0].length;
    }

    public static Matrix E(int n) throws InvalidAttributeValueException
    {
        if(n<1)
        {
            throw new InvalidAttributeValueException("A rank of unit matrix must be natural");
        }

        double[][] array = new double[n][n];
        int i, j;

        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                array[i][j] = i==j ? 1 : 0;
            }
        }

        return new Matrix(array);
    }

    // matrix manipulations
    public static Matrix transpond(Matrix A) throws InvalidAttributeValueException
    {
        double[][] newArray = new double[A.array[0].length][A.array.length];
        int i, j;

        for(i=0; i<A.array[0].length; i++)
        {
            for(j=0; j<A.array.length; j++)
            {
                newArray[i][j] = A.array[j][i];
            }
        }

        return new Matrix(newArray);
    }

    public static Matrix inv(Matrix A) throws InvalidAttributeValueException
    {
        if(!Matrix.isMatrixSquare(A))
        {
            throw new InvalidAttributeValueException("Matrix must be square");
        }

        double determinant = Matrix.det(A);
        if(Math.abs(determinant) < 1e-9)
        {
            throw new InvalidAttributeValueException("Matrix can't be inversed (det equals 0)");
        }

        int N = A.array.length;
        double[][] adjArray = new double[N][N]; // Матриця алгебраїчних доповнень C

        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                adjArray[i][j] = Matrix.cofactor(A, i, j);
            }
        }
        Matrix cofactorMatrix = new Matrix(adjArray);
        Matrix adjA = Matrix.transpond(cofactorMatrix);

        return Matrix.divide(adjA, determinant);
    }

    public static Matrix pinv(Matrix A) throws InvalidAttributeValueException {
        int m = A.array.length;
        int n = A.array[0].length;
        if (m == n) {
            return Matrix.inv(A);
        }

        Matrix AT = Matrix.transpond(A);

        if (m > n) {
            Matrix ATA = Matrix.matMul(AT, A); // n x n

            double detATA = Matrix.det(ATA);
            if (Math.abs(detATA) < 1e-9) {
                throw new InvalidAttributeValueException("Cannot compute pseudoinverse: Matrix (A^T * A) is singular (det is zero).");
            }

            Matrix ATA_inv = Matrix.inv(ATA); // (n x n)^(-1)
            return Matrix.matMul(ATA_inv, AT); // n x m
        }

        if (m < n) {
            // Формула: A+ = A^T * (A * A^T)^(-1)
            Matrix AAT = Matrix.matMul(A, AT); // m x m

            double detAAT = Matrix.det(AAT);
            if (Math.abs(detAAT) < 1e-9) {
                throw new InvalidAttributeValueException("Cannot compute pseudoinverse: Matrix (A * A^T) is singular (det is zero).");
            }

            Matrix AAT_inv = Matrix.inv(AAT); // (m x m)^(-1)
            return Matrix.matMul(AT, AAT_inv); // n x m
        }

        throw new InvalidAttributeValueException("Unknown error during pseudoinverse calculation.");
    }

    public static Matrix pow(Matrix A, int n) throws InvalidAttributeValueException
    {
        if(!Matrix.isMatrixSquare(A))
        {
            throw new InvalidAttributeValueException("Matrix must be square");
        }
        if(n<=0)
        {
            throw new InvalidAttributeValueException("Degree must be natural");
        }

        return n==1 ? A : Matrix.matMul(A, pow(A, n-1));
    }

    public static Matrix exp(Matrix A) throws InvalidAttributeValueException
    {
        if(!Matrix.isMatrixSquare(A))
        {
            throw new InvalidAttributeValueException("Matrix must be square");
        }

        Matrix res = Matrix.E(A.array.length);
        for(int i=0; i<10; i++)
        {
            res = Matrix.add(res, Matrix.divide(Matrix.pow(A, i), Combinatorics.factorial(i)));
        }
        return new Matrix(res.getArray());
    }

    public static Matrix sin(Matrix A) throws InvalidAttributeValueException
    {
        if(!Matrix.isMatrixSquare(A))
        {
            throw new InvalidAttributeValueException("Matrix must be square");
        }

        Matrix res = Matrix.E(A.array.length);
        for(int i=0; i<10; i++)
        {
            res = Matrix.add(res, Matrix.divide(Matrix.pow(A, 2*i+1), Math.pow(-1, i)*Combinatorics.factorial(2*i+1)));
        }
        return new Matrix(res.getArray());
    }

    public static Matrix cos(Matrix A) throws InvalidAttributeValueException
    {
        if(!Matrix.isMatrixSquare(A))
        {
            throw new InvalidAttributeValueException("Matrix must be square");
        }

        Matrix res = Matrix.E(A.array.length);
        for(int i=0; i<10; i++)
        {
            res = Matrix.add(res, Matrix.divide(Matrix.pow(A, 2*i), Math.pow(-1, i)*Combinatorics.factorial(2*i)));
        }
        return new Matrix(res.getArray());
    }

    // private methods
    // algebraic matrix complement
    private static double cofactor(Matrix A, int row, int col) throws InvalidAttributeValueException {
        int N = A.array.length;
        if (N == 1) return A.array[0][0];

        double[][] minorArray = getMinorArray(A, row, col);
        return Math.pow(-1, row + col) * Matrix.det(new Matrix(minorArray));
    }

    // getting array of minor
    private static double[][] getMinorArray(Matrix A, int skipRow, int skipCol) {
        int N = A.array.length;
        double[][] minor = new double[N-1][N-1];
        int minorRow = 0;

        for (int i = 0; i < N; i++) {
            if (i == skipRow) continue;

            int minorCol = 0;
            for (int j = 0; j < N; j++) {
                if (j == skipCol) continue;

                minor[minorRow][minorCol] = A.array[i][j];
                minorCol++;
            }
            minorRow++;
        }
        return minor;
    }
}
