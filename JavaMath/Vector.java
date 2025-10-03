package JavaMath;

import javax.management.InvalidAttributeValueException;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Function;
import java.util.stream.IntStream;

public class Vector implements Sequence
{
    double[] array;

    public Vector(int len) throws InvalidAttributeValueException
    {
        if(len<=0)
        {
            throw new InvalidAttributeValueException("Length of vector must be natural");
        }

        array = new double[len];
    }
    public Vector(double[] arr)
    {
        this.array = arr.clone();
    }
    public Vector(Vector other)
    {
        this.array = other.array.clone();
    }

    public double[] getArray()
    {
        return this.array;
    }

    public void setArray(double[] array)
    {
        this.array = array;
    }

    public void apply(Function<Double, Double> f)
    {
        this.array = Arrays.stream(this.array).parallel().map(x -> f.apply(x)).toArray();
    }

    public Matrix toMatrix(VectorAsMatrixType type) throws InvalidAttributeValueException {
        double[][] newArray;
        if(type.equals(VectorAsMatrixType.ROW))
        {
            newArray = new double[1][this.array.length];
            IntStream.range(0, this.array.length).forEach(i -> {newArray[0][i]=this.array[i];});
        }
        else if(type.equals(VectorAsMatrixType.COLUMN))
        {
            newArray = new double[this.array.length][1];
            IntStream.range(0, this.array.length).forEach(i -> {newArray[i][0]=this.array[i];});
        }
        else
        {
            newArray = new double[1][1];
        }

        return new Matrix(newArray);
    }

    @Override
    public double norm()
    {
        Vector copy = new Vector(this.array.clone());
        copy.apply(x -> {return x*x;});
        double S = Arrays.stream(copy.array).parallel().sum();

        return Math.sqrt(S);
    }

    @Override
    public double minElement() {
        return Arrays.stream(this.array).parallel().min().getAsDouble();
    }

    @Override
    public double maxElement() {
        return Arrays.stream(this.array).parallel().max().getAsDouble();
    }

    public static double scalarProduct(Vector v1, Vector v2) throws InvalidAttributeValueException {
        if(v1.array.length!=v2.array.length)
        {
            throw new InvalidAttributeValueException("Vectors with different number of elements");
        }
        AtomicReference<Double> res = new AtomicReference<>((double) 0);

        IntStream.range(0, v1.array.length).parallel().map(i -> {
            res.updateAndGet(v -> new Double((double) (v + v1.array[i] * v2.array[i])));
            return i;
        });

        return res.get();
    }
}
