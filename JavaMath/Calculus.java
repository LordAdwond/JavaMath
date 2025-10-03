package JavaMath;

import javax.management.InvalidAttributeValueException;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Function;
import java.util.stream.IntStream;

public class Calculus
{
    // calculus tasks
    public static double derivate(Function<Double, Double> f, double x)
    {
        double h = Double.MIN_VALUE;
        return (f.apply(x+h)-f.apply(x))/h;
    }

    public static double integrate(Function<Double, Double> f, double a, double b)
    {
        if(a==b)
        {
            return 0;
        }
        if(a>b)
        {
            return Calculus.integrate(f, b, a);
        }

        int n = 10000;
        double h = (b-a)/n;
        AtomicReference<Double> res = new AtomicReference<>((double) 0.5*(f.apply(a)+f.apply(b)));

        IntStream.range(1, n-1).parallel().map(i -> {
            res.updateAndGet(v -> new Double(v+f.apply(a+i*h)));
            return i;
        });

        return res.get();
    }

    public static double interpolate(double[] Xk, double[] Yk, double x) throws InvalidAttributeValueException
    {
        if(Xk.length!=Yk.length)
        {
            throw new InvalidAttributeValueException("A number of Xs and a number of Ys are not same");
        }

        double result = 0, S=0;
        int i, j;

        for(i=0; i<Xk.length; i++)
        {
            S=1;
            for(j=0; j<Xk.length; j++)
            {
                if(i!=j)
                {
                    S *= (x-Xk[j])/(Xk[i]-Xk[j]);
                }
                result += Yk[i]*S;
            }
        }

        return result;
    }

    public static double interpolate(Point[] Pk, double x) throws InvalidAttributeValueException
    {
        double result = 0, S=0;
        int i, j;

        for(i=0; i<Pk.length; i++)
        {
            if(Pk[i].dim()!=2)
            {
                throw new InvalidAttributeValueException("All points must have two numbers");
            }
        }

        for(i=0; i<Pk.length; i++)
        {
            S=1;
            for(j=0; j<Pk.length; j++)
            {
                if(i!=j)
                {
                    S *= (x-Pk[j].component(0))/(Pk[i].component(0)-Pk[j].component(0));
                }
            }
            result += Pk[i].component(1)*S;
        }

        return result;
    }

    //discrete transforms
    // Discrete Fourier Transform
    public static Complex[] DFT1d(double[] f) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        if(f.length==0)
        {
            throw new InvalidAttributeValueException("Array of values of the function can't be empty");
        }

        int N = f.length;
        Complex I = new Complex(0, -1);
        Complex[] res = new Complex[f.length];
        int i, j;

        for(i=0; i<f.length; i++)
        {
            for(j=0; j<f.length; j++)
            {
                res[i] = Complex.add(res[i], Complex.multiply(f[j], Complex.exp(Complex.multiply(2*Math.PI*i*j/N, I))));
            }
            res[i] = Complex.divide(N, res[i]);
        }

        return res;
    }

    public static Complex[][] DFT2d(double[][] f) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        Matrix A = new Matrix(f);
        int N = f.length, M=f[0].length;
        Complex I = new Complex(0, -1);
        Complex[][] res = new Complex[f.length][f[0].length];
        int i, j, k, l;

        for(i=0; i<N; i++)
        {
            for(j=0; j<M; j++)
            {

                for(k=0; k<N; k++)
                {
                    for(l=0; l<M; l++)
                    {
                        res[i][j] = Complex.add(res[i][j], Complex.multiply(f[k][l], Complex.exp(Complex.multiply(2*Math.PI*(i*k+j*l)/N, I))));
                    }
                }
                res[i][j] = Complex.divide(N*M, res[i][j]);

            }
        }

        return res;
    }

    public static Complex[][] DFT2d(Complex[][] f) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        int N = f.length, M=f[0].length;
        Complex I = new Complex(0, -1);
        Complex[][] res = new Complex[f.length][f[0].length];
        int i, j, k, l;

        for(i=0; i<N; i++)
        {
            for(j=0; j<M; j++)
            {

                for(k=0; k<N; k++)
                {
                    for(l=0; l<M; l++)
                    {
                        res[i][j] = Complex.add(res[i][j], Complex.multiply(f[k][l], Complex.exp(Complex.multiply(2*Math.PI*(i*k+j*l)/N, I))));
                    }
                }
                res[i][j] = Complex.divide(N*M, res[i][j]);

            }
        }

        return res;
    }

    public static Complex[][] DFT2d(Matrix f) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        return DFT2d(f.getArray());
    }

    public static Complex[] DFT1d(Complex[] f) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        if(f.length==0)
        {
            throw new InvalidAttributeValueException("Array of values of the function can't be empty");
        }

        int N = f.length;
        Complex I = new Complex(0, -1);
        Complex[] res = new Complex[f.length];
        int i, j;

        for(i=0; i<f.length; i++)
        {
            for(j=0; j<f.length; j++)
            {
                res[i] = Complex.add(res[i], Complex.multiply(f[j], Complex.exp(Complex.multiply(2*Math.PI*i*j/N, I))));
            }
            res[i] = Complex.divide(N, res[i]);
        }

        return res;
    }

    public static Complex[] DFT1d(Vector v) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        return DFT1d(v.getArray());
    }

    public static Complex[] IDFT1d(double[] f) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        if(f.length==0)
        {
            throw new InvalidAttributeValueException("Array of values of the function can't be empty");
        }

        int N = f.length;
        Complex I = new Complex(0, 1);
        Complex[] res = new Complex[f.length];
        int i, j;

        for(i=0; i<f.length; i++)
        {
            for(j=0; j<f.length; j++)
            {
                res[i] = Complex.add(res[i], Complex.multiply(f[j], Complex.exp(Complex.multiply(2*Math.PI*i*j/N, I))));
            }
            res[i] = Complex.divide(N, res[i]);
        }

        return res;
    }

    public static Complex[] IDFT1d(Complex[] f) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        if(f.length==0)
        {
            throw new InvalidAttributeValueException("Array of values of the function can't be empty");
        }

        int N = f.length;
        Complex I = new Complex(0, 1);
        Complex[] res = new Complex[f.length];
        int i, j;

        for(i=0; i<f.length; i++)
        {
            for(j=0; j<f.length; j++)
            {
                res[i] = Complex.add(res[i], Complex.multiply(f[j], Complex.exp(Complex.multiply(2*Math.PI*i*j/N, I))));
            }
            res[i] = Complex.divide(N, res[i]);
        }

        return res;
    }

    public static Complex[] IDFT1d(Vector v) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        return IDFT1d(v.getArray());
    }

    public static Complex[][] IDFT2d(double[][] f) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        Matrix A = new Matrix(f);
        int N = f.length, M=f[0].length;
        Complex I = new Complex(0, 1);
        Complex[][] res = new Complex[f.length][f[0].length];
        int i, j, k, l;

        for(i=0; i<N; i++)
        {
            for(j=0; j<M; j++)
            {

                for(k=0; k<N; k++)
                {
                    for(l=0; l<M; l++)
                    {
                        res[i][j] = Complex.add(res[i][j], Complex.multiply(f[k][l], Complex.exp(Complex.multiply(2*Math.PI*(i*k+j*l)/N, I))));
                    }
                }
                res[i][j] = Complex.divide(N*M, res[i][j]);

            }
        }

        return res;
    }

    public static Complex[][] IDFT2d(Complex[][] f) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        int N = f.length, M=f[0].length;
        Complex I = new Complex(0, 1);
        Complex[][] res = new Complex[f.length][f[0].length];
        int i, j, k, l;

        for(i=0; i<N; i++)
        {
            for(j=0; j<M; j++)
            {

                for(k=0; k<N; k++)
                {
                    for(l=0; l<M; l++)
                    {
                        res[i][j] = Complex.add(res[i][j], Complex.multiply(f[k][l], Complex.exp(Complex.multiply(2*Math.PI*(i*k+j*l)/N, I))));
                    }
                }
                res[i][j] = Complex.divide(N*M, res[i][j]);

            }
        }

        return res;
    }

    public static Complex[][] IDFT2d(Matrix f) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        return IDFT2d(f.getArray());
    }

    // Discrete Laplace Transform
    public static double[] DLT1d(double[] f) throws InvalidAttributeValueException
    {
        if(f.length==0)
        {
            throw new InvalidAttributeValueException("Array of values of the function can't be empty");
        }

        int N = f.length;
        double[] res = new double[N];
        int i, j;

        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                res[i] += f[j]*Math.exp(-2*Math.PI*i*j/N);
            }
            res[i] = res[i]/N;
        }

        return res;
    }

    public static double[] DLT1d(Vector v) throws InvalidAttributeValueException
    {
        return DLT1d(v.getArray());
    }

    public static double[] IDLT1d(double[] f) throws InvalidAttributeValueException
    {
        if(f.length==0)
        {
            throw new InvalidAttributeValueException("Array of values of the function can't be empty");
        }

        int N = f.length;
        double[] res = new double[N];
        int i, j;

        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                res[i] += f[j]*Math.exp(2*Math.PI*i*j/N);
            }
            res[i] = res[i]/N;
        }

        return res;
    }

    public static double[] IDLT1d(Vector v) throws InvalidAttributeValueException
    {
        return IDLT1d(v.getArray());
    }

    public static double[][] DLT2d(double[][] f) throws InvalidAttributeValueException
    {
        Matrix A = new Matrix(f);
        int N = f.length, M=f[0].length;
        double[][] res = new double[f.length][f[0].length];
        int i, j, k, l;

        for(i=0; i<N; i++)
        {
            for(j=0; j<M; j++)
            {

                for(k=0; k<N; k++)
                {
                    for(l=0; l<M; l++)
                    {
                        res[i][j] = Math.exp(-2*Math.PI*(i*k+j*l)/(N*M));
                    }
                }
                res[i][j] /= N*M;

            }
        }

        return res;
    }

    public static double[][] DLT2d(Matrix f) throws InvalidAttributeValueException
    {
        return DLT2d(f.getArray());
    }

    public static double[][] IDLT2d(double[][] f) throws InvalidAttributeValueException
    {
        Matrix A = new Matrix(f);
        int N = f.length, M=f[0].length;
        double[][] res = new double[f.length][f[0].length];
        int i, j, k, l;

        for(i=0; i<N; i++)
        {
            for(j=0; j<M; j++)
            {

                for(k=0; k<N; k++)
                {
                    for(l=0; l<M; l++)
                    {
                        res[i][j] = Math.exp(2*Math.PI*(i*k+j*l)/(N*M));
                    }
                }
                res[i][j] /= N*M;

            }
        }

        return res;
    }

    public static double[][] IDLT2d(Matrix f) throws InvalidAttributeValueException
    {
        return IDLT2d(f.getArray());
    }
}
