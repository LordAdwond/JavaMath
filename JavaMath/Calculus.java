package JavaMath;

import javax.management.InvalidAttributeValueException;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Function;
import java.util.stream.IntStream;

public class Calculus
{
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
}
