package JavaMath;

import javax.management.InvalidAttributeValueException;

import static JavaMath.Combinatorics.*;

public class Complex {
    double real = 0;
    double image = 0;

    public Complex() {};
    public Complex(double Real)
    {
        this.real = Real;
    }
    public Complex(double Real, double Imag)
    {
        this.real = Real;
        this.image = Imag;
    }

    public boolean equals(Complex other)
    {
        return this.real == other.real && this.image == other.image;
    }
    public String toString()
    {
        return String.format("(%f, %f)", this.real, this.image);
    }

    public double getReal()
    {
        return this.real;
    }
    public void setReal(double real)
    {
        this.real = real;
    }

    public double getImage()
    {
        return this.image;
    }
    public void setImage(double image)
    {
        this.image = image;
    }

    public void add(double value)
    {
        this.real += value;
    }
    public void add(Complex value)
    {
        this.real += value.real;
        this.image += value.image;
    }
    public void subtract(double value)
    {
        this.real -= value;
    }
    public void subtract(Complex value)
    {
        this.real -= value.real;
        this.image -= value.image;
    }
    public void multiply(double value)
    {
        this.real *= value;
    }
    public void multiply(Complex value)
    {
        double a = this.real * value.real - this.image * value.image;
        double b = this.real * value.image + this.image * value.real;

        this.real = a;
        this.image = b;
    }
    public void divide(double value) throws InvalidAttributeValueException
    {
        if(value==0)
        {
            throw new InvalidAttributeValueException("Complex number can't be divided by zero");
        }

        this.real /= value;
        this.image /= value;
    }
    public void divide(Complex value) throws InvalidAttributeValueException
    {
        if(value.real==0 && value.image==0)
        {
            throw new InvalidAttributeValueException("Complex number can't be divided by zero");
        }

        double a = (this.real * value.real + this.image * value.image)/(Math.pow(value.real, 2)+Math.pow(value.image, 2));
        double b = (this.real * value.image - this.image * value.real)/(Math.pow(value.real, 2)+Math.pow(value.image, 2));

        this.real = a;
        this.image = b;
    }
    public double module()
    {
        return Math.sqrt(Math.pow(this.real, 2) + Math.pow(this.image, 2));
    }
    public double argument()
    {
        if(this.real>0 && this.image>0)
        {
            return Math.atan(this.image/this.real);
        }
        else if(this.real<0)
        {
            return Math.PI + Math.atan(this.image/this.real);
        }
        else if(this.real>0 && this.image<0)
        {
            return 2*Math.PI + Math.atan(this.image/this.real);
        }
        else if(this.real==0)
        {
            return Math.PI/2;
        }

        return 0;
    }

    public static Complex add(Complex z1, Complex z2)
    {
        return new Complex(z1.real+ z2.real, z1.image+ z2.image);
    }

    public static Complex add(double r, Complex z)
    {
        return new Complex(z.real+r, z.image);
    }

    public static Complex subtract(Complex z1, Complex z2)
    {
        return new Complex(z1.real-z2.real, z1.image-z2.image);
    }

    public static Complex subtract(double r, Complex z)
    {
        return new Complex(z.real-r, z.image);
    }

    public static Complex multiply(Complex z1, Complex z2)
    {
        double a = z1.real * z2.real - z1.image * z2.image;
        double b = z1.real * z2.image + z1.image * z2.real;

        return new Complex(a, b);
    }

    public static Complex multiply(double r, Complex z)
    {
        return new Complex(r*z.real, r*z.image);
    }

    public static Complex divide(Complex z1, Complex z2) throws InvalidAttributeValueException
    {
        if(z2.real==0 && z2.image==0)
        {
            throw new InvalidAttributeValueException("Complex number can't be divided by zero");
        }

        double a = (z1.real * z2.real + z1.image * z2.image)/(Math.pow(z2.real, 2)+Math.pow(z2.image, 2));
        double b = (z1.real * z2.image - z1.image * z2.real)/(Math.pow(z2.real, 2)+Math.pow(z2.image, 2));

        return new Complex(a, b);
    }

    public static Complex divide(double r, Complex z) throws InvalidAttributeValueException
    {
        if(r==0)
        {
            throw new InvalidAttributeValueException("Complex number can't be divided by zero");
        }

        return new Complex(z.real/r, z.image/r);
    }

    public static Complex pow(Complex z, int n) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        Complex init = new Complex(1, 0);
        int i = 0;

        if(n>0)
        {
            for(i=1; i<=n; i++)
            {
                init.multiply(z);
            }

            return init;
        }

        if(n<0)
        {
            n = -n;
            for(i=0; i<=n; i++)
            {
                init.divide(z);
            }

            return init;
        }

        return init;
    }

    public static Complex exp(Complex z) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        Complex res = new Complex(1, 0);
        Complex power = new Complex();

        for(int i=1; i<=10; i++)
        {
            power = Complex.pow(z, i);
            power.divide((double) factorial(i));
            res.add(power);
        }

        return res;
    }

    public static Complex sin(Complex z) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        Complex i = new Complex(0, 1), inv_i = new Complex(0, -1);
        Complex doubled_i = new Complex(0, 2);
        Complex res = new Complex(1, 0);

        i.multiply(z);
        inv_i.multiply(z);
        res = Complex.exp(i);
        res.subtract(Complex.exp(inv_i));
        res.divide(doubled_i);

        return res;
    }

    public static Complex cos(Complex z) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        Complex i = new Complex(0, 1), inv_i = new Complex(0, -1);
        Complex res = new Complex(1, 0);

        i.multiply(z);
        inv_i.multiply(z);
        res = Complex.exp(i);
        res.add(Complex.exp(inv_i));
        res.divide(2);

        return res;
    }

    public static Complex tan(Complex z) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        Complex res = Complex.sin(z);
        res.divide(Complex.cos(z));
        return res;
    }

    public static Complex cot(Complex z) throws InvalidAttributeValueException, CloneNotSupportedException
    {
        Complex res = Complex.cos(z);
        res.divide(Complex.sin(z));
        return res;
    }

    public static Complex log(Complex z)
    {
        return new Complex(Math.log(z.module()), z.argument());
    }

    public static Complex[] sqrt(Complex z)
    {
        Complex[] res = new Complex[2];
        res[0] = new Complex(Math.sqrt(0.5*(z.module()+z.real)), Math.sqrt(0.5*(z.module()-z.real))*Math.signum(z.image));
        res[1] = new Complex(-Math.sqrt(0.5*(z.module()+z.real)), -Math.sqrt(0.5*(z.module()-z.real))*Math.signum(z.image));

        return res.clone();
    }
}
