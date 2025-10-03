package JavaMath;

import javax.management.InvalidAttributeValueException;

public class Combinatorics {
    public static int factorial(int n)
    {
        if(n>0)
        {
            return n*factorial(n-1);
        }

        if(n<0)
        {
            return n*factorial(n+1);
        }

        return 1;
    }

    public static int C(int n, int k) throws InvalidAttributeValueException {
        if(n<1 || k<1 || k>n)
        {
            throw new InvalidAttributeValueException("Incorrect input values");
        }

        return factorial(n) / (factorial(k)*factorial(n-k));
    }

    public static int A(int n, int k) throws InvalidAttributeValueException {
        if(n<1 || k<1 || k>n)
        {
            throw new InvalidAttributeValueException("Incorrect input values");
        }

        return factorial(n) / factorial(n-k);
    }

    public static int CWithRepetitions(int n, int k) throws InvalidAttributeValueException {
        if(n<1 || k<1)
        {
            throw new InvalidAttributeValueException("Incorrect input values");
        }

        return factorial(n+k-1) / (factorial(k)*factorial(n-1));
    }

    public static int AWithRepetitions(int n, int k) throws InvalidAttributeValueException {
        if(n<1 || k<1)
        {
            throw new InvalidAttributeValueException("Incorrect input values");
        }

        return (int) Math.pow(n, k);
    }
}
