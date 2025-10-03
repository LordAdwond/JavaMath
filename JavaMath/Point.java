package JavaMath;

import javax.management.InvalidAttributeValueException;

public class Point {
    double[] coordinates;

    public Point()
    {
        this.coordinates = new double[0];
    }

    public Point(double[] arr)
    {
        this.coordinates = arr.clone();
    }

    public double[] getCoordinates()
    {
        return this.coordinates;
    }

    public void setCoordinates(double[] coordinates)
    {
        this.coordinates = coordinates;
    }

    public int dim()
    {
        return this.coordinates.length;
    }

    public double component(int k) throws InvalidAttributeValueException
    {
        if(k<1 || k>this.dim())
        {
            throw new InvalidAttributeValueException("Invalid number of component");
        }

        return this.coordinates[k-1];
    }

    @Override
    public String toString()
    {
        String point = "(";
        for(int i=0; i<this.coordinates.length-1; i++)
        {
            point = point + Double.toString(this.coordinates[i]) + "; ";
        }
        point = point + Double.toString(this.coordinates[this.coordinates.length-1]) + ")";

        return point;
    }

    public Vector toVector()
    {
        return new Vector(this.coordinates);
    }
}
