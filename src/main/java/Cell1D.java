import java.util.ArrayList;

public class Cell1D
{
        public double start;
        public double end;
        public double[] referenceWeights = {
                (322.0 - 13 * Math.sqrt(70)) / 900,
                (322.0 + 13 * Math.sqrt(70)) / 900,
                128.0 / 225,
                (322.0 + 13 * Math.sqrt(70)) / 900,
                (322.0 - 13 * Math.sqrt(70)) / 900};
        public double[] referencePoints = {                             //auf [-1,1]
                -Math.sqrt(5.0 + 2 * Math.sqrt(10.0 / 7)) / 3,
                -Math.sqrt(5.0 - 2 * Math.sqrt(10.0 / 7)) / 3,
                0,
                Math.sqrt(5.0 - 2 * Math.sqrt(10.0 / 7)) / 3,
                Math.sqrt(5.0 + 2 * Math.sqrt(10.0 / 7)) / 3};
        public double[] points;
        public double[] weights;
        public ArrayList<LagrangeBasisFunction1D> shapefunctions;

        public Cell1D(double start, double end)
        {
                this.start = start;
                this.end = end;
                points = new double[referencePoints.length];
                weights = new double[referencePoints.length];
                for(int i = 0; i < points.length; i++)
                {
                       points[i] = referencePoints[i]*length()/2+center();
                       weights[i] = referenceWeights[i]*length()/2;
                }
        }

        public double length()
        {
                return (end - start);
        }

        public double center()
        {
                return (end + start) / 2;
        }

        public boolean isInCell(double pos)                //in [0,1]
        {
                return (pos >= start && pos <= end);
        }

        public double positionOnReferenceCell(double pos)
        {
                return (pos - start) / length();
        }

        public double positionOnGrid(double pospp)
        {
                return pospp * length() + start;
        }

        public double jacobiDeterminant(double pos)
        {
                return 1.0 / length();
        }

        public void distributeFunctions(int polynomialDegree)
        {
                for(int  i = 0; i < polynomialDegree + 1; i++)
                {
                        shapefunctions.add(new LagrangeBasisFunction1D(polynomialDegree, i, this));
                }
        }
        public void print()
        {
                System.out.println("Cell: start" + start + ", end: "+ end);
        }
}
