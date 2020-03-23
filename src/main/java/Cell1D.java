import java.util.ArrayList;

public class Cell1D
{
        public double xStart;
        public double xEnd;
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

        public Cell1D(double xStart, double xEnd)
        {
                this.xStart = xStart;
                this.xEnd = xEnd;
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
                return (xEnd - xStart);
        }

        public double center()
        {
                return (xEnd + xStart) / 2;
        }

        public boolean isInCell(double pos)                //in [0,1]
        {
                return (pos >= xStart && pos <= xEnd);
        }

        public double positionOnReferenceCell(double pos)
        {
                return (pos - xStart) / length();
        }

        public double positionOnGrid(double pospp)
        {
                return pospp * length() + xStart;
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
                System.out.println("Cell: start" + xStart + ", end: "+xEnd);
        }
}
