public class LagrangeBasisFunction1D extends Function1D
{
	public Cell1D cell;
	public int polynomialDegree;
	public int localFunctionNumber;
	public double degreeOfFreedom;
	public LagrangeBasisFunction1D(int polynomialDegree, int localFunctionNumber, Cell1D cell)
	{
		this.cell = cell;
		this.localFunctionNumber = localFunctionNumber;
		this.polynomialDegree = polynomialDegree;
		this.degreeOfFreedom = 1.0/polynomialDegree*localFunctionNumber; //equidistant Lagrange dofs
	}
	public double valueOnReferenceCell(double pos)
	{
		switch (this.polynomialDegree)
		{
			case 1:
				switch (localFunctionNumber)
				{
					case 0:
						return 1.0-pos;
					case 1:
						return pos;
				}
				break;
			case 2:
				switch (localFunctionNumber)
				{
					case 0:
						return 1.0 - 3.0*pos + 2.0*pos*pos;
					case 1:
						return (pos - pos*pos)*4.0;
					case 2:
						return 2.0*pos*pos-pos;
				}
				break;
			default:
				throw new UnsupportedOperationException();
		}
		return 0;
	}
	public double derivativeOnReferenceCell(double pos)
	{

		switch (this.polynomialDegree)
		{
			case 1:
				switch (localFunctionNumber)
				{
					case 0:
						return -1.0;
					case 1:
						return 1.0;
				}
				break;
			case 2:
				switch (localFunctionNumber)
				{
					case 0:
						return -3.0 + 4.0*pos;
					case 1:
						return (1.0 - 2.0*pos)*4.0;
					case 2:
						return 4.0*pos-1.0;
				}
				break;
			default:
				throw new UnsupportedOperationException();
		}
		return 0;
	}

	@Override
	public double value(double pos)
	{
		if(cell.isInCell(pos))
			return valueOnReferenceCell(cell.positionOnReferenceCell(pos));
		else
			return 0;
	}

	@Override
	public double derivative(double pos)
	{
		if(cell.isInCell(pos))
			return derivativeOnReferenceCell(cell.positionOnReferenceCell(pos))*cell.jacobiDeterminant(pos);
		return 0;
	}
	public void print()
	{
		System.out.println("Shapefunction: Degree " +polynomialDegree+", Local Number "+localFunctionNumber+
			", on Cell:\n\t\t");
		cell.print();
	}
}
