public class AngularTPCell extends TPCell
{
	public AngularTPCell(double xStart, double yStart, double xEnd, double yEnd, int polynomialDegree)
	{
		super(xStart, yStart, xEnd, yEnd, polynomialDegree);
	}
	public AngularTPCell(TPCell cell)
	{
		super(cell.cellx, cell.celly,cell.polynomialDegree);
	}
	public AngularTPCell(Cell1D cellx, Cell1D celly, int polynomialDegree)
	{
		super(cellx, celly, polynomialDegree);
	}
}
