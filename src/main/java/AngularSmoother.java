public class AngularSmoother
{
	AngularGrid grid;
	String[] args;
	public AngularSmoother(AngularGrid g, String[] args)
	{
		this.grid = g;
		this.args = args;
	}
	public static AngularSmoother fromSmoother(AngularSmoother s)
	{
		return new AngularSmoother(s.grid,s.args);
	}
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		throw new UnsupportedOperationException();
	}
}

