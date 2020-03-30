public class Smoother
{
	Grid grid;
	String[] args;
	public Smoother(Grid g, String[] args)
	{
		this.grid = g;
		this.args = args;
	}
	public static Smoother fromSmoother(Smoother s)
	{
		return new Smoother(s.grid,s.args);
	}
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		throw new UnsupportedOperationException();
	}
}

