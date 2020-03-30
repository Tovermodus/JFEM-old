import com.sun.source.tree.BreakTree;

public class RichardsonSmoother extends Smoother
{
	double omega;
	public RichardsonSmoother(Grid g, String[] args)
	{
		super(g, args);
		this.omega = Double.parseDouble(args[0]);
	}
	@Override
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		DoubleTensor residuum = rightHandSide.sub(grid.A.mvmul(iterate));
		System.out.println(residuum.vectorNorm());
		iterate = iterate.add(residuum).mul(omega);
		return iterate;

	}
}
