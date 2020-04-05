public class AllSmoother extends AngularSmoother
{
	LambdaSmoother l;
	DownwindGaussSeidelSmoother d;
	BlockJacobiSmoother b;

	public AllSmoother(AngularGrid g, String[] args)
	{
		super(g, args);
		l = new LambdaSmoother(g,args);
		d = new DownwindGaussSeidelSmoother(g,args);
		b = new BlockJacobiSmoother(g,args);
	}
	@Override
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		double resnorm = rightHandSide.sub(grid.A.mvmul(iterate)).vectorNorm();
		DoubleTensor itcopy = new DoubleTensor(iterate);
		DoubleTensor lambIt = l.smooth(itcopy,rightHandSide);
		double lanorm = rightHandSide.sub(grid.A.mvmul(lambIt)).vectorNorm();
		if(lanorm < resnorm)
		{
			resnorm = lanorm;
			iterate = new DoubleTensor(lambIt);
		}
		itcopy = new DoubleTensor(iterate);
		lambIt = d.smooth(itcopy,rightHandSide);
		lanorm = rightHandSide.sub(grid.A.mvmul(lambIt)).vectorNorm();
		if(lanorm < resnorm)
		{
			resnorm = lanorm;
			iterate = new DoubleTensor(lambIt);
		}
		itcopy = new DoubleTensor(iterate);
		lambIt = b.smooth(itcopy,rightHandSide);
		lanorm = rightHandSide.sub(grid.A.mvmul(lambIt)).vectorNorm();
		if(lanorm < resnorm)
		{
			iterate = new DoubleTensor(lambIt);
			resnorm = lanorm;
		}
		return iterate;
	}
}
