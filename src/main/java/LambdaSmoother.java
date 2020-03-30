import com.google.common.collect.Iterables;

public class LambdaSmoother extends AngularSmoother
{
	double isoScatter;
	public LambdaSmoother(AngularGrid g, String[] args)
	{
		super(g, args);
		isoScatter = Double.parseDouble(args[0]);
	}

	@Override
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		System.out.println(isoScatter);
		DoubleTensor scatteredIterate = new DoubleTensor(iterate.size());
		/*for(int i = 0; i < grid.directions.length; i++)
		{
			for(int j = 0; j < grid.g.shapeFunctions.size(); j++)
				for(int k = 0; k < grid.directions.length; k++)
					scatteredIterate.add(j*grid.directions.length+i,
						grid.direction_weights[k]*isoScatter*iterate.at(j*grid.directions.length+k)/3);
		}*/
		scatteredIterate = grid.Scat.mvmul(iterate).mul(-1);
		AngularFESpaceFunction solfunc =
			new AngularFESpaceFunction(grid.g.shapeFunctions,
				grid.directions,
				iterate);
		solfunc.plot(100,grid.directions,"/home/tovermodus/plot");
		solfunc =
			new AngularFESpaceFunction(grid.g.shapeFunctions,
				grid.directions,
				scatteredIterate);
		solfunc.plot(100,grid.directions,"/home/tovermodus/scatter");
		DoubleTensor newRight = scatteredIterate.mul(1).add(rightHandSide);
		solfunc =
			new AngularFESpaceFunction(grid.g.shapeFunctions,
				grid.directions,
				newRight);
		solfunc.plot(100,grid.directions,"/home/tovermodus/right");
		iterate = grid.TransAbs.solve(scatteredIterate.add(rightHandSide));
		return iterate;
	}
}
