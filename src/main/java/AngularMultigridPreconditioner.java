import java.lang.reflect.Constructor;
import java.util.ArrayList;

public class AngularMultigridPreconditioner implements VectorMultiplication
{
	AngularGridHierarchy angularGridHierarchy;

	ArrayList<AngularSmoother> smoothers;
	public AngularMultigridPreconditioner(AngularGridHierarchy angularGridHierarchy, Class<? extends AngularSmoother>  smoother,String[] args)
	{
		this.angularGridHierarchy = angularGridHierarchy;

		smoothers = new ArrayList<>();
		Constructor<? extends AngularSmoother> constr = null;
		try
		{
			constr = smoother.getConstructor(AngularGrid.class, String[].class);
		} catch (NoSuchMethodException e)
		{
			e.printStackTrace();
		}
		for (AngularGrid grid : angularGridHierarchy.angularGrids)
		{
			try
			{
				assert constr != null;
				smoothers.add(constr.newInstance(grid, args));
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		}
	}
	@Override
	public DoubleTensor mvmul(DoubleTensor vector)
	{
		return multiGridV(angularGridHierarchy.prolongationMatrices.size(),vector,smoothers);
	}

	public DoubleTensor multiGridV(int level, DoubleTensor rightHandSide, ArrayList<AngularSmoother> angularSmoothers)
	{
		int preIters = 1;
		int postIters = 1;
		DoubleTensor iterate = new DoubleTensor((int)rightHandSide.size());
		AngularGrid g = angularGridHierarchy.angularGrids.get(level);
		AngularSmoother angularSmoother = angularSmoothers.get(level);
		//rightHandSide.print_formatted("rhs");
		if(level == 0)
		{
			//g.A.solve(rightHandSide).print_formatted("level0solve");
			return g.A.solve(rightHandSide);
		}
		for(int k = 0; k < preIters; k++)
		{
			iterate = angularSmoother.smooth(iterate, rightHandSide);
		}
		DoubleTensor restrictedRightHandSide =
			angularGridHierarchy.prolongationMatrices.get(level - 1).tvmul(rightHandSide.sub(g.A.mvmul(iterate)));
		DoubleTensor correction = multiGridV(level - 1,
			restrictedRightHandSide,angularSmoothers);
		iterate = iterate.add(angularGridHierarchy.prolongationMatrices.get(level - 1).mvmul(correction));
		for(int k = 0; k < postIters; k++)
		{
			iterate = angularSmoother.smooth(iterate, rightHandSide);
		}
		//iterate.print_formatted("return");
		return iterate;

	}

}
