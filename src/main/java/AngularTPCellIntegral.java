public class AngularTPCellIntegral
{
	AngularScalarFunction weight;
	String name;
	public static String TRANSPORT="Transport";
	public static String SCATTERING="Scattering";
	public static String ABSORPTION="Absorption";
	protected AngularTPCellIntegral()
	{
	}

	public AngularTPCellIntegral(AngularScalarFunction weight, String name)
	{
		this.weight = weight;
		this.name = name;
	}

	public AngularTPCellIntegral(String name)
	{
		this(AngularScalarFunction.oneFunction(), name);
	}
	public double evaluateAngularCellIntegral(TPCell cell,
	                                          ScalarShapeFunction function1, ScalarShapeFunction function2,
	                                          DoubleTensor direction1, DoubleTensor direction2,
	                                          double directionWeight1, double directionWeight2)
	{
		if(name.equals(TRANSPORT))
		{
			if(direction1.sub(direction2).vectorNorm()>=1e-14)
				return 0.;
			double ret = 0;

			for(int i = 0; i < cell.cellx.weights.length; i++)
				for(int j = 0; j < cell.celly.weights.length; j++)
				{
					DoubleTensor point = DoubleTensor.vectorFromValues(cell.cellx.points[i],
						cell.celly.points[j]);
					double integralWeight = cell.cellx.weights[i]*cell.celly.weights[j];

					ret += direction1.inner(function1.derivative(point))*function2.value(point)*weight.value(point,direction1)*directionWeight1*integralWeight;
				}
			return ret;
		}
		if(name.equals(SCATTERING))
		{
			double ret = 0;
			for(int i = 0; i < cell.cellx.weights.length; i++)
				for(int j = 0; j < cell.celly.weights.length; j++)
				{
					DoubleTensor point = DoubleTensor.vectorFromValues(cell.cellx.points[i],
						cell.celly.points[j]);
					double integralWeight = cell.cellx.weights[i]*cell.celly.weights[j];
					ret -= directionWeight1*directionWeight2*weight.value(point,direction1,
						direction2)*function1.value(point)*function2.value(point)*integralWeight;
				}
			return ret;
		}
		if(name.equals(ABSORPTION))
		{
			if(direction1.sub(direction2).vectorNorm()>=1e-14)
				return 0.;
			double ret = 0;
			for(int i = 0; i < cell.cellx.weights.length; i++)
				for(int j = 0; j < cell.celly.weights.length; j++)
				{
					DoubleTensor point = DoubleTensor.vectorFromValues(cell.cellx.points[i],
						cell.celly.points[j]);
					double integralWeight = cell.cellx.weights[i]*cell.celly.weights[j];
					ret += directionWeight1*weight.value(point,direction2,
						direction1)*function1.value(point)*function2.value(point)*integralWeight;
				}
			return ret;
		}
		throw new UnsupportedOperationException();
	}
}
