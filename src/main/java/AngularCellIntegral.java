public abstract class AngularCellIntegral
{
	AngularScalarFunction weight;
	String name;
	public static String TRANSPORT="Transport";
	public static String SCATTERING="Scattering";
	public static String ABSORPTION="Absorption";
	protected AngularCellIntegral()
	{
	}

	public AngularCellIntegral(AngularScalarFunction weight, String name)
	{
		this.weight = weight;
		this.name = name;
	}

	public AngularCellIntegral(String name)
	{
		this(AngularScalarFunction.oneFunction(), name);
	}
	public abstract double evaluateAngularCellIntegral(Cell cell,
	                                          ScalarShapeFunction function1, ScalarShapeFunction function2,
	                                          DoubleTensor direction1, DoubleTensor direction2,
	                                          double directionWeight1, double directionWeight2);
}
