public abstract class AngularFaceIntegral
{
	public static  String JUMPJUMP = "JumpJump";
	AngularScalarFunction weight;
	String name;
	public static String JUMP_NORMALAVERAGE_JUMP = "JumpNormalaverageJump";
	public static String JUMP_NORMALAVERAGE_JUMP_NORMALAVERAGE = "JumpNormalaverageJumpNormalaverage";
	public static String VALUE_VALUE = "ValueValue";
	public static String UPWIND = "Upwind";
	protected AngularFaceIntegral()
	{
	}

	public AngularFaceIntegral(AngularScalarFunction weight, String name)
	{
		this.weight = weight;
		this.name = name;
	}


	public AngularFaceIntegral(String name)
	{
		this(AngularScalarFunction.oneFunction(), name);
	}
	public abstract double evaluateAngularFaceIntegral(Face face,
	                                                   ScalarShapeFunction function1, ScalarShapeFunction function2,
	                                                   DoubleTensor direction1, DoubleTensor direction2,
	                                                   double directionWeight1, double directionWeight2);
}
