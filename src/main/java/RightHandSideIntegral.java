public abstract class RightHandSideIntegral
{
	ScalarFunction rightHandSide;
	ScalarFunction weight;
	public RightHandSideIntegral(ScalarFunction rightHandSide, ScalarFunction weight)
	{
		this.rightHandSide = rightHandSide;
		this.weight = weight;
	}

	public RightHandSideIntegral(ScalarFunction rightHandSide)
	{
		this(rightHandSide, ScalarFunction.oneFunction());
	}

	public abstract double evaluateRightHandSideIntegral(Cell k, ScalarShapeFunction v);
}
