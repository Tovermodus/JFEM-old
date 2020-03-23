public abstract class BoundaryFaceIntegral
{
	ScalarFunction rightHandSide;
	ScalarFunction weight;
	public BoundaryFaceIntegral(ScalarFunction rightHandSide,ScalarFunction weight)
	{
		this.rightHandSide = rightHandSide;
		this.weight = weight;
	}

	public abstract double evaluateBoundaryFaceIntegral(Face f, ScalarShapeFunction v);
}
