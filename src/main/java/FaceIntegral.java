

public abstract class FaceIntegral
{
	ScalarFunction weight;
	TensorFunction vectorWeight;
	String name;
	public FaceIntegral(ScalarFunction weight, String name)
	{
		this.weight = weight;
		this.name = name;
	}
	public FaceIntegral(TensorFunction vectorWeight, String name)
	{
		this.vectorWeight = vectorWeight;
		this.name = name;
	}
	public FaceIntegral(String name)
	{
		this(ScalarFunction.oneFunction(), name);
	}

	abstract double evaluateFaceIntegral(Face face, ScalarShapeFunction function1,
	                                     ScalarShapeFunction function2);

}
