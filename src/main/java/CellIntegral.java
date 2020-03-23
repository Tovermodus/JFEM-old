
public abstract class CellIntegral
{
	ScalarFunction weight;
	TensorFunction vectorWeight;
	String name;
	public CellIntegral(ScalarFunction weight, String name)
	{
		this.weight = weight;
		this.name = name;
	}
	public CellIntegral(TensorFunction vectorWeight, String name)
	{
		this.vectorWeight = vectorWeight;
		this.name = name;
	}
	public CellIntegral(String name)
	{
		this(ScalarFunction.oneFunction(), name);
	}
	public abstract double evaluateCellIntegral(Cell cell, ScalarShapeFunction shapeFunction1,
	                                            ScalarShapeFunction shapeFunction2);

}
