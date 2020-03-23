

public class LagrangeNodeFunctional extends ScalarNodeFunctional
{
	DoubleTensor point;
	public LagrangeNodeFunctional(DoubleTensor point)
	{
		this.point = point;
	}
	@Override
	public double evaluate(ScalarFunction func)
	{
		return func.value(point);
	}
}
