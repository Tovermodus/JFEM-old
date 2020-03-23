
public abstract class ScalarFunction
{
	public abstract double value(DoubleTensor pos);

	public abstract DoubleTensor derivative(DoubleTensor pos);
	public static ScalarFunction oneFunction()
	{
		return constantFunction(1);
	}
	public static ScalarFunction constantFunction(double constant)
	{
		ScalarFunction oneF = new ScalarFunction()
		{
			@Override
			public double value(DoubleTensor pos)
			{
				return constant;
			}

			@Override
			public DoubleTensor derivative(DoubleTensor pos)
			{
				return new DoubleTensor(pos.size());
			}
		};
		return oneF;
	}
}
