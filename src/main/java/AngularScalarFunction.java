public abstract class AngularScalarFunction
{
	public abstract double value(DoubleTensor posInSpace, DoubleTensor posOnSphere1, DoubleTensor posOnSphere2);
	public double value(DoubleTensor posInSpace, DoubleTensor posOnSphere)
	{
		return this.value(posInSpace,posOnSphere, new DoubleTensor(posOnSphere.size()));
	}

	public static AngularScalarFunction oneFunction()
	{
		return constantFunction(1);
	}
	public static AngularScalarFunction constantFunction(double constant)
	{
		AngularScalarFunction oneF = new AngularScalarFunction()
		{
			@Override
			public double value(DoubleTensor posInSpace, DoubleTensor posOnSphere1,
			                    DoubleTensor posOnSphere2)
			{
				return constant;
			}

		};
		return oneF;
	}
}
