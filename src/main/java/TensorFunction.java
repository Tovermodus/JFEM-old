public abstract class TensorFunction
{
	public int dimension;
	public abstract DoubleTensor value(DoubleTensor pos);

	public abstract DoubleTensor derivative(DoubleTensor pos);

}
