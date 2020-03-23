import java.util.ArrayList;

public abstract class Face
{
	ArrayList<Cell> cells;
	ArrayList<ScalarShapeFunction> shapeFunctions;
	TensorFunction normal;
	boolean isBoundaryFace;
	public Face(TensorFunction normal)
	{
		cells = new ArrayList<>();
		shapeFunctions = new ArrayList<>();
		this.normal = normal;
	}
	abstract Cell getNormalDownstreamCell(DoubleTensor pos);
	abstract Cell getNormalUpstreamCell(DoubleTensor pos);
	abstract DoubleTensor center();
	abstract boolean isOnFace(DoubleTensor pos);
	double jumpInValue(ScalarShapeFunction func, DoubleTensor pos)
	{
		return func.valueInCell(pos,getNormalUpstreamCell(pos)) - func.valueInCell(pos,
			getNormalDownstreamCell(pos));
	}

	DoubleTensor jumpInDerivative(ScalarShapeFunction func, DoubleTensor pos)
	{
		return func.derivativeInCell(pos,getNormalUpstreamCell(pos)).sub(
			func.derivativeInCell(pos, getNormalDownstreamCell(pos)));
	}

	double averageInValue(ScalarShapeFunction func, DoubleTensor pos)
	{
		return 0.5*(func.valueInCell(pos,getNormalUpstreamCell(pos))+func.valueInCell(pos,
			getNormalDownstreamCell(pos)));
	}

	DoubleTensor averageInDerivative(ScalarShapeFunction func, DoubleTensor pos)
	{
		return func.derivativeInCell(pos,getNormalUpstreamCell(pos)).add(
			func.derivativeInCell(pos,getNormalDownstreamCell(pos))).mul(0.5);
	}
	DoubleTensor normalAverageInValue(ScalarShapeFunction func, DoubleTensor pos)
	{
		return normal.value(pos).mul(0.5*jumpInValue(func, pos));
	}
	double normalAverageInDerivative(ScalarShapeFunction func, DoubleTensor pos)
	{
		return normal.value(pos).inner(jumpInDerivative(func, pos));
	}

}
