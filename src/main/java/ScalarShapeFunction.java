import java.util.ArrayList;
import java.util.Map;

public abstract class ScalarShapeFunction extends ScalarFunction
{
	ScalarNodeFunctional nodeFunctional;
	ArrayList<Cell> cells;
	int globalIndex;
	void setGlobalIndex(int index)
	{
		this.globalIndex = index;
	}
	abstract double valueInCell(DoubleTensor pos, Cell cell);
	abstract DoubleTensor derivativeInCell(DoubleTensor pos, Cell cell);
	abstract Map<Integer, Double> prolongate(ArrayList<ScalarShapeFunction> refinedFunctions);
}
