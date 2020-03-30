import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class FESpaceFunction extends ScalarFunction
{
	Map<ScalarShapeFunction, Double> values;
	public FESpaceFunction(ScalarShapeFunction functions [], double coefficients [])
	{
		assert(functions.length == coefficients.length);
		values = new HashMap<>();
		for(int i = 0; i < functions.length; i++)
		{
			values.put(functions[i], coefficients[i]);
		}
	}
	public FESpaceFunction(ArrayList<ScalarShapeFunction> functions, DoubleTensor coefficients)
	{
		assert(functions.size() == coefficients.size());
		values = new HashMap<>();
		for(int i = 0; i < functions.size(); i++)
		{
			values.put(functions.get(i), coefficients.at(i));
		}
	}
	@Override
	public double value(DoubleTensor pos)
	{
		return values.entrySet().stream().parallel().mapToDouble(entry->entry.getKey().value(pos)*entry.getValue()).sum();
	}
	@Override
	public DoubleTensor derivative(DoubleTensor pos)
	{
		return null;
	}
}
