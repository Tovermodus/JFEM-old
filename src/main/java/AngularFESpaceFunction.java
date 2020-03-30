import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class AngularFESpaceFunction extends AngularScalarFunction
{
	Map<Pair<ScalarShapeFunction,DoubleTensor>, Double> values;
	public AngularFESpaceFunction(ScalarShapeFunction functions [], DoubleTensor directions[],
	                              double coefficients [])
	{
		assert(functions.length*directions.length == coefficients.length);
		values = new HashMap<>();
		for(int i = 0; i < functions.length; i++)
		{
			for(int j = 0; j < directions.length; j++)
			values.put(new Pair<>(functions[i],directions[j]), coefficients[i*directions.length+j]);
		}
	}
	public AngularFESpaceFunction(ArrayList<ScalarShapeFunction> functions, DoubleTensor directions[],
	                              DoubleTensor coefficients)
	{
		assert(functions.size()*directions.length == coefficients.size());
		values = new HashMap<>();
		for(int i = 0; i < functions.size(); i++)
		{
			for(int j = 0; j < directions.length; j++)
				values.put(new Pair<>(functions.get(i),directions[j]),
					coefficients.at(i*directions.length+j));
		}
	}
	@Override
	public double value(DoubleTensor posInSpace, DoubleTensor posOnSphere1, DoubleTensor posOnSphere2)
	{
		return values.entrySet().stream().parallel().filter(entry->entry.getKey().getSecond().equals(posOnSphere1)).mapToDouble(entry->entry.getKey().getFirst().value(posInSpace)*entry.getValue()).sum();
	}
}
