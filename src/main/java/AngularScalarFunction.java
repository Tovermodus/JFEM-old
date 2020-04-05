import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public abstract class AngularScalarFunction
{
	public abstract double value(DoubleTensor posInSpace, DoubleTensor posOnSphere1, DoubleTensor posOnSphere2);
	public double value(DoubleTensor posInSpace, DoubleTensor posOnSphere)
	{
		return this.value(posInSpace,posOnSphere, posOnSphere);
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
	public void plot(int pointres, DoubleTensor[] directions, String filenameWithoutEnding,double xStart,
	                 double yStart, double xEnd, double yEnd)
	{
		double[][][] values = new double[pointres][pointres][directions.length];
		for (int k = 0; k < pointres; k++)
		{
			for (int j = 0; j < pointres; j++)
			{
				for (int l = 0; l < directions.length; l++)
				{
						DoubleTensor pos =
							DoubleTensor.vectorFromValues(xStart+(xEnd - xStart) / (pointres - 1) * k,
								yStart+(yEnd-yStart) / (pointres - 1) * j);
						values[j][k][l] += value(pos,directions[l]);
				}

			}
		}
		try
		{
			for (int k = 0; k < directions.length; k++)
			{
				BufferedWriter plotWriter =
					new BufferedWriter(new FileWriter(filenameWithoutEnding + k + ".dat"));
				for (int i = 0; i < pointres; i++)
				{
					for (int j = 0; j < pointres; j++)
					{
						plotWriter.write(Double.toString(values[i][j][k]) + " ");
					}
					plotWriter.newLine();
				}
				plotWriter.flush();
				plotWriter.close();
			}
		} catch (
			IOException e)
		{
			e.printStackTrace();
		}
	}
}
