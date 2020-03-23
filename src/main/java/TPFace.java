import javax.swing.*;

public class TPFace extends Face
{
	Cell1D cell1d;
	double otherCoordinate;
	int normaldirection;
	public TPFace(Cell1D cell1d, double otherCoordinate, int normaldirection)
	{
		super(new TensorFunction()
		{
			@Override
			public DoubleTensor value(DoubleTensor pos)
			{
				if(normaldirection == 0)
					return DoubleTensor.vectorFromValues(1,0);
				else
					return DoubleTensor.vectorFromValues(0,1);
			}

			@Override
			public DoubleTensor derivative(DoubleTensor pos)
			{
				return DoubleTensor.squareMatrixFromValues(0,0,0,0);
			}
		});
		this.cell1d = cell1d;
		this.otherCoordinate = otherCoordinate;
		this.normaldirection = normaldirection;
	}

	public Cell getUpStreamCell(DoubleTensor pos, DoubleTensor direction)
	{
		if(normal.value(pos).inner(direction)>0)
			return getNormalUpstreamCell(pos);
		else
			return getNormalDownstreamCell(pos);
	}
	public Cell getDownStreamCell(DoubleTensor pos, DoubleTensor direction)
	{
		if(normal.value(pos).inner(direction)<0)
			return getNormalUpstreamCell(pos);
		else
			return getNormalDownstreamCell(pos);
	}
	@Override
	Cell getNormalUpstreamCell(DoubleTensor pos)
	{
		for(Cell cell:cells)
		{
			if(cell.isInCell(pos))
			{
				if (cell.center().at(normaldirection) < otherCoordinate)
					return cell;
			}
		}
		return null;
	}

	@Override
	Cell getNormalDownstreamCell(DoubleTensor pos)
	{
		for(Cell cell:cells)
		{

			if(cell.isInCell(pos))
			{
				if (cell.center().at(normaldirection) > otherCoordinate)
					return cell;
			}
		}
		return null;
	}

	@Override
	DoubleTensor center()
	{
		if(normaldirection == 0)
			return DoubleTensor.vectorFromValues(otherCoordinate,cell1d.center());
		else
			return DoubleTensor.vectorFromValues(cell1d.center(),otherCoordinate);

	}

	@Override
	boolean isOnFace(DoubleTensor pos)
	{
		return pos.at(normaldirection) == otherCoordinate && cell1d.isInCell(pos.at(1-normaldirection));
	}

	public void print()
	{
		if(this.normaldirection == 0)
		{
			System.out.println(otherCoordinate+"×["+cell1d.xStart+","+cell1d.xEnd+"]");
		}
		else
		{
			System.out.println("["+cell1d.xStart+","+cell1d.xEnd+"]x"+otherCoordinate);
		}

	}




}
