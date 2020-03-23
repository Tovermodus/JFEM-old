

import java.util.ArrayList;

public class TPCell extends Cell
{
	Cell1D cellx;
	Cell1D celly;
	int polynomialDegree;
	public TPCell(double xStart, double yStart, double xEnd, double yEnd, int polynomialDegree)
	{
		super();
		cellx = new Cell1D(xStart, xEnd);
		celly = new Cell1D(yStart, yEnd);
		this.polynomialDegree = polynomialDegree;
	}
	public TPCell(Cell1D cellx,Cell1D celly, int polynomialDegree)
	{
		this.cellx = cellx;
		this.celly = celly;
		this.polynomialDegree = polynomialDegree;
	}

	@Override
	void addFace(Face face)
	{
		this.faces.add(face);
		face.cells.add(this);
	}

	@Override
	void distributeFunctions(ArrayList<ScalarShapeFunction> globalshapeFunctions)
	{
		for(int i = 0; i <= polynomialDegree; i++)
		{
			for(int j = 0; j <= polynomialDegree; j++)
				shapeFunctions.add(new TPShapeFunction(new LagrangeBasisFunction1D(polynomialDegree,
					i, cellx), new LagrangeBasisFunction1D(polynomialDegree,j,celly),this));
		}
		for(ScalarShapeFunction shapeFunction:shapeFunctions)
		{
			for(Face face:faces)
			{
				face.shapeFunctions.add(shapeFunction);
			}
			globalshapeFunctions.add(shapeFunction);
			shapeFunction.setGlobalIndex(globalshapeFunctions.size()-1);

		}
	}


	@Override
	boolean isInCell(DoubleTensor pos)
	{
		return pos.x()>=cellx.xStart && pos.x()<= cellx.xEnd && pos.y()>=celly.xStart && pos.y()<=celly.xEnd ;
	}

	@Override
	DoubleTensor center()
	{
		return DoubleTensor.vectorFromValues(0.5*(cellx.xStart+cellx.xEnd),0.5*(celly.xStart+celly.xEnd));
	}


}



