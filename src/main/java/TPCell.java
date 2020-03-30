

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
		return pos.x()>=cellx.start && pos.x()<= cellx.end && pos.y()>=celly.start && pos.y()<=celly.end;
	}

	@Override
	DoubleTensor center()
	{
		return DoubleTensor.vectorFromValues(0.5*(cellx.start +cellx.end),0.5*(celly.start +celly.end));
	}

	@Override
	ArrayList<Cell> refine(ArrayList<Face> refinedFaces)
	{
		ArrayList<Cell> refinedCells = new ArrayList<>();
		TPCell cell1 = new TPCell(cellx.start,celly.start,cellx.center(),celly.center(),polynomialDegree);
		TPFace face1 = new TPFace(new Cell1D(cellx.start,cellx.center()),celly.center(),1);
		TPFace face2 = new TPFace(new Cell1D(celly.start,celly.center()),cellx.center(),0);
		TPCell cell2 = new TPCell(cellx.center(),celly.start,cellx.end,celly.center(),polynomialDegree);
		TPCell cell3 = new TPCell(cellx.start,celly.center(),cellx.center(),celly.end,polynomialDegree);
		TPCell cell4 = new TPCell(cellx.center(),celly.center(),cellx.end,celly.end,polynomialDegree);
		TPFace face3 = new TPFace(new Cell1D(cellx.center(),cellx.end),celly.center(),1);
		TPFace face4 = new TPFace(new Cell1D(celly.center(),celly.end),cellx.center(),0);
		cell1.faces.add(face1);
		cell1.faces.add(face2);
		cell2.faces.add(face2);
		cell2.faces.add(face3);
		cell3.faces.add(face1);
		cell3.faces.add(face4);
		cell4.faces.add(face3);
		cell4.faces.add(face4);
		face1.cells.add(cell1);
		face1.cells.add(cell3);
		face2.cells.add(cell1);
		face2.cells.add(cell2);
		face3.cells.add(cell2);
		face3.cells.add(cell4);
		face4.cells.add(cell3);
		face4.cells.add(cell4);
		face1.isBoundaryFace = false;
		face2.isBoundaryFace = false;
		face3.isBoundaryFace = false;
		face4.isBoundaryFace = false;
		refinedCells.add(cell1);
		refinedCells.add(cell2);
		refinedCells.add(cell3);
		refinedCells.add(cell4);
		refinedFaces.add(face1);
		refinedFaces.add(face2);
		refinedFaces.add(face3);
		refinedFaces.add(face4);
		refined = true;
		return refinedCells;
	}
	public void print()
	{
		System.out.println("["+cellx.start+","+cellx.end+"]×["+celly.start+","+celly.end+"]");
	}


}



