import java.util.ArrayList;

public abstract class Cell
{
	ArrayList<ScalarShapeFunction> shapeFunctions;
	ArrayList<Face> faces;
	boolean refined;
	public Cell()
	{
		faces = new ArrayList<>();
		shapeFunctions = new ArrayList<>();
	}
	abstract void addFace(Face face);
	abstract void distributeFunctions(ArrayList<ScalarShapeFunction> shapeFunctions);
	abstract boolean isInCell(DoubleTensor pos);
	abstract DoubleTensor center();
	abstract ArrayList<Cell> refine(ArrayList<Face> refinedFaces);

}
