public class AngularTPFace extends TPFace
{
	public AngularTPFace(Cell1D cell1d, double otherCoordinate, int normaldirection)
	{
		super(cell1d, otherCoordinate, normaldirection);
	}
	public AngularTPFace(TPFace face)
	{
		super(face.cell1d, face.otherCoordinate, face.normaldirection);
	}
}
