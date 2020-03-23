public class TPBoundaryFaceIntegral extends BoundaryFaceIntegral
{
	public TPBoundaryFaceIntegral(ScalarFunction rightHandSide, ScalarFunction weight)
	{
		super(rightHandSide, weight);
	}

	@Override
	public double evaluateBoundaryFaceIntegral(Face f, ScalarShapeFunction v)
	{
		if(f instanceof TPFace)
			return this.evaluateBoundaryFaceIntegral((TPFace) f, v);
		throw new UnsupportedOperationException();
	}

	public double evaluateBoundaryFaceIntegral(TPFace f, ScalarShapeFunction v)
	{
		if(!f.isBoundaryFace)
			return 0;
		double ret = 0;
		for(int i = 0; i < f.cell1d.weights.length; i++)
		{
			double weight = f.cell1d.weights[i];
			DoubleTensor point = new DoubleTensor(2);
			point.set(f.normaldirection,f.cell1d.points[i]);
			point.set(1-f.normaldirection,f.otherCoordinate);
			ret += v.value(point)*rightHandSide.value(point)*this.weight.value(point)*weight;
		}
		return ret;
	}
}
