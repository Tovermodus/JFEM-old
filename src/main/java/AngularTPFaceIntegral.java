public class AngularTPFaceIntegral extends AngularFaceIntegral
{
	protected AngularTPFaceIntegral()
	{
	}

	public AngularTPFaceIntegral(AngularScalarFunction weight, String name)
	{
		super(weight,name);
	}


	public AngularTPFaceIntegral(String name)
	{
		super(name);
	}

	@Override
	public double evaluateAngularFaceIntegral(Face face, ScalarShapeFunction function1, ScalarShapeFunction function2, DoubleTensor direction1, DoubleTensor direction2, double directionWeight1, double directionWeight2)
	{

		if(face instanceof TPFace)
			return evaluateAngularFaceIntegral((TPFace) face,function1,function2,direction1, direction2,
				directionWeight1,directionWeight2);
		throw new UnsupportedOperationException();
	}

	public double evaluateAngularFaceIntegral(TPFace face,
	                                          ScalarShapeFunction function1, ScalarShapeFunction function2,
	                                          DoubleTensor direction1, DoubleTensor direction2,
	                                          double directionWeight1, double directionWeight2)
	{
		if(name.equals(VALUE_VALUE))
		{
			if(direction1.sub(direction2).vectorNorm()>=1e-14|| !face.isBoundaryFace || ((TPFace)face).getDownStreamCell(face.center(),direction1) == null)
				return 0.;
			double ret = 0;
			for(int i = 0; i < face.cell1d.weights.length; i++)
			{
				double weight = face.cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-face.normaldirection,face.cell1d.points[i]);
				point.set(face.normaldirection,face.otherCoordinate);
				ret += function1.value(point)*function2.value(point)
					*Math.abs(direction1.inner(face.normal.value(point)))
					*weight*directionWeight1*this.weight.value(point,direction1);
			}
			return ret;
		}
		if(name.equals(JUMP_NORMALAVERAGE_JUMP_NORMALAVERAGE))
		{

			if(direction1.sub(direction2).vectorNorm()>=1e-14||face.isBoundaryFace)
				return 0.;
			double ret = 0;
			for(int i = 0; i < face.cell1d.weights.length; i++)
			{
				double integralWeight = face.cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-face.normaldirection,face.cell1d.points[i]);
				point.set(face.normaldirection,face.otherCoordinate);
				ret += face.normalAverageInValue(function1,point).inner(face.normalAverageInValue(function2,point))
					*Math.abs(direction1.inner(face.normal.value(point)))
					*integralWeight*directionWeight1*this.weight.value(point,direction1);
			}
			return ret;
		}
		if(name.equals(JUMP_NORMALAVERAGE_JUMP))
		{

			if(direction1.sub(direction2).vectorNorm()>=1e-14||face.isBoundaryFace)
				return 0.;
			double ret = 0;
			for(int i = 0; i < face.cell1d.weights.length; i++)
			{
				double weight = face.cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-face.normaldirection,face.cell1d.points[i]);
				point.set(face.normaldirection,face.otherCoordinate);
				ret += -2*face.normalAverageInValue(function1,point).inner(direction1)
					*face.jumpInValue(function2,point)*weight*directionWeight1*this.weight.value(point,direction1);
			}
			return ret;
		}
		if(name.equals(JUMPJUMP))
		{
			if(direction1.sub(direction2).vectorNorm()>=1e-14||(face.isBoundaryFace))
				return 0.;
			double ret = 0;
			for (int i = 0; i < face.cell1d.weights.length; i++)
			{
				double weight = face.cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-face.normaldirection, face.cell1d.points[i]);
				point.set(face.normaldirection, face.otherCoordinate);
				//if(((TPFace)face).getDownStreamCell(point,direction1) == null)
				//	return 0;
				ret += face.jumpInValue(function1, point) * face.jumpInValue(function2, point)
					* Math.abs(direction1.inner(face.normal.value(point)))
					* this.weight.value(point,direction1) * weight * directionWeight1;
			}
			return ret;
		}
		if(name.equals(UPWIND))
		{

			if(direction1.sub(direction2).vectorNorm()>=1e-14|| face.getDownStreamCell(face.center(),
				direction1) == null)
				return 0.;
			double ret = 0;
			for(int i = 0; i < face.cell1d.weights.length; i++)
			{
				double weight = face.cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-face.normaldirection, face.cell1d.points[i]);
				point.set(face.normaldirection, face.otherCoordinate);
				//if(((TPFace)face).getDownStreamCell(point,direction1) == null)
				//	return 0;
				//direction1.print_formatted();
				//face.normal.value(point).print_formatted();
				ret += face.jumpInValue(function1, point) * face.jumpInValue(function2, point) *
				this.weight.value(point,
					direction1) * weight * Math.abs(direction1.inner(face.normal.value(point)));
//				ret += (face.averageInValue(function1, point)*direction1.inner(face.normal.value(point))+
//					0.5*face.jumpInValue(function1,point)*Math.abs(direction1.inner(face.normal.value(point))))*function2.value(point)
//				*weight*directionWeight1;
			}
			return ret;
		}
		throw new UnsupportedOperationException();
	}
}
