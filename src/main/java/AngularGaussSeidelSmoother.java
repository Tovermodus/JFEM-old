import no.uib.cipr.matrix.LowerTriangPackMatrix;
import org.ujmp.core.Matrix;

public class AngularGaussSeidelSmoother extends AngularSmoother
{
	DoubleTensor r;
	DoubleTensor l;
	DoubleTensor linv;
	LowerTriangPackMatrix low;
	public AngularGaussSeidelSmoother(AngularGrid g, String[] args)
	{
		super(g, args);
		g.A.print_formatted();
		l = new DoubleTensor(g.A.getM(),g.A.getN(),g.A.isSparse);
		r = new DoubleTensor(g.A.getM(),g.A.getN(),g.A.isSparse);
		if(g.A.isSparse)
		{
			for(int i = 0; i < g.A.sparseEntries; i++)
			{
				if(g.A.sparseYs[i] >= g.A.sparseXs[i])
				{
					l.add(g.A.sparseYs[i],g.A.sparseXs[i],g.A.sparseValues[i]);
				}
				else
					r.add(g.A.sparseYs[i],g.A.sparseXs[i],g.A.sparseValues[i]);
			}

		}
		//linv = l.inverse();
		low = MatrixFormats.getLowerTriangMTJMatrix(l);
		//l.print_formatted();
		//MatrixFormats.printMTJMatrix(low);
	}

	@Override
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		//System.out.println("kj");
		//iterate = l.solveBiCGStab(rightHandSide.sub(r.mvmul(iterate)),1e-11);//iterate.fromUJMPVector(lm.solve
		//iterate = linv.mvmul(rightHandSide.sub(r.mvmul(iterate)));


		iterate =
			MatrixFormats.fromMTJvector(low.solve(MatrixFormats.getMTJvector(rightHandSide.sub(r.mvmul
			(iterate))),
				MatrixFormats.getMTJvector(iterate)));

		// (rightHandSide
		// .sub(r
		// .mvmul
		// (iterate))
		// .getUJMPVector()));
		//System.out.println("sad");
		return iterate;
	}
}
