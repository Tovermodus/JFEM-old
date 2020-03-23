import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.GMRES;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import no.uib.cipr.matrix.sparse.LinkedSparseMatrix;
import org.ojalgo.matrix.store.MatrixStore;
import org.ujmp.core.DenseMatrix;
import org.ujmp.core.Matrix;
import org.ujmp.core.SparseMatrix;
import org.ojalgo.matrix.store.Primitive64Store;
import org.ojalgo.matrix.store.SparseStore;
import org.ojalgo.matrix.task.iterative.GaussSeidelSolver;
import java.util.ArrayList;
import java.util.DuplicateFormatFlagsException;
import java.util.Iterator;
import java.util.Optional;

import  no.uib.cipr.matrix.sparse.BiCGstab;
public class DoubleTensor
{
	Matrix mat;
	boolean isSparse;
	int dimension;
	private DoubleTensor()
	{}
	public DoubleTensor(int size)
	{
		dimension = 1;
		this.isSparse=false;
		mat = DenseMatrix.Factory.zeros(size,1);
	}
	public DoubleTensor(int sizeX, int sizeY, boolean isSparse)
	{
		dimension = 2;
		this.isSparse=isSparse;
		if(!isSparse)
			mat = DenseMatrix.Factory.zeros(sizeX, sizeY);
		else
			mat = SparseMatrix.Factory.zeros(sizeX,sizeY);
		this.isSparse = isSparse;
	}
	private SparseStore<Double> getOJALGOMatrix()
	{

		if(dimension != 2|| !isSparse)
			throw new UnsupportedOperationException();
		SparseStore<Double> m = SparseStore.PRIMITIVE64.make(getM(),getN());
		Iterator<long[]> coordIterator = mat.allCoordinates().iterator();
		Iterator<Object> valueIterator = mat.allValues().iterator();
		while(coordIterator.hasNext() && valueIterator.hasNext())
		{
			long[] coord = coordIterator.next();
			double val = mat.getAsDouble(coord);
			m.set((int)coord[0],(int)coord[1],val);
		}
		return m;
	}
	private LinkedSparseMatrix getMTJMatrix()
	{
		if(dimension != 2|| !isSparse)
			throw new UnsupportedOperationException();
		LinkedSparseMatrix m = new LinkedSparseMatrix(getM(),getN());
		Iterator<long[]> coordIterator = mat.allCoordinates().iterator();
		Iterator<Object> valueIterator = mat.allValues().iterator();
		while(coordIterator.hasNext() && valueIterator.hasNext())
		{
			long[] coord = coordIterator.next();
			double val = mat.getAsDouble(coord);
			m.set((int)coord[0],(int)coord[1],val);
		}
		return m;
	}
	private Primitive64Store getOJALGOVector()
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		Primitive64Store m = Primitive64Store.FACTORY.make(size(),1);
		for(int i = 0; i < size(); i++)
		{
			m.set(i,at(i));
		}
		return m;
	}
	private DenseVector getMTJvector()
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		DenseVector m = new DenseVector(size());
		for(int i = 0; i < size(); i++)
		{
			m.set(i,at(i));
		}
		return m;
	}
	private DoubleTensor fromOJALGOVector(Primitive64Store v)
	{
		DoubleTensor m = new DoubleTensor(v.size());
		for(int i = 0; i < v.size(); i++)
		{
			m.set(i,v.get(i));
		}
		return m;
	}
	private DoubleTensor fromMTJvector(DenseVector v)
	{
		DoubleTensor m = new DoubleTensor(v.size());
		for(int i = 0; i < v.size(); i++)
		{
			m.set(i,v.get(i));
		}
		return m;
	}
	public DoubleTensor(DoubleTensor t)
	{
		mat = t.mat.clone();
		isSparse = t.isSparse;
		dimension = t.dimension;
	}
	public int size()
	{
		return (int)(mat.getColumnCount()*mat.getRowCount());
	}
	public int getM()
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		return (int) mat.getRowCount();
	}
	public int getN()
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		return (int) mat.getColumnCount();
	}
	public double at(int vectorIndex)
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		return mat.getAsDouble(vectorIndex,0);
	}
	public double at(int matrixIndex1, int matrixIndex2)
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		return mat.getAsDouble(matrixIndex1,matrixIndex2);
	}
	public void set(int vectorIndex, double value)
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		mat.setAsDouble(value,vectorIndex,0);
	}
	public void set(int matrixIndex1,int matrixIndex2, double value)
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		mat.setAsDouble(value,matrixIndex1,matrixIndex2);
	}
	public void add(int matrixIndex1,int matrixIndex2, double value)
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		set(matrixIndex1,matrixIndex2,at(matrixIndex1,matrixIndex2)+value);
	}
	public void add(int vectorIndex, double value)
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		set(vectorIndex,at(vectorIndex)+value);
	}
	public void print_formatted()
	{
		if(dimension == 2)
			for(int i = 0; i < getM(); i++)
			{
				System.out.print("[");
				for(int j = 0; j < getN(); j++)
				{
					if (Math.abs(at(i, j)) <= 1e-14)
						System.out.print("<>........  ");
					else
					{
						if(at(i,j)>=0)
							System.out.print("+");
						System.out.print(String.format("%6.3e", at(i, j)) + "  ");
					}
				}
				System.out.println("]");
			}
		if(dimension == 1)
		{
			System.out.print("[");
			for (int j = 0; j < size(); j++)
			{
				if (Math.abs(at(j)) <= 1e-14)
					System.out.print("<>........  ");
				else
				{
					if(at(j)>=0)
						System.out.print("+");
					System.out.print(String.format("%6.3e", at(j)) + "  ");
				}
			}
			System.out.println("]");
		}
	}
	public void print_log()
	{
		if(dimension != 2)
			return;
		for(int i = 0; i < getM(); i++)
		{
			for(int j = 0; j < getN(); j++)
			{
				if (Math.abs(at(i, j)) <= 1e-14)
					System.out.print("   \t");
				else
					System.out.print((int) Math.log10(Math.abs(at(j, i))) + "\t");
			}
			System.out.println();
		}
	}
	public static DoubleTensor vectorFromValues(double... values)
	{

		DoubleTensor ret = new DoubleTensor(values.length);
		for(int i = 0; i < values.length; i++)
			ret.set(i,values[i]);
		return ret;
	}
	public static DoubleTensor squareMatrixFromValues(double... values)
	{
		if(Math.pow((int)Math.sqrt(values.length),2) != values.length)
			throw new UnsupportedOperationException("nope");
		DoubleTensor ret = new DoubleTensor((int)Math.sqrt(values.length),(int)Math.sqrt(values.length),false);
		for(int i = 0; i < values.length; i++)
			ret.set(i,values[i]);
		return ret;
	}
	public double x()
	{
		if(dimension == 1)
			return at(0);
		else throw new UnsupportedOperationException("nope");
	}
	public double y()
	{
		if(dimension == 1)
			return at(1);
		else throw new UnsupportedOperationException("nope");
	}
	public double inner(DoubleTensor other)
	{       double ret = 0;
		if(this.size() == other.size())
		{
			for(int i = 0; i < this.size(); i++)
			{
				ret += this.at(i)*other.at(i);
			}
			return ret;
		}
		throw new UnsupportedOperationException();
	}
	public DoubleTensor add(DoubleTensor other)
	{
		if(this.dimension == other.dimension && this.size() == other.size())
		{
			DoubleTensor ret =  new DoubleTensor();
			ret.dimension = this.dimension;
			ret.isSparse = this.isSparse && other.isSparse;
			ret.mat = this.mat.plus(other.mat);
			return ret;
		}
		throw new UnsupportedOperationException();
	}
	public DoubleTensor mul(double scalar)
	{
		DoubleTensor ret = new DoubleTensor();
		ret.mat = this.mat.times(scalar);
		ret.dimension = dimension;
		ret.isSparse = isSparse;
		return ret;
	}
	public DoubleTensor mvmul(DoubleTensor vector)
	{
		if(vector.dimension != 1 || this.dimension != 2)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor();
		ret.isSparse = false;
		ret.dimension = 1;
		ret.mat = this.mat.mtimes(vector.mat);
		return ret;
	}

	public DoubleTensor mmmul(DoubleTensor other)
	{
		if(other.dimension != 2 || this.dimension != 2)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor();
		ret.dimension = 2;
		ret.mat = this.mat.mtimes(other.mat);
		ret.isSparse = ret.mat.isSparse();
		return ret;
	}
	public DoubleTensor transpose()
	{
		if(this.dimension != 2)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor();
		ret.dimension = 2;
		ret.mat = this.mat.transpose();
		ret.isSparse = ret.mat.isSparse();
		return ret;
	}
	public double vectorNorm()
	{
		if(this.dimension != 1)
			throw new UnsupportedOperationException();
		return Math.sqrt(this.inner(this));
	}
	public DoubleTensor sub(DoubleTensor other)
	{
		if(this.dimension == other.dimension && this.size() == other.size())
		{
			DoubleTensor ret =  new DoubleTensor();
			ret.dimension = this.dimension;
			ret.isSparse = this.isSparse && other.isSparse;
			ret.mat = this.mat.minus(other.mat);
			return ret;
		}
		throw new UnsupportedOperationException();
	}
	public DoubleTensor solveSymm(DoubleTensor rightHandSide)
	{
		if(rightHandSide.dimension!=1)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor();
		ret.mat = mat.solveSymm(rightHandSide.mat);
		ret.dimension = 1;
		ret.isSparse = false;
		return ret;
	}
	public DoubleTensor solve(DoubleTensor rightHandSide)
	{
		if(rightHandSide.dimension!=1)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor();
		ret.mat = mat.solve(rightHandSide.mat);
		ret.dimension = 1;
		ret.isSparse = false;
		return ret;
	}
	public DoubleTensor solveCG(DoubleTensor rightHandSide, double residualTolerance)
	{
		System.out.println(this.mat.isSparse());

		if(getN() != getN() || dimension != 2 || rightHandSide.dimension != 1)
			throw new UnsupportedOperationException();
		DoubleTensor z;
		DoubleTensor newResiduum;
		DoubleTensor iterate = new DoubleTensor(getM());
		DoubleTensor residuum = rightHandSide.sub(this.mvmul(iterate));
		DoubleTensor defect = new DoubleTensor(residuum);
		double α = 0;
		double β = 0;
		for(int iter = 0; iter < getM() && residuum.vectorNorm() > residualTolerance; iter++)
		{
			z = this.mvmul(defect);
			α = residuum.inner(residuum)/defect.inner(z);
			iterate = iterate.add(defect.mul(α));
			newResiduum = residuum.sub(z.mul(α));
			β = newResiduum.inner(newResiduum)/residuum.inner(residuum);
			defect = newResiduum.add(defect.mul(β));
			residuum = newResiduum;
			if(iter%10 == 0)
				System.out.println("iter "+iter+" " +residuum.vectorNorm());
		}
		return iterate;
	}
	public DoubleTensor solveGMRES(DoubleTensor rightHandSide)
	{
		DenseVector rhs = rightHandSide.getMTJvector();
		DenseVector v = new DenseVector(getM());
		DenseVector x = new DenseVector(getM());
		LinkedSparseMatrix m = getMTJMatrix();
		BiCGstab bcg = new BiCGstab(v);
		GMRES gm = new GMRES(v);
		try
		{
			x = (DenseVector) gm.solve(m, rhs, x);
		} catch (IterativeSolverNotConvergedException e)
		{
			e.printStackTrace();
		}
		System.out.println(x.size());
		return fromMTJvector(x);
	}
	public DoubleTensor solveBiCGStab(DoubleTensor rightHandSide, double residualTolerance)
	{


		if(getN() != getN() || dimension != 2 || rightHandSide.dimension != 1)
			throw new UnsupportedOperationException();
		DoubleTensor v = new DoubleTensor(getM());
		DoubleTensor s = new DoubleTensor(getM());
		DoubleTensor t = new DoubleTensor(getM());
		DoubleTensor iterate = new DoubleTensor(getM());
		DoubleTensor startResiduum = rightHandSide.sub(this.mvmul(iterate));
		DoubleTensor residuum = new DoubleTensor(startResiduum);
		DoubleTensor p = new DoubleTensor(residuum);
		double α = 0;
		double ω = 0;
		double ρ = residuum.inner(residuum);
		double ρlast = residuum.inner(residuum);
		double β = 0;
		for(int iter = 0; iter < getM() && residuum.vectorNorm() > residualTolerance; iter++)
		{
			v = this.mvmul(p);
			α = ρ/v.inner(startResiduum);
			s = residuum.sub(v.mul(α));
			t = this.mvmul(s);
			ω = t.inner(s)/t.inner(t);
			iterate = iterate.add(p.mul(α)).add(s.mul(ω));
			residuum = s.sub(t.mul(ω));
			ρlast = ρ;
			ρ = residuum.inner(startResiduum);
			β = α/ω*ρ/ρlast;
			p = residuum.add(p.mul(β)).sub(v.mul(ω*β));
			if(iter%10 == 0)
			System.out.println("iter "+iter+" "+residuum.vectorNorm());
		}
		return iterate;


	}
}
