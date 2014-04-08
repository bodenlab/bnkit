/**
 * 
 */
package bn.example;

import bn.BNet;
import bn.CPT;
import bn.EnumDistrib;
import bn.EnumVariable;
import bn.Enumerable;
import bn.JPT;
import bn.Predef;
import bn.alg.ApproxInfer;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.Variable;
import bn.file.BNBuf;

/**
 * @author mikael
 */
public class SimpleExample {

	
	/**
	 * Example taken from Russell and Norvig (2003; p. 493-494)
	 * @param args
	 */
	public static void main(String[] args) {
		EnumVariable B = Predef.Boolean("Burglary");
		EnumVariable E = Predef.Boolean("Earthquake");
		EnumVariable A = Predef.Boolean("Alarm");
		EnumVariable J = Predef.Boolean("John calls");
		EnumVariable M = Predef.Boolean("Mary calls");
		
		CPT b = new CPT(B);
		CPT e = new CPT(E);
		CPT a = new CPT(A,    B, E);
		CPT j = new CPT(J,    A);
		CPT m = new CPT(M,    A);
		
		b.put(new EnumDistrib(Enumerable.bool, 0.001, 0.999));
		b.print();
		e.put(new EnumDistrib(Enumerable.bool, 0.002, 0.998));
		e.print();
		a.put(new EnumDistrib(Enumerable.bool, 0.95, 0.05), true, true);
		a.put(new EnumDistrib(Enumerable.bool, 0.94, 0.06), true, false);
		a.put(new EnumDistrib(Enumerable.bool, 0.29, 0.71), false, true);
		a.put(new EnumDistrib(Enumerable.bool, 0.001, 0.999), false, false);
		a.print();
		j.put(new EnumDistrib(Enumerable.bool, 0.90, 0.10), true);
		j.put(new EnumDistrib(Enumerable.bool, 0.05, 0.95), false);
		j.print();
		m.put(new EnumDistrib(Enumerable.bool, 0.70, 0.30), true);
		m.put(new EnumDistrib(Enumerable.bool, 0.01, 0.99), false);
		m.print();
		
		BNet bn = new BNet();
		bn.add(b, e, a, j, m);

		j.setInstance(true);
		m.setInstance(true);
		//e.setInstance(false);

		System.out.println("Variable elimination--------------");
		VarElim ve = new VarElim();
		ve.instantiate(bn);
		Query q = ve.makeQuery(B);
		JPT jpt = ve.infer(q);
		jpt.display();
		j.setInstance(false);
		jpt = ve.infer(q);
		jpt.display();
		
		j.setInstance(true);
		m.setInstance(true);
		System.out.println("Approximate Inference------------");
		ApproxInfer ai = new ApproxInfer();
		ai.instantiate(bn);
		Query qq = ai.makeQuery(B);
		JPT jpt1 = ai.infer(qq);
		jpt1.display();
		j.setInstance(false);
		jpt1 = ai.infer(qq);
		jpt1.display();
		
		
//		BNBuf.save(bn, "data/bn_simple.xml");
		
		/*

		
		*/
		
	}
	
}
