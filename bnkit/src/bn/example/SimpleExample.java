/**
 * 
 */
package bn.example;

import bn.BNet;
import bn.node.CPT;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Enumerable;
import bn.Predef;
import dat.Variable.Assignment;
import bn.alg.ApproxInference;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
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
            System.out.println("Query is \'" + q + "\'");
            CGTable r0 = (CGTable)ve.infer(q);
            r0.display();

            j.setInstance(false);
            q = ve.makeQuery(B);
            System.out.println("Query is \'" + q + "\'");
            r0 = (CGTable)ve.infer(q);
            r0.display();
            
            j.resetInstance();
            m.resetInstance();
            a.setInstance(true);
            q = ve.makeQuery(J,M,E);
            System.out.println("Query is \'" + q + "\'");
            r0 = (CGTable)ve.infer(q);
            r0.display();

            Query q_mpe = ve.makeMPE();
            System.out.println("Query is \'" + q_mpe + "\'");
            CGTable r1 = (CGTable)ve.infer(q_mpe); 
            Assignment[] assign = r1.getMPE();
            for (Assignment assign1 : assign) {
                System.out.println("\t" + assign1.var + " = " + assign1.val);
            }
            r1.display();
            
            b.resetInstance();
            j.setInstance(true);
            m.setInstance(true);
            System.out.println("Approximate Inference------------");
            ApproxInference ai = new ApproxInference();
            ai.instantiate(bn);
            ai.setIterations(1000);
            Query qq = ai.makeQuery(B);
            CGTable cgt = (CGTable)ai.infer(qq);
            cgt.display();

            //BNBuf.save(bn, "data/bn_simple.xml");
		
	}
	
}
