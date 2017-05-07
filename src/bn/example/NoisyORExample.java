package bn.example;

import java.util.ArrayList;

import bn.BNet;
import bn.BNode;
import bn.node.CPT;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Enumerable;
import bn.JPT;
import bn.node.NoisyOR;
import bn.Predef;
import bn.node.SmartNoisyOR;
import bn.alg.CGVarElim;
import bn.alg.EM;
import bn.alg.LearningAlg;
import bn.alg.Query;
import bn.alg.QueryResult;

public class NoisyORExample {
	
	public static void main1 (String [] args) {
		
		EnumVariable v_gene1 = Predef.Nominal(new String [] {"low",  "medium", "high"}, "Gene1");
		EnumVariable v_gene2 = Predef.Nominal(new String [] {"low",  "medium", "high"}, "Gene2");
		EnumVariable v_gene3 = Predef.Nominal(new String [] {"low",  "medium", "high"}, "Gene3");
		EnumVariable v_gene4 = Predef.Nominal(new String [] {"low",  "medium", "high"}, "Gene4");
		EnumVariable v_disease = Predef.Nominal(new String [] {"negative", "some symptoms", "terminal"}, "Disease");
		CPT gene1 = new CPT(v_gene1);
		CPT gene2 = new CPT(v_gene2);
		CPT gene3 = new CPT(v_gene3);
		CPT gene4 = new CPT(v_gene4);
		NoisyOR disease = new NoisyOR(v_disease, new EnumVariable[] {v_gene1, v_gene2, v_gene3, v_gene4}, new String [] {"high", "high", "high", "low"});
		//CPT disease = new CPT(v_disease, new EnumVariable[] {v_gene1, v_gene2, v_gene3, v_gene4});
		disease.put(new Object[] {"high", "low", "low", "medium"}, new EnumDistrib(v_disease.getDomain(), new double[]{0.1, 0.3, 0.6}));
		disease.put(new Object[] {"medium", "high", "medium", "medium"}, new EnumDistrib(v_disease.getDomain(), new double[]{0.8, 0.15, 0.05}));
		disease.put(new Object[] {"low", "high", "low", "medium"}, new EnumDistrib(v_disease.getDomain(), new double[]{0.75, 0.15, 0.1}));
		disease.put(new Object[] {"medium", "medium", "high", "medium"}, new EnumDistrib(v_disease.getDomain(), new double[]{0.8, 0.15, 0.05}));
		disease.put(new Object[] {"medium", "medium", "low", "low"}, new EnumDistrib(v_disease.getDomain(), new double[]{0.15, 0.15, 0.70}));
		disease.print();
		Enumerable en = disease.getVariable().getDomain();
		System.out.println(en.getIndex("terminal"));
		System.out.println(en.get(2));
		System.out.println(disease.get(new Object [] {"high", "medium", "low", "low"},"terminal"));
	}
	
	public static void main2 (String [] args) {
		EnumVariable v1 = Predef.Boolean();
        EnumVariable v2 = Predef.Boolean();
        EnumVariable v3 = Predef.Boolean();
        ArrayList<Object> thisList = new ArrayList<Object>();

        NoisyOR NoisyOR1 = new NoisyOR(v1, new EnumVariable[]{v2, v3}, new Object [] {true, true});
        NoisyOR1.put(new Object[]{true, false}, new EnumDistrib(v1.getDomain(), new double[]{0.78, 0.22}));
        NoisyOR1.put(new Object[]{false, true}, new EnumDistrib(v2.getDomain(), new double[]{0.6, 0.4}));
        NoisyOR1.put(new Object[]{false, false}, new EnumDistrib(v2.getDomain(), new double[]{0.15, 0.85}));
        NoisyOR1.print();
        int [] indicies = NoisyOR1.getTable().getIndices();
        for (int i : indicies) {
        	System.out.println(i);
        }
        Object [] key = NoisyOR1.getTable().getKey(1);
        for (Object k : key) {
        	System.out.println(k);
        }
        System.out.println(NoisyOR1.get(new Object[] {true, true}, false));
        //training example
        Object[][] values={
				// A      B      C      D  
				{true,  false, true,  true},
				{true,  false,  true,  true},
				{true,  false, true,  true},
				{true, false,  true,  false},
				// true false true -> poscount=3, negcount=1
				{false,  false,  true,  false},
				{false, false, true, false},
				{false, false, true, true},
				// true true true -> poscount=1, negcount=2
				{true,  false, false, true},
				{true, false,  false, true},
				{true,  false,  false, true},
				{true,  false,  false, true},
				{true,  false,  false, true},
				{true,  false,  false, false},
				// true false false ->poscount=4, negcount=1
				{false, true,  false, false},
				{false, true, false, false},
				{false, true, false, false},
				{false, true, false, false},
				{false, true, false, true},
				// false, true, false -> poscount = 1, negcount = 4
				{false, false, false, false},
				{false, false, false, false},
				{false, false, false, false},
				{false, false, false, true}};
				// false false false -> poscount = 1, negcount = 4
			// construct a BN with two root nodes
        EnumVariable va = Predef.Boolean("A");
        EnumVariable vb = Predef.Boolean("B");
        EnumVariable vc = Predef.Boolean("C");
        EnumVariable vd = Predef.Boolean("D");
        CPT a = new CPT(va);
        CPT b = new CPT(vb);
        CPT c = new CPT(vc);
        NoisyOR d = new NoisyOR(vd, new EnumVariable [] {va, vb, vc}, new Object [] {true, true, true});
        BNet bn = new BNet();
        bn.add(a,b,c,d);
        LearningAlg em = new EM(bn);
        ArrayList<BNode> nodes = new ArrayList<BNode>();
        nodes.add(a);nodes.add(b);nodes.add(c);nodes.add(d);
        em.train(values, nodes);
        d.print();
        //model trained - now to try inference
        CGVarElim ve = new CGVarElim();
        a.setInstance(true);
        b.setInstance(true);
        c.setInstance(true);
        ve.instantiate(bn);
        QueryResult qr = ve.infer(d);
        //QueryResult qr = inf.infer(q);
        JPT jpt=qr.getJPT();
		jpt.display();
		//d.print();
	}
	
	public static void main (String [] args) {
		Object[][] values={
				// A      B      C      D    E 
				{true,  false, true,  true},
				{true,  false,  true,  true},
				{true,  false, true,  true},
				{true, false,  true,  false},
				// true false true -> poscount=3, negcount=1
				{true,  true,  true,  false},
				{true, true, true, false},
				{true, true, true, true},
				// true true true -> poscount=1, negcount=2
				{true,  false, false, true},
				{true, false,  false, true},
				{true,  false,  false, true},
				// true false false ->poscount=3, negcount=0
				{false, true,  false, false},
				{false, true, false, false},
				{false, true, false, false},
				{false, true, false, false},
				{false, true, false, true},
				// false, true, false -> poscount = 1, negcount = 4
				{false, false, false, false},
				{false, false, false, false},
				{false, false, false, false},
				{false, false, false, true}};
		EnumVariable va = Predef.Boolean("A");
        EnumVariable vb = Predef.Boolean("B");
        EnumVariable vc = Predef.Boolean("C");
        EnumVariable vd = Predef.Boolean("D");
        CPT a = new CPT(va);
        CPT b = new CPT(vb);
        CPT c = new CPT(vc);
        SmartNoisyOR d = new SmartNoisyOR(vd, new EnumVariable [] {va, vb, vc}, new Object [] {true, true, true});
        BNet bn = new BNet();
        bn.add(a,b,c,d);
        LearningAlg em = new EM(bn);
        ArrayList<BNode> nodes = new ArrayList<BNode>();
        nodes.add(a);nodes.add(b);nodes.add(c);nodes.add(d);
        em.train(values, nodes);
        a.print();
        b.print();
        c.print();
        d.print();
		
	}

}
