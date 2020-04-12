package bn.example;

import bn.BNet;
import bn.JPT;
import bn.Predef;
import bn.alg.*;
import bn.file.BNBuf;
import bn.node.CPT;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Enumerable;
import dat.Variable;

/*
    Example application around Brad's question about causal queries:
    Intelligence = I, Grade = G, Conversation (a.k.a interview) = C, Job (get or not) = J.
 */
public class CausalQuery {

    public static void main(String[] args) {

        // define a few enumerable domains...
        Enumerable ilevels = new Enumerable(new String[] {"Low", "Medium", "High"});
        Enumerable grades = new Enumerable(new String[] {"A", "B", "C", "Fail"});
        Enumerable intervs = new Enumerable(new String[] {"Impress", "Good", "Flunk"});

        // define variables, all "enumerable" using the above domains, + predefined Boolean
        EnumVariable INTEL = new EnumVariable(ilevels, "Intel");
        EnumVariable GRADE = new EnumVariable(grades, "Grade");
        EnumVariable IVIEW = new EnumVariable(intervs, "Interv");
        EnumVariable GOTJOB = Predef.Boolean("GotJob");

        // give the "arrows" between, and parameters of all BN nodes, i.e. define conditional probability tables
        CPT intel = new CPT(INTEL);                                         // CPT for P(INTEL)
        intel.put(new EnumDistrib(ilevels, 0.333, 0.333, 0.333));  // distribution of levels of intelligence (only if not set)
        CPT grade = new CPT(GRADE, INTEL);                                  // CPT for P(GRADE|INTEL)
        grade.put(new EnumDistrib(grades, 0.05, 0.20, 0.30, 0.45), "Low");
        grade.put(new EnumDistrib(grades, 0.20, 0.30, 0.30, 0.20), "Medium");
        grade.put(new EnumDistrib(grades, 0.35, 0.40, 0.20, 0.05), "High");
        CPT iview = new CPT(IVIEW, INTEL);	                                // CPT for P(IVIEW|INTEL)
        iview.put(new EnumDistrib(intervs, 0.10, 0.30, 0.60), "Low");
        iview.put(new EnumDistrib(intervs, 0.25, 0.50, 0.25), "Medium");
        iview.put(new EnumDistrib(intervs, 0.60, 0.30, 0.10), "High");
        CPT gotjob = new CPT(GOTJOB, GRADE, IVIEW);                         // CPT for P(GOTJOB|GRADE,IVIEW)
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.99, 0.01), "A", "Impress");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.75, 0.25), "B", "Impress");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.50, 0.50), "C", "Impress");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.25, 0.75), "Fail", "Impress");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.90, 0.10), "A", "Good");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.65, 0.35), "B", "Good");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.30, 0.70), "C", "Good");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.10, 0.90), "Fail", "Good");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.50, 0.50), "A", "Flunk");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.25, 0.75), "B", "Flunk");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.10, 0.90), "C", "Flunk");
        gotjob.put(new EnumDistrib(Enumerable.bool, 0.01, 0.99), "Fail", "Flunk");

        // define the scope of the Bayesian network
        BNet bn=new BNet();
        bn.add(intel, grade, iview, gotjob);

        // print out all tables as specified above
        intel.print();
        grade.print();
        iview.print();
        gotjob.print();

        // set up inference
		VarElim inf = new VarElim();
        // instantiate GRADE to "C"
        grade.setInstance("C");
        // must let BN know to prepare for query/queries
        inf.instantiate(bn);

        // specify what variables that we want to query, i.e. query : P(GOTJOB, INTEL, IVIEW | variables instantiated above)
        Query q = inf.makeQuery(new Variable[] {GOTJOB, INTEL, IVIEW});
        QueryResult qr = inf.infer(q);
        JPT jpt=qr.getJPT();
        jpt.display(); // print out the results, that is a joint probability table with all possible value combinations

        // same query but MAP (aka MPE most probable explanation)
        q = inf.makeMPE(new Variable[] {GOTJOB, INTEL, IVIEW});
        CGTable cgt = (CGTable)inf.infer(q);
        Variable.Assignment[] as = cgt.getMPE();
        for (Variable.Assignment a : as)
            System.out.println("\t" + a);

        // P(INTEL|GRADE=C)
        q = inf.makeQuery(new Variable[] {INTEL});
        qr = inf.infer(q);
        jpt=qr.getJPT();
        jpt.display();

        // P(IVIEW|GRADE=C)
        q = inf.makeQuery(new Variable[] {IVIEW});
        qr = inf.infer(q);
        jpt=qr.getJPT();
        jpt.display();

        grade.setInstance("A");
        inf.instantiate(bn);

        // P(INTEL|GRADE=A)
        q = inf.makeQuery(new Variable[] {INTEL});
        qr = inf.infer(q);
        jpt=qr.getJPT();
        jpt.display();

        // P(IVIEW|GRADE=A)
        q = inf.makeQuery(new Variable[] {IVIEW});
        qr = inf.infer(q);
        jpt=qr.getJPT();
        jpt.display();

    }

}
