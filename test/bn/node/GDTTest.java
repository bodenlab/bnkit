package bn.node;

import bn.BNet;
import bn.BNode;
import bn.Predef;
import bn.alg.EM;
import bn.prob.EnumDistrib;
import bn.prob.GaussianDistrib;
import dat.Continuous;
import dat.EnumVariable;
import dat.Variable;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

class GDTTest {
    Variable x3 = Predef.Real("X3");
    EnumVariable x2 = Predef.Boolean("X2");
    EnumVariable x1 = Predef.Nominal(new String[] {"A","B", "C"}, "X1");

    @Test
    void fromJSON() {
        GDT gdt = new GDT(x3, x2, x1);
        gdt.put(new GaussianDistrib(2.1, 0.5), true, "A");
        gdt.put(new GaussianDistrib(1.3, 0.5), false, "C");
        gdt.put(new GaussianDistrib(6.4, 0.5), false, "B");
        gdt.put(new GaussianDistrib(7.1, 0.5), false, "A");
        System.out.println(gdt.getStateAsText());
        System.out.println(gdt.toJSON());
        GDT gdt2 = GDT.fromJSON(gdt.toJSON(), x3, x2, x1);
        assertEquals(gdt2.getStateAsText(), gdt.getStateAsText());
        assertEquals(gdt2.toJSON().toString(), gdt.toJSON().toString());

    }

    @Test
    void testTrainGMM() {
        EnumVariable X = Predef.Boolean("X");
        Variable<Continuous> Y = Predef.Real("Y");
        GDT gdt = new GDT(Y, X);
        gdt.randomize(0);
        CPT cpt = new CPT(X);
        cpt.randomize(0);
        BNet bn = new BNet();
        bn.add(new BNode[] {gdt, cpt});
        EM em = new EM(bn);

        // Sample data (you should replace this with your actual data)
        double[][] data = {
                {2.0},
                {3.5},
                {1.5},
                {7.0},
                {8.0},
                {2.2},
                {7.7}
        };

        cpt.getDistrib().set(new double[] {0.5, 0.5});
        for (int i = 0; i < 10; i ++) {
            for (int d = 0; d < data.length; d++) {
                for (Object component : cpt.getDistrib().getDomain().getValues()) {
                    double weight = cpt.getDistrib().get(component);
                    gdt.countInstance(new Object[] {component}, data[d][0], weight);
                }
                gdt.maximizeInstance();
            }
        }
    }
}