package bn.node;

import bn.Predef;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Enumerable;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class CPTTest {

    EnumVariable x3 = Predef.Number(3, "X3");
    EnumVariable x2 = Predef.Boolean("X2");
    EnumVariable x1 = Predef.Nominal(new String[] {"A","B", "C"}, "X1");

    @Test
    void fromJSON() {
        CPT cpt = new CPT(x3, x2, x1);
        cpt.put(EnumDistrib.random(x3.getDomain()), true, "A");
        cpt.put(EnumDistrib.random(x3.getDomain()), false, "C");
        cpt.put(EnumDistrib.random(x3.getDomain()), false, "B");
        cpt.put(EnumDistrib.random(x3.getDomain()), false, "A");
        System.out.println(cpt.getStateAsText());
        System.out.println(cpt.toJSON());
        CPT cpt2 = CPT.fromJSON(cpt.toJSON(), x3, x2, x1);
        assertEquals(cpt2.getStateAsText(), cpt.getStateAsText());
        assertEquals(cpt2.toJSON().toString(), cpt.toJSON().toString());
    }
}